import glob

localrules: MetaWrap_refine_check

checkpoint MetaWrap_refine_check:
    input:
        expand("Binning/{binner}/bins",binner=BIN_STEPS )
    output:
        "status/MetaWrap_refine.ready"
    shell:
        "touch {output}"

mw_min_comp = config['binning']['MetaWrap']['min_completeness']
mw_max_cont = config['binning']['MetaWrap']['max_contamination']
binner_out_path_list = [OUTPUTDIR + "/Binning/" + binner + "/bins" for binner in BIN_STEPS
                        if (glob.glob(OUTPUTDIR + "/Binning/" + binner + "/bins/*.fasta")
                            and any(os.stat(bin).st_size > 0 for bin in glob.glob(OUTPUTDIR + "/Binning/" + binner +
                                                                                  "/bins/*.fasta")))]
if glob.glob(OUTPUTDIR + "/Binning/DASTool/bins/*.fasta"):
    binner_out_path_list.extend([OUTPUTDIR + "/Binning/DASTool/bins"])

rule MetaWrap_refine:
    input:
        binner_fasta_dirs=check_binner_bin_fastas,
        setup_check=f"{DBPATH}/checkm/checkm_db.ready"
    output:
        "Binning/MetaWrap/bins.stats",
        directory("Binning/MetaWrap/bins"),
        out_dir=directory("Binning/MetaWrap/")
    params:
        binner_inputs_list=get_mw_binner_inputs(binner_out_path_list, mw_min_comp, mw_max_cont),
        mw_min_comp=mw_min_comp,
        mw_max_cont=mw_max_cont,
        prefix_1="metawrap_round_01",
        prefix_2="metawrap_round_02",
        prefix_3="metawrap_round_03",
        total_mem=int(float(int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]) * getThreads(BIGCORENO))
    threads:
        getThreads(BIGCORENO)
    resources:
        runtime = "48:00:00",
        mem = BIGMEMCORE
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning_MetaWrap.log"
    message:
        "MetaWrap: Running MetaWrap."
    shell:
        """
        # Make sure conda env MetaWrap is used
        export PATH=${{CONDA_PREFIX}}/bin:${{PATH}}
        # Set up input array
        final_out={params.prefix_1}
        bin_folder="metawrap_{params.mw_min_comp}_{params.mw_max_cont}_bins"
        binner_input_array=({input.binner_fasta_dirs})
        
        # Ensure no special chars in bin names, which lead to MetaWrap CheckM crash
        for i in ${{binner_input_array[@]}}; do
          for bin in ${{i}}/*.fasta; do
            fixed_name=$(python -c "print(\\"${{bin}}\\".replace('|', '___'))")
            if [[ "${{fixed_name}}" != "${{bin}}" ]]; then
              mv "${{bin}}" "${{fixed_name}}"
            fi
          done
        done
        
        if [ ${{#binner_input_array[@]}} -eq 0 ]; then
          mkdir -p {output.out_dir}/bins
          touch {output.out_dir}/bins.stats {output.out_dir}/bins.contigs
          exit 0
        fi

        # Round 1
        metawrap bin_refinement -o {output.out_dir}/{params.prefix_1} \
                                --quick \
                                -t {threads} \
                                -m {params.total_mem} \
                                {params.binner_inputs_list[0]} \
                                -c {params.mw_min_comp} \
                                -x {params.mw_max_cont} >> {log} 2>&1
                
        # Rename bin fastas back to .fasta
        for bin in {output.out_dir}/{params.prefix_1}/${{bin_folder}}/*.fa; do
          bin_id=$(basename -s .fa ${{bin}})
          mv ${{bin}} ${{bin%/*}}/${{bin_id}}.fasta
        done
        rm -rf {output.out_dir}/{params.prefix_1}/${{bin_folder}}/*/
        
        # Round 2
        if [[ ${{#binner_input_array[*]}} > 3 ]]; then
            mv {output.out_dir}/{params.prefix_1}/${{bin_folder}} {output.out_dir}/{params.prefix_1}/${{bin_folder}}_round_01
            metawrap bin_refinement -o {output.out_dir}/{params.prefix_2} \
                    -t {threads} \
                    -m {params.total_mem} \
                    {params.binner_inputs_list[1]} \
                    -c {params.mw_min_comp} \
                    -x {params.mw_max_cont} >> {log} 2>&1
            
            # Rename bin fastas back to .fasta
            for bin in {output.out_dir}/{params.prefix_2}/${{bin_folder}}/*.fa; do
              bin_id=$(basename -s .fa ${{bin}})
              mv ${{bin}} ${{bin%/*}}/${{bin_id}}.fasta
            done
            rm -rf {output.out_dir}/{params.prefix_1} {output.out_dir}/{params.prefix_2}/${{bin_folder}}/*/
            final_out={params.prefix_2}
        fi

        # Round 3
        if [[ ${{#binner_input_array[*]}} > 5 ]]; then
            mv {output.out_dir}/{params.prefix_2}/${{bin_folder}} {output.out_dir}/{params.prefix_2}/${{bin_folder}}_round_02
            metawrap bin_refinement -o {output.out_dir}/{params.prefix_3} \
                    -t {threads} \
                    -m {params.total_mem} \
                    {params.binner_inputs_list[2]} \
                    -c {params.mw_min_comp} \
                    -x {params.mw_max_cont} >> {log} 2>&1
            
            # Rename bin fastas back to .fasta
            for bin in {output.out_dir}/{params.prefix_3}/${{bin_folder}}/*.fa; do
              bin_id=$(basename -s .fa ${{bin}})
              mv ${{bin}} ${{bin%/*}}/${{bin_id}}.fasta
            done
            rm -rf {output.out_dir}/{params.prefix_2}  {output.out_dir}/{params.prefix_3}/${{bin_folder}}/*/
            final_out={params.prefix_3}
        fi

        # Clean up
        rm -rf {output.out_dir}/${{final_out}}/bins \
               {output.out_dir}/${{final_out}}/work_files
        mv {output.out_dir}/${{final_out}}/* {output.out_dir}/ && rm -rf {output.out_dir}/${{final_out}}
        for i in {output.out_dir}/metawrap_{params.mw_min_comp}_{params.mw_max_cont}*; do
          mv ${{i}} ${{i/metawrap_{params.mw_min_comp}_{params.mw_max_cont}_/}}
        done
        

        if [[ ! -d Visualization ]]; then
          mkdir Visualization
        fi
        
        mv {output.out_dir}/figures/* Visualization/
        rm -rf {output.out_dir}/figures
        
        # Revert no special chars in bin names, which lead to MetaWrap CheckM crash
        for i in ${{binner_input_array[@]}}; do
          for bin in ${{i}}/*.fasta; do
            unfixed_name=$(python -c "print(\\"${{bin}}\\".replace('___', '|'))")
            if [[ "${{unfixed_name}}" != "${{bin}}" ]]; then
              mv "${{bin}}" "${{unfixed_name}}"
            fi
          done
        done
        """


