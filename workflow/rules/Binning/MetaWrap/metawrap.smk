import glob

localrules: MetaWrap_refine_check

checkpoint MetaWrap_refine_check:
    input:
        expand("Binning/{binner}/bins",binner=BIN_STEPS + ['DASTool']),
        expand("Binning/{binner}",binner=BIN_STEPS + ['DASTool'])
    output:
        "status/MetaWrap_refine.ready"
    shell:
        "touch {output}"

rule MetaWrap_refine:
    input:
        binner_fasta_dirs = check_binner_bin_fastas,
        setup_check = f"{DBPATH}/checkm/checkm_db.ready"
    output:
        "Binning/MetaWrap/bins.stats",
        directory("Binning/MetaWrap/bins"),
        out_dir=directory("Binning/MetaWrap/")
    params:
        binner_inputs_list=get_mw_binner_inputs,
        mw_min_comp=config['binning']['MetaWrap']['min_completeness'],
        mw_max_cont=config['binning']['MetaWrap']['max_contamination'],
        prefix_1="metawrap_round_01",
        prefix_2="metawrap_round_02",
        prefix_3="metawrap_round_03",
        total_mem=int(float(int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]) * getThreads(MAIN_BIGCORENO, bigmem=True)),
        script_path = os.path.join(SRCDIR, "MW_scripts"),
        checkm_db = os.path.join(DBPATH, "checkm")
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning/MetaWrap.log"
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
        CHECKM_DATA_PATH={params.checkm_db} SOFT={params.script_path} \
        {params.script_path}/bin_refinement.sh -o {output.out_dir}/{params.prefix_1} \
                                               --quick \
                                               -t {threads} \
                                               -m {params.total_mem} \
                                               {params.binner_inputs_list[0]} \
                                               -c {params.mw_min_comp} \
                                               -x {params.mw_max_cont} >> {log} 2>&1
                
        # Rename bin fastas back to .fasta
        shopt -s nullglob
        for bin in {output.out_dir}/{params.prefix_1}/${{bin_folder}}/*.fa; do
          bin_id=$(basename -s .fa ${{bin}})
          mv ${{bin}} ${{bin%/*}}/${{bin_id}}.fasta
        done
        shopt -u nullglob
        rm -rf {output.out_dir}/{params.prefix_1}/${{bin_folder}}/*/
        
        # Round 2
        if [[ ${{#binner_input_array[*]}} > 3 ]]; then
            mv {output.out_dir}/{params.prefix_1}/${{bin_folder}} {output.out_dir}/{params.prefix_1}/${{bin_folder}}_round_01
            CHECKM_DATA_PATH={params.checkm_db} SOFT={params.script_path} \
            {params.script_path}/bin_refinement.sh -o {output.out_dir}/{params.prefix_2} \
                                                   --quick \
                                                   -t {threads} \
                                                   -m {params.total_mem} \
                                                   {params.binner_inputs_list[1]} \
                                                   -c {params.mw_min_comp} \
                                                   -x {params.mw_max_cont} >> {log} 2>&1
            
            # Rename bin fastas back to .fasta
            shopt -s nullglob
            for bin in {output.out_dir}/{params.prefix_2}/${{bin_folder}}/*.fa; do
              bin_id=$(basename -s .fa ${{bin}})
              mv ${{bin}} ${{bin%/*}}/${{bin_id}}.fasta
            done
            shopt -u nullglob
            
            rm -rf {output.out_dir}/{params.prefix_1} \
                   {output.out_dir}/{params.prefix_2}/${{bin_folder}}/*/ \
                   {output.out_dir}/{params.prefix_2}/${{bin_folder}}_round_01*
            final_out={params.prefix_2}
        fi
        
        # Round 3
        if [[ ${{#binner_input_array[*]}} > 5 ]]; then
            mv {output.out_dir}/{params.prefix_2}/${{bin_folder}} {output.out_dir}/{params.prefix_2}/${{bin_folder}}_round_02
            CHECKM_DATA_PATH={params.checkm_db} SOFT={params.script_path} \
            {params.script_path}/bin_refinement.sh -o {output.out_dir}/{params.prefix_3} \
                                                   --quick \
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
            rm -rf {output.out_dir}/{params.prefix_2} \
                   {output.out_dir}/{params.prefix_3}/${{bin_folder}}/*/ \
                   {output.out_dir}/{params.prefix_3}/${{bin_folder}}_round_02*
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

        if test -n "$(find {output.out_dir}/figures -maxdepth 1 -name '*' -print -quit)"; then
          mv {output.out_dir}/figures/* Visualization/
        fi
        
        rm -rf {output.out_dir}/figures {output.out_dir}/bins{{A,B,C}}*
        
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

# Compress single binner out
# rule compress_binner_out:
#     input:
#         # Make sure its done at the very end
#         check = "status/workflow.done",
#         binner_out = check_binner_dirs
#     output:
#         "Binning/{binner}.tar.gz"
#     threads:
#         # Ensure at least 10gb of mem
#         max(int(math.ceil(10 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 8)
#     resources:
#         runtime = SHORT_RUNTIME,
#         mem_mb = MAIN_MEMCORE_MB
#     log:
#         "logs/misc/compress_binner_{binner}.log"
#     message:
#         "compress_binner_out: Compressing output of {wildcards.binner}."
#     conda:
#         os.path.join(ENVDIR, "IMP_preprocessing.yaml")
#     shell:
#        """
#        echo '.Compressing {wildcards.binner} output.'
#        tar -cv {input.binner_out} 2> {log} | pigz --stdout \
#                                                   --best \
#                                                   --processes {threads} \
#                                                   --rsyncable \
#                                                   --verbose > {output} 2>> {log} && \
#        rm -rf {input.binner_out} >> {log} 2>&1
#        """

