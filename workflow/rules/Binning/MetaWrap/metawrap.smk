import glob

localrules: MetaWrap_refine_check

checkpoint MetaWrap_refine_check:
    input:
        check_mw_input,
        "Binning/binny/contig2bin.tsv" if 'binny' in IMP_STEPS else []
    output:
        "status/MetaWrap_refine.ready"
    shell:
        "touch {output}"

rule get_checkm2_db:
    output:
        checkm2_db = f"{DBPATH}/checkm2/uniref100.KO.1.dmnd"
    conda:
        os.path.join(ENVDIR, "IMP_metawrap.yaml")
    log:
        os.path.join(ENVDIR, "logs/checkm2_downloadDB.log")
    shadow: "copy-minimal"
    shell:
        "checkm2 database --download --path {output.checkm2_db} >> {log} 2>&1"

rule MetaWrap_refine:
    input:
        binner_fasta_dirs = check_binner_bin_fastas,
        checkm2_db = f"{DBPATH}/checkm2/uniref100.KO.1.dmnd"
    output:
        "Binning/MetaWrap/bins.stats",
        directory("Binning/MetaWrap/bins"),
        out_dir = directory("Binning/MetaWrap"),
        bin2cont = "Binning/MetaWrap/bins.contigs"
    params:
        binner_inputs_list=get_mw_binner_inputs,
        mw_min_comp=config['binning']['MetaWrap']['min_completeness'],
        mw_max_cont=config['binning']['MetaWrap']['max_contamination'],
        prefix_1="metawrap_round_01",
        prefix_2="metawrap_round_02",
        prefix_3="metawrap_round_03",
        total_mem=int(float(int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]) * getThreads(MAIN_BIGCORENO, bigmem=True)),
        script_path = os.path.join(SRCDIR, "MW_scripts"),
        imp_out_dir = OUTPUTDIR,
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning/MetaWrap.log"
    message:
        "MetaWrap: Running MetaWrap."
    shadow: "copy-minimal"
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
        
        # Copy checkmdb into work dir
        wdir=$(pwd -P)
        cm2_db_file=$(basename {input.checkm2_db})
        rsync -aq {input.checkm2_db} $wdir/

        # Ensure Visualization folder exists
        mkdir -p {params.imp_out_dir}/Visualization

        # Round 1
        echo "########################################################################################" >> {log}
        echo "Running MetaWrap bin refinement with the binner outputs: {params.binner_inputs_list[0]}" >> {log}
        echo "########################################################################################" >> {log}

        SOFT={params.script_path} \
        {params.script_path}/bin_refinement.sh -o {output.out_dir}/{params.prefix_1} \
                                               --checkm2-db-path $wdir/$cm2_db_file \
                                               -t {threads} \
                                               -m {params.total_mem} \
                                               {params.binner_inputs_list[0]} \
                                               -c {params.mw_min_comp} \
                                               -x {params.mw_max_cont} >> {log} 2>&1

        # If no bins are in the output directory, create empty output files and exit
        shopt -s nullglob
        bin_files=( "{output.out_dir}/{params.prefix_1}/${{bin_folder}}"/*.fa )
        shopt -u nullglob

        if [[ ${{#bin_files[@]}} -eq 0 ]]; then
          echo "No bins above thresholds found. Creating empty output files and exiting." >> {log}
          mkdir -p {output.out_dir}/bins >> {log} 2>&1
          touch {output.out_dir}/bins.stats {output.out_dir}/bins.contigs >> {log} 2>&1
          exit 0
        fi

        # Rename bin fastas back to .fasta
        shopt -s nullglob
        for bin in {output.out_dir}/{params.prefix_1}/${{bin_folder}}/*.fa; do
          bin_id=$(basename -s .fa ${{bin}})
          mv ${{bin}} ${{bin%/*}}/${{bin_id}}.fasta
        done
        shopt -u nullglob

        # Move plot to Visualization
        if [[ -f {output.out_dir}/{params.prefix_1}/figures/binning_results.png ]]; then
          mv {output.out_dir}/{params.prefix_1}/figures/binning_results.png \
             {params.imp_out_dir}/Visualization/metawrap_round_01_binning_results.png
        fi
        rm -rf {output.out_dir}/{params.prefix_1}/${{bin_folder}}/*/
        
        # Round 2
        if [[ ${{#binner_input_array[*]}} > 3 ]]; then
            echo "########################################################################################" >> {log}
            echo "Running MetaWrap bin refinement with the binner outputs: {params.binner_inputs_list[1]}" >> {log}
            echo "########################################################################################" >> {log}

            mv {output.out_dir}/{params.prefix_1}/${{bin_folder}} {output.out_dir}/{params.prefix_1}/${{bin_folder}}_round_01
            SOFT={params.script_path} \
            {params.script_path}/bin_refinement.sh -o {output.out_dir}/{params.prefix_2} \
                                                   --checkm2-db-path $wdir/$cm2_db_file \
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
            
            # Move plot to Visualization
            if [[ -f {output.out_dir}/{params.prefix_2}/figures/binning_results.png ]]; then
              mv {output.out_dir}/{params.prefix_2}/figures/binning_results.png \
                 {params.imp_out_dir}/Visualization/metawrap_round_02_binning_results.png
            fi

            rm -rf {output.out_dir}/{params.prefix_1} \
                   {output.out_dir}/{params.prefix_2}/${{bin_folder}}/*/ \
                   {output.out_dir}/{params.prefix_2}/${{bin_folder}}_round_01*
            final_out={params.prefix_2}
        fi
        
        # Round 3
        if [[ ${{#binner_input_array[*]}} > 5 ]]; then

            echo "\n########################################################################################" >> {log}
            echo "Running MetaWrap bin refinement with the binner outputs: {params.binner_inputs_list[2]}" >> {log}
            echo "########################################################################################\n" >> {log}

            mv {output.out_dir}/{params.prefix_2}/${{bin_folder}} {output.out_dir}/{params.prefix_2}/${{bin_folder}}_round_02
            SOFT={params.script_path} \
            {params.script_path}/bin_refinement.sh -o {output.out_dir}/{params.prefix_3} \
                                                   --checkm2-db-path $wdir/$cm2_db_file \
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

            # Move plot to Visualization
            if [[ -f {output.out_dir}/{params.prefix_3}/figures/binning_results.png ]]; then
              mv {output.out_dir}/{params.prefix_3}/figures/binning_results.png \
                 {params.imp_out_dir}/Visualization/metawrap_round_03_binning_results.png
            fi

            rm -rf {output.out_dir}/{params.prefix_2} \
                   {output.out_dir}/{params.prefix_3}/${{bin_folder}}/*/ \
                   {output.out_dir}/{params.prefix_3}/${{bin_folder}}_round_02*
            final_out={params.prefix_3}
        fi
        
        # Clean up
        rm -rf {output.out_dir}/${{final_out}}/bins \
               {output.out_dir}/${{final_out}}/work_files \
               $wdir/$cm2_db_file
        mv {output.out_dir}/${{final_out}}/* {output.out_dir}/ && rm -rf {output.out_dir}/${{final_out}}

        for i in {output.out_dir}/metawrap_{params.mw_min_comp}_{params.mw_max_cont}*; do
          mv ${{i}} ${{i/metawrap_{params.mw_min_comp}_{params.mw_max_cont}_/}}
        done
        
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
