localrules: MetaWrap_refine_check

checkpoint MetaWrap_refine_check:
    input:
        expand("Binning/{binner}/bins",binner=BIN_STEPS)
    output:
        "Binning/MetaWrap_refine.ready"
    shell:
        "touch {output}"


def get_mw_binner_inputs(binner_fasta_dirs):
    out_list = []
    prefixes = ['-A ', '-B ', '-C ']
    round_1_length = min(len(binner_fasta_dirs), 3)
    input_round_1 = ' '.join([prefix + binner_fasta_path for prefix, binner_fasta_path in
                              zip(prefixes, binner_fasta_dirs[:round_1_length])])
    out_list.append(input_round_1)
    if len(binner_fasta_dirs) > 3:
        input_round_2 = ' '.join([prefix + binner_fasta_path for prefix, binner_fasta_path in
                                  zip(prefixes, ['Binning/MetaWeap/metawrap_round_01/bins'] + binner_fasta_dirs[3:])])
        out_list.append(input_round_2)
    return out_list


rule MetaWrap_refine:
    input:
        binner_fasta_dirs=check_binner_bin_fastas + ["Binning/DASTool/selected_DASTool_bins"]
    output:
        f"Binning/MetaWrap/metawrap_{params.mw_min_comp}_{params.mw_max_cont}_bins.stats",
        out_dir=directory("Binning/MetaWrap")
    params:
        binner_inputs_list=get_mw_binner_inputs(binner_fasta_dirs),
        mw_min_comp=config[binning][MetaWrap][min_completeness],
        mw_max_cont=config[binning][MetaWrap][max_contamination]
    threads: getThreads(8)
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/binning_DAS_tool.log"
    message: "DAS_tool: Running DAS_Tool."
    shell:
        """
        prefix_1="metawrap_round_01"
        prefix_2="metawrap_round_02"
        final_out=${prefix_1}
        bin_folder="metawrap_{params.mw_min_comp}_{params.mw_max_cont}_bins"
        read -a binner_input_array <<< {binner_fasta_dirs}
    
        # Round 1
        metawrap bin_refinement -o {out_dir}/${prefix_1}/ \
                        -t {threads} \
                        -m {mem} \
                        {params.binner_inputs_list[0]} \
                        -c {params.mw_min_comp} \
                        -x {params.mw_max_cont}
        
        # Rename bin fastas back to .fasta
        for bin in {out_dir}/${prefix_1}/${bin_folder}/*.fa; do
          bin_id=$(echo ${bin##*/} | cut -d "." -f 2)
          mv ${bin} ${bin%/*}/metawrap.${bin_id}.fasta
        done
        
        if [[ ${#binner_input_array[*]} > 3 ]]; then
            # Round 2
            metawrap bin_refinement -o {out_dir}/${prefix_1}/ \
                    -t {threads} \
                    -m {mem} \
                    {params.binner_inputs_list[1]} \
                    -c {params.mw_min_comp} \
                    -x {params.mw_max_cont}
            
            # Rename bin fastas back to .fasta
            for bin in {out_dir}/${prefix_2}/${bin_folder}/*.fa; do
              bin_id=$(echo ${bin##*/} | cut -d "." -f 2)
              mv ${bin} ${bin%/*}/metawrap.${bin_id}.fasta
            done
            rm -rf {out_dir}/${prefix_1}
            final_out=${prefix_2}
        fi
        
        # Clean up
        rm -rf {out_dir}/${final_out}/bins {out_dir}/${final_out}/${bin_folder}/binsA {out_dir}/${final_out}/work_files
        mv {out_dir}/${final_out}/* {out_dir}/ && rm -rf {out_dir}/${final_out}
        for i in Binning/MetaWrap/metawrap_{params.mw_min_comp}_{params.mw_max_cont}*; do
          mv ${i} ${i/{params.mw_min_comp}_{params.mw_max_cont}_/}
        done
        """


