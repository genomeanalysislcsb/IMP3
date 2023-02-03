rule abawaca:
    input:
       assembly="Assembly/%s.assembly.merged.fa" % ASS,
       depth="Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        s2b="Binning/abawaca/scaffold2bin.tsv",
        out_dir=directory("Binning/abawaca"),
        out_bin_dir=directory("Binning/abawaca/bins")
    params:
        esom_source=os.path.join(ROOTDIR, "submodules/ESOM"),
        esom_min=config['binning']['abawaca']['esom_min'],
        esom_max=config['binning']['abawaca']['esom_max']
    threads:
        getThreads(32)
    resources:
        runtime = "24:00:00",
        mem_mb = MEMCORE_MB
    conda:
        ENVDIR + "/IMP_binning.yaml"
    log:
        "logs/binning_abawaca.log"
    message:
        "abawaca: Running abawaca."
    shell:
        """
        # Create out dir
        mkdir -p {output.out_bin_dir}
        cp {input.assembly} {output.out_dir}/assembly.fasta
        
        # Run ESOM
        echo 'Running ESOM.' >> {log} 2>&1
        perl {params.esom_source}/esomWrapper.pl -path {output.out_dir} \
                                                 -dir {output.out_dir}/ESOM \
                                                 -ext fasta \
                                                 -min {params.esom_min} \
                                                 -max {params.esom_max} >> {log} 2>&1
        
        # Run abawaca
        echo 'Running abawaca.' >> {log} 2>&1
        abawaca {output.out_dir}/ESOM/Tetra_esom_{params.esom_min}.names \
                {output.out_dir}/ESOM/Tetra_esom_{params.esom_min}.lrn \
                {output.out_dir}/ESOM/Tetra_esom_{params.esom_min}_{params.esom_max}_split.fasta\
                {output.out_dir}/out >> {log} 2>&1 && \
        mv {output.out_dir}/out/* {output.out_dir} && rm -rf {output.out_dir}/out &&\
        tail -n +2 {output.out_dir}/scaf2cluster.txt > {output.s2b} &&\
        {SRCDIR}/write_fastas_from_c2b_table.py {output.out_dir}/assembly.fasta \
                                                {output.s2b} \
                                                {output.out_bin_dir} \
                                                ',' \
                                                False \
                                                abawaca &&\
        rm -rf {output.out_dir}/assembly.fasta {output.out_dir}/final-clusters {output.out_dir}/clusters \
               {output.out_dir}/dp2cluster.txt {output.out_dir}/scaf2cluster.txt
        """