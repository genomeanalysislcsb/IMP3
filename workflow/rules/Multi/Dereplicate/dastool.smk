localrules: DASTool_check

checkpoint DASTool_check:
    input:
        lambda wildcards: expand("Dereplication/{binner}/{sample}/scaffold2bin.tsv",
        binner=config['dereplication']['cross_mapping_rebinning']['binners'].split(),
        sample=wildcards.sample)
    output:
        "Dereplication/{sample}/DAS_tool.ready"
    shell:
        "touch {output}"

rule DAS_tool_multi:
    input:
        contigs = lambda wildcards: expand("{raw_directory}/Assembly/{ass}.assembly.merged.fa",
                raw_directory=samples.path[wildcards.sample],
                ass=samples.assembly_type[wildcards.sample]),
        proteins="Dereplication/{sample}/prokka.renamed.faa", 
        tabs=check_bintabs_multi
    output:
        "Dereplication/{sample}/selected_DASTool_summary.txt",
        "Dereplication/{sample}/selected_DASTool_scaffolds2bin.txt",
        directory("Dereplication/{sample}/selected_DASTool_bins")
    threads: getThreads(4)
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/dereplication_DAS_tool.{sample}.log"
    message: "DAS_tool: Running DAS_Tool."
    shell:
        """
        TABS="{input.tabs}"
        TAB_tmp=${{TABS//Dereplication\//}}
        TAB_tmp2=${{TAB_tmp//\/{wildcards.sample}\/scaffold2bin.tsv/}}
        TABSl=${{TAB_tmp2// /,}}
        DAS_Tool -i ${{TABS// /,}} \
         -l $TABSl \
         -c {input.contigs} \
         --search_engine diamond \
         --proteins {input.proteins} --score_threshold 0.3 \
         --threads {threads} --write_bins 1 \
         -o Dereplication/{wildcards.sample}/selected > {log} 2>&1
        """

def input_for_protein_names(wildcards):
    if samples.bins[wildcards.sample]:
        inputs = expand("{raw_directory}/Binning/prokka.renamed.faa",
                  raw_directory=samples.path[wildcards.sample])
    else:
        inputs = expand("{raw_directory}/Analysis/annotation/{file}",
                           raw_directory=samples.path[wildcards.sample],
                           file=["prokka.faa","annotation.filt.gff"])
    return inputs

CP_PROTEIN_SHELL = """
                   cp {input} {output}
                   """
                
PREP_PROTEIN_SHELL = """
                     TMP_PROT=$(mktemp --tmpdir={TMPDIR} -t "proteins_XXXXXX.faa")
                     TMP_TAB=$(mktemp --tmpdir={TMPDIR} -t "proteins_anno_XXXX.tsv")
                     cut -f 1,9 {input[1]} | sort -k 1 | perl -F'\\t' -ane '$_=~/ID=(.+)_(.+);inference/; print "$F[0]\\t$1\\t$2\\n";' > $TMP_TAB
	                 sed 's/ .*//' {input[0]} >> $TMP_PROT
                     {SRCDIR}/prokka2dastoolHeader.pl  $TMP_PROT $TMP_TAB > {output} 2> {log}
                     """

rule get_protein_names:
    input:
        input_for_protein_names
    output:
        "Dereplication/{sample}/prokka.renamed.faa"
    resources:
        runtime = "4:00:00",
        mem = MEMCORE
    threads: 1
    log: "logs/dereplication_get_protein_names.{sample}.log"
    message: "get_protein_names: Preparing protein sequences for DAS tool."
    run:
        if len(input) > 1:
            shell(PREP_PROTEIN_SHELL)
        else:
            shell(CP_PROTEIN_SHELL)

