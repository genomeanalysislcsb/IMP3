# rules from EukDetect, https://github.com/allind/EukDetect
# re-worked to run with IMP-configuration by Anna Heintz-Buschart

# original license:
#MIT License

#Copyright (c) 2020 Abigail Lind

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

#targets

localrules: ED_alncmd, ED_runall, ED_printaln, ED_aln, ED_filter

rule ED_runall:
    input:
        expand("{output_dir}/aln/{sample}_aln_q30_lenfilter.sorted.bam", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"),
        expand("{output_dir}/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"),
        expand("{output_dir}/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam.bai", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"),
        expand("{output_dir}/filtering/{sample}_read_counts_and_mismatches.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"), 
        expand("{output_dir}/filtering/{sample}_all_hits_table.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"), 
        expand("{output_dir}/{sample}_filtered_hits_table.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"), 
        expand("{output_dir}/{sample}_filtered_hits_taxonomy.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect")
    output:
        touch("status/eukdetect.done")

rule ED_printaln:
    input:
        expand("{output_dir}/{sample}_alignment_commands.txt", output_dir="Analysis/taxonomy/EukDetect",sample=TYPES)

rule ED_aln:
    input:
        expand("{output_dir}/aln/{sample}_aln_q30_lenfilter.sorted.bam", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect")

rule ED_filter:
    input:
        expand("{output_dir}/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"),
        expand("{output_dir}/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam.bai", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"),
        expand("{output_dir}/filtering/{sample}_read_counts_and_mismatches.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"), 
        expand("{output_dir}/filtering/{sample}_all_hits_table.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"), 
        expand("{output_dir}/{sample}_filtered_hits_table.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect"),
        expand("{output_dir}/{sample}_filtered_hits_taxonomy.txt", sample=TYPES, output_dir="Analysis/taxonomy/EukDetect")

#rules

def alninputs(wildcards):
    if wildcards.sample == "mg" and len(MG)==2:
        input = ["Preprocessing/mg.r1.preprocessed.fq","Preprocessing/mg.r2.preprocessed.fq"]
    if wildcards.sample == "mt" and len(MT)==2:
        input = ["Preprocessing/mt.r1.preprocessed.fq","Preprocessing/mt.r2.preprocessed.fq"] 
    if wildcards.sample == "mg" and len(MG)==1:
        input = "Preprocessing/mg.se.preprocessed.fq"
    if wildcards.sample == "mt" and len(MT)==1:
        input = "Preprocessing/mt.se.preprocessed.fq"
    return input

def construct_input(wildcards,input):
    if len(input) == 2:
        outs = "-1 " + input[0] + " -2 " + input[1]
    else:
        outs = "-U " + input[0]
    return outs

rule ED_alncmd:
    input:
        alninputs
    output:
        expand("{output_dir}/{{sample}}_alignment_commands.txt", output_dir="Analysis/taxonomy/EukDetect")
    params:
        bam = "Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter.sorted.bam",
        minreadlen = config["eukdetect"]["min_readlen"],
        db = DBPATH + "/" + config["eukdetect"]["database_dir"] + "/all_buscos_v4.fna",
        inputstring = lambda wildcards,input: construct_input(wildcards,input)
    shell:
        """
        echo bowtie2 --quiet --omit-sec-seq --no-discordant --no-unal -x {params.db} {params.inputstring}  \|
         perl -lane {{'}}$l =0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > {params.minreadlen} or /^@/{{'}} \| 
         samtools view -q 30 -bS - \| 
         samtools sort -o {params.bam} -  >> {output}
        """

rule ED_runaln:
    input:
        alninputs
    output:
        "Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter.sorted.bam"
    params:
        inputstring = lambda wildcards,input: construct_input(wildcards,input),
        minreadlen = config["eukdetect"]["min_readlen"],
        db = DBPATH + "/" + config["eukdetect"]["database_dir"] + "/all_buscos_v4.fna"
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_aln.{sample}.log"
    message: "runaln: Running mapping for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        bowtie2 --quiet --omit-sec-seq --no-discordant --no-unal -x {params.db} {params.inputstring} 2>> {log} |\
         perl -lane '$l =0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > {params.minreadlen} or /^@/' |\
         samtools view -q 30 -bS - 2>> {log} |\
         samtools sort -o {output} - 2>> {log}
        """


rule ED_bam2fastq:
    input:
        bam = "Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter.sorted.bam"
    output:
        temp("Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter.fq")
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_bam2fastq.{sample}.log"
    message: "bam2fastq: Running bam2fastq for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        bedtools bamtofastq -i {input.bam} -fq {output} &>> {log}
        """

rule ED_find_low_complexity:
    input:
        fq = "Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter.fq"
    output:
        temp("Analysis/taxonomy/EukDetect/aln/{sample}_low_complexity.txt")
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
#    log: "logs/analysis/eukdetect_find_low_complexity.{sample}.log"
    message: "find_low_complexity: Finding low complexity for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        cat {input.fq} | kz | awk '{{ if ($4<0.5) print $1 }}' > {output}
        """

rule ED_remove_low_complexity:
    input:
        inbam = "Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter.sorted.bam",
        drop_complex ="Analysis/taxonomy/EukDetect/aln/{sample}_low_complexity.txt"
    output:
        temp("Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter_complexity.sorted.bam")
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_remove_low_complexity.{sample}.log"
    message: "remove_low_complexity: Removing low complexity for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        samtools view -h {input.inbam} 2>> {log} | grep -v -f {input.drop_complex} 2>> {log} |\
         samtools view -bS - > {output} 2>> {log}
        """


rule ED_fixmate:
    input:
        inbam = "Analysis/taxonomy/EukDetect/aln/{sample}_aln_q30_lenfilter_complexity.sorted.bam"  
    output:
        temp("Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexity_fixmate.bam")
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_fixmate.{sample}.log"
    message: "fixmate: Running fixmate for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        samtools sort -n {input.inbam} 2>> {log} | samtools fixmate -m - {output} 2>> {log}
        """

rule ED_markdup:
    input:
        "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexity_fixmate.bam"
    output:
        temp("Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexity_rmdup.bam")
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_markdup.{sample}.log"
    message: "markdup: Running markdup for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        samtools sort {input} 2>> {log} | samtools markdup -r -s - {output} 2>> {log}
        """

rule ED_rmsort:
    input:
        "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexity_rmdup.bam"
    output:
        "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam"
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_rmsort.{sample}.log"
    message: "rmsort: Running rmsort for eukdetect on {wildcards.sample} reads."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        samtools sort -o {output} {input} &>> {log}
        """

if not 'assembly' in IMP_STEPS and not 'analysis' in IMP_STEPS:
    rule ED_index:
        input:
            bam = "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam"
        output:
            bai = "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam.bai"
        threads: 1
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        log: "logs/analysis/eukdetect_index.{sample}.log"
        message: "index: Indexing bam for eukdetect on {wildcards.sample} reads."
        conda: ENVDIR + "/eukdetect.yaml"
        shell:
            """
            samtools index {input.bam} {output.bai} &>> {log}
            """

rule ED_countreads:
    input:
        bam = "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam",
        bai = "Analysis/taxonomy/EukDetect/filtering/{sample}_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam.bai"
    output:
        counts = "Analysis/taxonomy/EukDetect/filtering/{sample}_read_counts_and_mismatches.txt"
    params:
        countscript = expand("{eukdetect_dir}/eukdetect/bam_to_pid.py", eukdetect_dir=config["eukdetect"]["eukdetect_dir"]),
        ref = DBPATH + "/" + config["eukdetect"]["database_dir"] + "/all_buscos_v4.fna"
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_countreads.{sample}.log"
    message: "countreads: Counting {wildcards.sample} reads for eukdetect."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        python {params.countscript} {input.bam} {params.ref} > {output.counts} 2>> {log}
        """

rule ED_taxonomize:
    input:
        counts = "Analysis/taxonomy/EukDetect/filtering/{sample}_read_counts_and_mismatches.txt"
    output:
        allhits = "Analysis/taxonomy/EukDetect/filtering/{sample}_all_hits_table.txt",
        primary_hits = "Analysis/taxonomy/EukDetect/{sample}_filtered_hits_table.txt",
        primary_taxonomy ="Analysis/taxonomy/EukDetect/{sample}_filtered_hits_taxonomy.txt"
    params:
        script = expand("{eukdetect_dir}/eukdetect/count_primary_and_secondary.py", eukdetect_dir=config["eukdetect"]["eukdetect_dir"]),
        taxid_link = DBPATH + "/" + config["eukdetect"]["database_dir"] + "/busco_taxid_link.txt",
        inherited_markers = DBPATH + "/" + config["eukdetect"]["database_dir"] + "/specific_and_inherited_markers_per_taxid.txt",
        taxdb = DBPATH + "/" + config["eukdetect"]["database_dir"] + "/taxa.sqlite"
    threads: 1
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/analysis/eukdetect_taxonomize.{sample}.log"
    message: "taxonomize: Taxonomy of {wildcards.sample} reads for eukdetect."
    conda: ENVDIR + "/eukdetect.yaml"
    shell:
        """
        python {params.script} --dbfile {params.taxdb} --taxid_link {params.taxid_link} \
         --inherited_markers {params.inherited_markers} --readcounts {input.counts} \
         --primarytab {output.primary_hits} --primarytax {output.primary_taxonomy} --alltab {output.allhits}
        """


#rule filterpass:
#input:
#counts = "{output_dir}/filtering/{sample}_read_counts_and_mismatches.txt",
#script = expand("{eukdetect_dir}/eukdetect/taxonomy_counts_with_coverage_and_pid_filterpass.py", eukdetect_dir=config["eukdetect_dir"]),
#taxid_link = expand("{tax_dir}/busco_taxid_link.txt", tax_dir=config["database_dir"]),
#inherited_markers = expand("{tax_dir}/specific_and_inherited_markers_per_taxid.txt", tax_dir=config["database_dir"]),
#taxdb = expand("{tax_dir}/taxa.sqlite", tax_dir=config["database_dir"])
#output:
#taxonomy ="{output_dir}/{sample}_hit_taxonomy_filterpass.txt",
#readcounts = "{output_dir}/{sample}_stats_per_filtered_taxid.txt"
#shell:
#"python {input.script} {input.taxid_link} {input.taxdb} {input.inherited_markers} {input.counts} {output.taxonomy} {output.readcounts}"
