from snakemake.utils import min_version
min_version("6.0")

# Copy some params since Snakemake yaml doesn't allow param mapping
config["binning"]["binny"]["mem"] = {"big_mem_per_core_gb": int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]],
                                     "mem_per_core": int(MEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]],
                                     "big_mem_avail": True}
config["binning"]["binny"]["tmp_dir"] = config["tmp_dir"]
config["binning"]["binny"]["sample"] = config["sample"]
config["binning"]["binny"]["db_path"] = os.path.join(ROOTDIR, "submodules/binny/database")
config["binning"]["binny"]["snakemake_env"] = os.path.join(ROOTDIR, "conda/snakemake_env")
config["binning"]["binny"]["prokka_env"] = os.path.join(ROOTDIR, "workflow/envs/IMP_annotation.yaml")
config["binning"]["binny"]["mantis_env"] = os.path.join(ROOTDIR, "workflow/envs/IMP_mantis.yaml")
config["binning"]["binny"]["outputdir"] = os.path.join( "Binning/binny")
config["binning"]["binny"]["raws"] = {"assembly": f"Assembly/{ASS}.assembly.merged.fa",
                                      "contig_depth": f"Stats/mg/{ASS}.assembly.contig_depth.txt"}

# Get conda env dir
config["binning"]["binny"]["conda_source"] = workflow.conda_prefix if workflow.conda_prefix else srcdir("../conda")

config["binning"]["binny"]["write_contig_data"] = "True"

binny_contig_data = os.path.join(config["binning"]["binny"]["outputdir"], "contig_data.tsv.gz")
binny_s2b = "Binning/binny/scaffold2bin.tsv"

module binny_workflow:
    snakefile:
        os.path.join(ROOTDIR, "submodules/binny/Snakefile")
    config:
        config["binning"]["binny"]

# use rule * from binny_workflow as binny_*

# localrules: prepare_input_data, binny_prepare_input_data

use rule mantis_checkm_marker_sets from binny_workflow as binny_mantis_checkm_marker_sets with:
    input:
        proteins="Analysis/annotation/prokka.faa"
    output:
        "Binning/binny/intermediary/mantis_out/consensus_annotation.tsv",
        out_dir=directory("Binning/binny/intermediary/mantis_out")
    params:
        binny_cfg=os.path.join(ROOTDIR, "submodules/binny/config/binny_mantis.cfg")
    threads:
        getThreads(64)
    resources:
        runtime="24:00:00",
        mem_mb = MEMCORE_MB
    log:
        "logs/binning_binny_mantis_checkm.log"
    message:
        "Mantis: Annotate {input} with CheckM markers"


use rule binny from binny_workflow as binny_binny with:
    input:
        mgdepth='Stats/mg/%s.assembly.contig_depth.txt' % ASS,
        raw_gff="Analysis/annotation/annotation.filt.gff",
        assembly="Assembly/%s.assembly.merged.fa" % ASS,
        hmm_markers="Binning/binny/intermediary/mantis_out/consensus_annotation.tsv"
    output:
        directory('Binning/binny/bins'),
        os.path.join("Binning/binny/contig_data.tsv.gz"),
        os.path.join("Binning/binny/binny.done")
    params:
        binny_out=config["binning"]["binny"]["outputdir"],
        sample=config['sample'],
        py_functions=ROOTDIR + "/submodules/binny/workflow/scripts/binny_functions.py",
        t2p=ROOTDIR + "/submodules/binny/database/pfam/tigrfam2pfam.tsv",
        marker_sets=ROOTDIR + "/submodules/binny/database/taxon_marker_sets_lineage_sorted.tsv",
        gff=os.path.join(config["binning"]["binny"]["outputdir"], "intermediary/annotation_CDS_RNA_hmms_checkm.gff"),
        min_completeness=config["binning"]["binny"]["bin_quality"]["min_completeness"],
        start_completeness=config["binning"]["binny"]["bin_quality"]["start_completeness"],
        purity=config["binning"]["binny"]["bin_quality"]["purity"],
        kmers=config["binning"]["binny"]["kmers"],
        mask_disruptive_sequences=config["binning"]["binny"]["mask_disruptive_sequences"],
        extract_scmags=config["binning"]["binny"]["extract_scmags"],
        coassembly_mode=config["binning"]["binny"]["coassembly_mode"],
        min_cutoff=config["binning"]["binny"]["min_cont_length_cutoff"],
        max_cutoff=config["binning"]["binny"]["max_cont_length_cutoff"],
        min_cutoff_marker=config["binning"]["binny"]["min_cont_length_cutoff_marker"],
        max_cutoff_marker=config["binning"]["binny"]["max_cont_length_cutoff_marker"],
        nx_val=config["binning"]["binny"]["NX_value"],
        max_n_contigs=config["binning"]["binny"]["max_n_contigs"],
        max_marker_lineage_depth_lvl=config["binning"]["binny"]["max_marker_lineage_depth_lvl"],
        distance_metric=config["binning"]["binny"]["distance_metric"],
        max_embedding_tries=config["binning"]["binny"]["embedding"]["max_iterations"],
        include_depth_initial=config["binning"]["binny"]["clustering"]["include_depth_initial"],
        include_depth_main=config["binning"]["binny"]["clustering"]["include_depth_main"],
        hdbscan_min_samples_range=config["binning"]["binny"]["clustering"]["hdbscan_min_samples_range"],
        hdbscan_epsilon_range=config["binning"]["binny"]["clustering"]["hdbscan_epsilon_range"],
        write_contig_data="True"
    threads:
        getThreads(64)
    resources:
        runtime = "48:00:00",
        mem_mb = MEMCORE_MB
    log:
        "logs/binning_binny.log"

rule binny_contig2bin:
    input:
        binny_contig_data
    output:
        binny_s2b
    resources:
        runtime = "01:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    message:
        "binny: Preparing table for DAS_Tool from binny output."
    run:
        import gzip
        with gzip.open(binny_contig_data,'rt') as gz: # {input.contig_data}
            with open(binny_s2b, 'w') as of:  # {output.s2b}
                header = next(gz)
                for line in gz:
                    line = line.strip('\n ')
                    line_contig = line.split('\t')[0]
                    line_bin = line.split('\t')[1]
                    if line_bin != 'N':
                        of.write(f'{line_contig}\t{line_bin}\n')
