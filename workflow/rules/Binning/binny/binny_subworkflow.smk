from snakemake.utils import min_version
min_version("6.0")

# Copy some params since Snakemake yaml doesn't allow param mapping
config["binning"]["binny"]['mem'] = {"big_mem_per_core_gb": int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]],
                                     "normal_mem_per_core_gb": int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]],
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
binny_s2b = "Binning/binny/contig2bin.tsv"


# First we check if theres enough contigs/the assembly is large enough to bin
BINNY_IMPOSSIBLE = False

checkpoint binny_check_assembly:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
    output:
        status_check = 'status/binning_binny_assembly.checked'
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(14)
    message:
        "binny: Checking assembly size."
    run:
        # Load assembly
        cont_dict = {}
        
        with open(input.assembly, 'r') as f:
            current_header = None
            header_line = False
            for line in f:
                if line.startswith(">"):
                    header_line = True
                else:
                    header_line = False
                line = line.strip(">\n\t ")
                # Get header
                if header_line:
                    # Transform sublists of sequences from previous header to one list.
                    if current_header:
                        cont_dict[current_header] = "".join(cont_dict[current_header])
                    cont_dict[line] = []
                    current_header = line
                # Get sequences
                else:
                    cont_dict[current_header].append(line)
            # Transform sublists of sequences from last header to one list.
            cont_dict[current_header] = "".join(cont_dict[current_header])

        # Check if assembly is large enough
        min_size = config["binning"]["binny"]["max_cont_length_cutoff_marker"]
        contigs_greater_min = [contig for contig in cont_dict if len(cont_dict[contig]) > min_size]
        n_contigs_greater_min = len(contigs_greater_min)
        contigs_greater_min_len_bp = sum([len(cont_dict[contig]) for contig in contigs_greater_min])

        if n_contigs_greater_min < 100 or contigs_greater_min_len_bp < 100000:
            BINNY_IMPOSSIBLE = True
        else:
            BINNY_IMPOSSIBLE = False
        
        # Write status file
        with open(output.status_check, 'w') as f:
            f.write("impossible" if BINNY_IMPOSSIBLE else "possible")


module binny_workflow:
    snakefile:
        os.path.join(ROOTDIR, "submodules/binny/Snakefile")
    config:
        config["binning"]["binny"]

use rule mantis_checkm_marker_sets from binny_workflow as binny_mantis_checkm_marker_sets with:
    input:
        status_check = 'status/binning_binny_assembly.checked',
        proteins="Analysis/annotation/bakta/annotation.faa"
    output:
        "Binning/tmp/consensus_annotation.tsv",
        out_dir = directory("Binning/tmp")
    params:
        binny_cfg=os.path.join(ROOTDIR, "submodules/binny/config/binny_mantis.cfg")
    threads:
        getThreads(64)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/binning/binny_mantis_checkm.log"
    message:
        "Mantis: Annotate {input} with CheckM markers"


use rule binny from binny_workflow as binny_binny with:
    input:
        status_check = 'status/binning_binny_assembly.checked',
        mgdepth = 'Stats/mg/%s.assembly.contig_depth.txt' % ASS,
        raw_gff = "Analysis/annotation/bakta/annotation.filt.gff",
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        hmm_markers = "Binning/tmp/consensus_annotation.tsv",
        tmp_dir = "Binning/tmp"
    output:
        os.path.join("Binning/binny/contig_data.tsv.gz"),
        os.path.join("Binning/binny/binny.done"),
        out_dir = directory('Binning/binny'),
        out_bin_dir = directory('Binning/binny/bins')
    params:
        binny_out=config["binning"]["binny"]["outputdir"],
        sample=config['sample'],
        py_functions=ROOTDIR + "/submodules/binny/workflow/scripts/binny_functions.py",
        t2p=ROOTDIR + "/submodules/binny/database/pfam/tigrfam2pfam.tsv",
        marker_sets=ROOTDIR + "/submodules/binny/database/taxon_marker_sets_lineage_sorted.tsv",
        gff = "Binning/tmp/annotation_CDS_RNA_hmms_checkm.gff",
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
        runtime = LONG_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/binning/binny.log"

# Check if we run binny
def binny_check(wildcards):
    with checkpoints.binny_check_assembly.get(**wildcards).output[0].open() as f:
        if f.read().strip() == "impossible":
            return 'status/binning_binny_assembly.checked'
        else:
            return [binny_contig_data, 'Binning/binny/bins']

rule binny_contig2bin:
    input:
        check = binny_check
    output:
        s2b = binny_s2b
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "binny: Preparing table for DAS_Tool from binny output."
    run:
        import gzip
        import os

        # If binny is impossible, we just touch and empty output file
        if input.check == 'status/binning_binny_assembly.checked':
            print("Binny impossible, touching empty file.")
            with open(output.s2b, 'w') as f:
                f.write('')
        else:
            print("Binny possible, preparing table.")
            with gzip.open(binny_contig_data,'rt') as gz: # {input.contig_data}
                with open(output.s2b, 'w') as of:  # {output.s2b}
                    header = next(gz)
                    for line in gz:
                        line = line.strip('\n ')
                        line_contig = line.split('\t')[0]
                        line_bin = line.split('\t')[1]
                        if line_bin != 'N':
                            of.write(f'{line_contig}\t{line_bin}\n')
