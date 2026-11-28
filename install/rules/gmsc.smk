##################################################
# GMSC (global microbial smORF catalogue, Duan et al. 2024) for imp3 smorf: GMSC-mapper, the 90 % family
# representatives with habitat / taxonomy / quality / domain annotation, and the search index
##################################################
rule gmsc_db:
    output:
        faa = "gmsc/GMSC10.90AA.faa.gz",
        path = "gmsc/gmsc-mapper.path",
        ready = touch("gmsc/gmsc.ready")
    params:
        tool = config["gmsc"].get("tool", "diamond")
    log:
        "gmsc/install.log"
    threads:
        8
    conda:
        os.path.join(ENV_DIR, "gmsc_mapper.yaml")
    message:
        "GMSC: GMSC-mapper database and {params.tool} index"
    shell:
        """
        gmsc-mapper downloaddb --dbdir gmsc --all > {log} 2>&1
        gmsc-mapper createdb -i {output.faa} -o gmsc -m {params.tool} >> {log} 2>&1
        command -v gmsc-mapper > {output.path}
        """
