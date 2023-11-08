rule trimmomatic_adapters:
    output:
        touch("adapters/adapters.done")
    shell:
        """
        wget --no-check-certificate {config[trimmomatic][pkg_url]} -O Trimmomatic-Src-0.32.zip
        unzip Trimmomatic-Src-0.32.zip
        cp -r trimmomatic-0.32/adapters .
        rm Trimmomatic-Src-0.32.zip && rm -rf trimmomatic-0.32
        touch {output}
        """

rule sortmerna_databases:
    output:
        expand("sortmerna/{files}.fasta", files=config["sortmerna"]["files"])
    params:
        tmp_dir=TMPDIR,
        pkg_url=config["sortmerna"]["pkg_url"],
    shell:
        """
        export TMPDIR="{params.tmp_dir}"
        TMPD=$(mktemp -d -t "XXXXXX")

        wget {params.pkg_url} --no-check-certificate -O $TMPD/sortmerna.tgz
        tar -xzf  $TMPD/sortmerna.tgz --strip-components=1 -C $TMPD

        mkdir -p sortmerna
        mv $TMPD/rRNA_databases/*.fasta sortmerna/.
        rm -rf $TMPD
        """


rule sortmerna_databases_index:
    input:
        expand("sortmerna/{files}.fasta", files=config["sortmerna"]["files"])
    output:
        expand(
            "sortmerna/idx/{files}.{ext}",
            files=config["sortmerna"]["files"],
            ext=['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats'])
    run:
        fastaindexed = expand(
            "sortmerna/idx/{files}",
            files=config["sortmerna"]["files"])
        ref = ':'.join('%s,%s' % (a, b) for a, b in zip(input, fastaindexed))
        shell("mkdir -p sortmerna")
        shell("indexdb_rna --ref {ref}")

if 'url' in config['filtering'] and config['filtering']['url']:
    rule filtering_database_from_url:
        output:
            'filtering/{filter}.fa'.format(filter=FILTER)
            #expand("filtering/{filter}.{ext}", filter=FILTER, ext=['fa'])
        params:
            tmp_dir=TMPDIR,
            filter = FILTER
        shell:
            """
            export TMPDIR="{params.tmp_dir}"
            TMPD=$(mktemp -d -t "XXXXXX")
            wget {config[filtering][url]} --no-check-certificate -O $TMPD/{params.filter}.fa.gz
            gunzip $TMPD/{params.filter}.fa.gz
            mkdir -p filtering
            mv $TMPD/{params.filter}.fa filtering/.
            rm -rf $TMPD
            """

if 'path' in config['filtering'] and config['filtering']['path']:
    rule filtering_database_from_path:
        output:
            'filtering/{filter}.fa'.format(filter=FILTER)
        params:
            filter = FILTER
        shell:
            """
            cp {config[filtering][path]} filtering/{params.filter}.fa
            mkdir -p filtering
            mv $TMPD/{params.filter}.fa filtering/.
            """


rule filtrering_database_index:
    input:
        "{fasta}"
    output:
        "{fasta}.bwt.2bit.64"
    shell:
        """
        bwa-mem2 index {wildcards.fasta}
        """


rule kegg_information:
    output:
        "ec2pathway.txt",
        "pathway2hierarchy.txt"
    shell:
        """
        wget --no-check-certificate {config[kegg][db_ec2pthy]} -O ec2pathway.txt.tmp
        grep "path:ec" ec2pathway.txt.tmp | sed -e 's/path:ec//g' | sed -e 's/ec://g' > {output[0]}
        rm ec2pathway.txt.tmp

        python {SRCDIR}/make.pwy.hierarchy.kegg.py > {output[1]}
        """

rule prokka_databases:
    output:
        expand("{db}", db=config["prokka"]["databases"]),
        "hmm/essentials_genes.hmm"
    params:
        tmp_dir=TMPDIR
    shell:
        """
        ### prokka by default will look databases where is located the binary.
        ### we have to softlink to put the binary somewhere and the databases somewhere else.
        ### the author don't want to add an option to locate the database.
        if [[ "{DBPATH}" = /* ]]
        then
            PP={DBPATH};
        else
            PP=$PWD/{DBPATH};
        fi
        cd $(dirname $(which prokka))/.. && ln -fs $PP db

        echo "Softlinking $(dirname $(which prokka))/../db to $PP"
        export TMPDIR="{params.tmp_dir}"
        TMPDIR=$(mktemp -d -t "XXXXXX")
        wget {config[prokka][pkg_url]} --no-check-certificate -O $TMPDIR/prokka.tgz
        tar -xzf $TMPDIR/prokka.tgz --strip-components=1 -C $TMPDIR
        mkdir -p {DBPATH}
        cp -r $TMPDIR/db/* {DBPATH}/.
        rm -rf $TMPDIR
        wget {config[prokka][essentials_genes_db_url]} --no-check-certificate -O {DBPATH}/hmm/essentials_genes.hmm
        prokka --setupdb
        """

rule checkm_data:
    output:
        "checkM/selected_marker_sets.tsv"
    shell:
        """
        wget -c --no-check-certificate {config[checkm][data_url]} -O checkm_data.tgz
        mkdir -p {DBPATH}/checkM
        tar -xvzf checkm_data.tgz -C {DBPATH}/checkM
        rm checkm_data.tgz
        """
