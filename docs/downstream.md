# Downstream: IMP explorer and sRNA-IMP

## IMP explorer

[IMP explorer](https://github.com/genomeanalysislcsb/imp-explorer) collects the IMP3 results of all samples of a
study into one release of tables (QC, taxonomy, MAGs and their abundance, functions) and serves them as an
interactive web app: composition, ordination, group tests, MAG quality and activity, functional potential vs
expression, and a search across layers.

```bash
imp-explorer export --imp-root /path/to/IMP/outputs --metadata metadata.tsv --study "My study" --out release_v0.1
imp-explorer run --data release_v0.1
```

Documentation: <https://genomeanalysislcsb.github.io/imp-explorer/>

## sRNA-IMP

[sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP) analyses small-RNA libraries (intra- and extracellular)
of the same samples on the IMP3 assembly, annotation and bins: human vs microbial small RNAs, tRNA / rRNA fragments,
novel sRNA loci, CRISPR RNAs, and comparisons with the metatranscriptome and metagenome.

Documentation: <https://genomeanalysislcsb.github.io/sRNA-IMP/>
