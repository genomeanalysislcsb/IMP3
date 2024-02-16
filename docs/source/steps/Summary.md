(step-summary)=

# IMP3 steps: Summary

IMP3 will generate an HTML report at the end of the pipeline to provide a comprehensive overview of the results generated. To see an example HTML report, please download the [archived_reports.zip](https://zenodo.org/records/161321/files/archived_reports.zip?download=1) from [this reference page](https://zenodo.org/records/161321), extract the folder an run the index.html to view the an example report in your browser. 

* Quality Control Metrics:
    * Visualizations and statistics regarding the quality of the input data, including read quality scores, GC content, and the results of preprocessing steps like trimming and filtering.

* Assembly Statistics:
    * Details about the assembly process, including metrics like the number of contigs produced, the N50 size, the total length of the assembly, and possibly comparisons of these metrics across different samples or assemblies.

* Binning Results:
    * Summary of binning outcomes, including the number of bins generated, their estimated completeness and contamination (often assessed by tools such as CheckM), and the taxonomic composition of these bins. This may include interactive charts or tables allowing users to explore the taxonomic diversity within the dataset.

* Taxonomic Classification:
    * A breakdown of the taxonomic composition of the microbial community identified in the samples. This could be presented through various interactive visualizations like pie charts, bar graphs, or hierarchical trees, showing the relative abundance of different taxonomic groups.

* Functional Annotation:
    * Summaries of the functional analysis, indicating the main functional categories or pathways identified within the assembled contigs. This could include data on gene families, metabolic pathways, and resistance genes, among others, often presented in an interactive format for detailed exploration.

* Comparative Analysis (if applicable):
    * For studies involving multiple samples or conditions, the report might include comparative analyses highlighting differences in microbial composition or function across the samples. This could be visualized through heatmaps, scatter plots, or other comparative graphs.

* Methodological Details:
    * A detailed description of the methods and parameters used in the pipeline, including software tools, databases, and specific settings, ensuring transparency and reproducibility of the analysis.

* Interactive Features:

    * The HTML format allows for interactive elements, such as the ability to zoom in on graphs, filter table views, or click through to more detailed information about specific contigs, bins, or taxa.

* Executive Summary:

   * An overview section summarizing the key findings and insights gained from the analysis, designed to provide a quick understanding of the most important results.