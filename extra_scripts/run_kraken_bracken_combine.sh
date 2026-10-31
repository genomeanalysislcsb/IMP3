#!/bin/bash

# Define constants and configurations
DBPATH="/path/to/kraken/database"
ENVDIR="/path/to/conda/envs"
ROOTDIR="/path/to/project/root"
CMBREPDIR="/path/to/combine_kraken_reports.py"
OUTDIR="/path/to/output/directory" # Corrected from "${/path/to/output/directory}" to "/path/to/output/directory"
SAMPLES_LIST="/path/to/samples_list.txt" # Each line in this file should contain the sample name

KRAKENCONF=0.501
KRAKENTHREADS=60
KRAKENMINHIT=3

LEVELS="D P C O F G S" # Define all levels for Bracken to process
BRACKEN_KMER=35 # Replace with the actual k-mer length or method to determine it
OMES="mg mt" # Define the two omes to loop over

# Load Conda environment for IMP_taxonomy, if required
# source activate ${ENVDIR}/IMP_taxonomy.yaml

# Loop through each sample name in the samples list
while read -r SAMPLE; do
    for OME in $OMES; do
        echo "Processing $OME for sample: $SAMPLE"

        # Adjusted file paths to include $OME in the path/filename
        PE_R1="Preprocessing/${SAMPLE}.${OME}.r1.preprocessed.fq"
        PE_R2="Preprocessing/${SAMPLE}.${OME}.r2.preprocessed.fq"
        SE="Preprocessing/${SAMPLE}.${OME}.se.preprocessed.fq"
        KRAKEN_PE_OUTPUT="${OUTDIR}/Analysis/taxonomy/kraken/${SAMPLE}.${OME}.pe.kraken.output"
        KRAKEN_PE_REPORT="${OUTDIR}/Analysis/taxonomy/kraken/${SAMPLE}.${OME}.pe.kraken.report"
        KRAKEN_SE_OUTPUT="${OUTDIR}/Analysis/taxonomy/kraken/${SAMPLE}.${OME}.se.kraken.output"
        KRAKEN_SE_REPORT="${OUTDIR}/Analysis/taxonomy/kraken/${SAMPLE}.${OME}.se.kraken.report"
        COMBINED_REPORT="${OUTDIR}/Analysis/taxonomy/kraken/${SAMPLE}.${OME}.kraken.report"
        LOG_DIR="${OUTDIR}/logs/analysis/taxonomy"
        BRACKEN_DIR="${OUTDIR}/Analysis/taxonomy/bracken/${SAMPLE}.${OME}"

        # Ensure output directories exist
        mkdir -p ${OUTDIR}/Analysis/taxonomy/kraken ${BRACKEN_DIR} ${LOG_DIR}

        # Run Kraken on PE reads
        kraken2 --db ${DBPATH} \
                --confidence ${KRAKENCONF} \
                --threads ${KRAKENTHREADS} \
                --paired \
                --use-names \
                --minimum-hit-groups ${KRAKENMINHIT} \
                --output ${KRAKEN_PE_OUTPUT} \
                --report ${KRAKEN_PE_REPORT} \
                ${PE_R1} ${PE_R2} &>> ${LOG_DIR}/${SAMPLE}.${OME}.pe.log
        gzip -f ${KRAKEN_PE_OUTPUT}

        # Run Kraken on SE reads
        kraken2 --db ${DBPATH} \
                --confidence ${KRAKENCONF} \
                --threads ${KRAKENTHREADS} \
                --use-names \
                --minimum-hit-groups ${KRAKENMINHIT} \
                --output ${KRAKEN_SE_OUTPUT} \
                --report ${KRAKEN_SE_REPORT} \
                ${SE} &>> ${LOG_DIR}/${SAMPLE}.${OME}.se.log
        gzip -f ${KRAKEN_SE_OUTPUT}

        # Combine Kraken reports
        python ${CMBREPDIR} \
            -r ${KRAKEN_PE_REPORT} ${KRAKEN_SE_REPORT} \
            -o ${COMBINED_REPORT} \
            --only-combined \
            --no-headers \
            --has-minimizer &> ${LOG_DIR}/combine_kraken_reports.${SAMPLE}.${OME}.log


        # Loop over all levels for Bracken
        for LEVEL in $LEVELS; do
            echo "Running Bracken for level $LEVEL on $OME of sample $SAMPLE"
            BRACKEN_OUTPUT="${BRACKEN_DIR}/${SAMPLE}.${OME}.${LEVEL}.bracken"

            bracken -d ${DBPATH} \
                    -i ${COMBINED_REPORT} \
                    -o ${BRACKEN_OUTPUT} \
                    -r ${BRACKEN_KMER} \
                    -l ${LEVEL} \
                    -t ${KRAKENTHREADS} &>> ${LOG_DIR}/bracken_reads.${SAMPLE}.${OME}.log

            # Check if Bracken output exists and isn't empty
            if [ ! -s "${BRACKEN_OUTPUT}" ]; then
                echo "Error: Bracken output for ${OME} of ${SAMPLE} at level ${LEVEL} is empty."
                exit 1
            fi
        done

        echo "Finished processing ${OME} for ${SAMPLE}"
    done
done < "${SAMPLES_LIST}"

echo "All samples and omes processed."
