#!/bin/bash

# Define constants and configurations
DBNAME="k2_nt_20231129"
DBPATH="/dev/shm/${DBNAME}"
ENV="kraken"
ROOTDIR="/mnt/data"
CMBREPDIR="${ROOTDIR}/local_tools/KrakenTools/combine_kreports.py"
OUTDIR="${ROOTDIR}/oskar/projects/covalux/Kraken_Bracken"
SAMPLE_DIR="/home/oskar.hickl/isilon/projects/ecosystem_biology/datatransfer/Laetitia_preprocessed_reads"
SAMPLES_LIST="${ROOTDIR}/oskar/projects/covalux/sample_list.txt" # Each line in this file should contain the sample name

FAILED_FILES_LIST="${OUTDIR}/failed_files.txt"

KRAKENCONF=0.501
KRAKENTHREADS=60
KRAKENMINHIT=3

LEVELS="D P C O F G S" # Define all levels for Bracken to process
BRACKEN_KMER=50 # Replace with the actual k-mer length or method to determine it
OMES="mg mt" # Define the two omes to loop over

# Load Conda environment for IMP_taxonomy, if required
# mamba init
# source ~/.bashrc
# mamba activate ${ENV}

# Loop through each sample name in the samples list
while read -r SAMPLE; do
    for OME in $OMES; do
        echo "Processing $OME for sample: $SAMPLE"

        # Adjusted file paths to include $OME in the path/filename
        PE_R1="${SAMPLE_DIR}/${SAMPLE}/${OME}.r1.preprocessed.fq"
        PE_R2="${SAMPLE_DIR}/${SAMPLE}/${OME}.r2.preprocessed.fq"
        SE="${SAMPLE_DIR}/${SAMPLE}/${OME}.se.preprocessed.fq"

        # Ensure read files exist
        if [ ! -s "${PE_R1}" ] || [ ! -s "${PE_R2}" ] || [ ! -s "${SE}" ]; then
            echo "Warning: One or more input files for ${OME} of ${SAMPLE} do not exist or are empty."
            # Skip sample ome combo
            continue
        fi

        KRAKEN_DIR="${OUTDIR}/${SAMPLE}/kraken/${OME}"
        KRAKEN_PE_OUTPUT="${KRAKEN_DIR}/${OME}.pe.kraken.output"
        KRAKEN_PE_REPORT="${KRAKEN_DIR}/${OME}.pe.kraken.report"
        KRAKEN_SE_OUTPUT="${KRAKEN_DIR}/${OME}.se.kraken.output"
        KRAKEN_SE_REPORT="${KRAKEN_DIR}/${OME}.se.kraken.report"
        COMBINED_REPORT="${OUTDIR}/${SAMPLE}/kraken/${OME}/${OME}.kraken.report"
        LOG_DIR="${OUTDIR}/logs"
        BRACKEN_DIR="${OUTDIR}/${SAMPLE}/bracken/${OME}"

        # Check if output exists and isnt empty and skip kraken if so
        if [ -s "${COMBINED_REPORT}" ]; then
            echo "Kraken report for ${OME} of ${SAMPLE} already exist and are not empty. Skipping Kraken."
        else
            # Ensure output directories exist
            mkdir -p ${KRAKEN_DIR} ${BRACKEN_DIR} ${LOG_DIR}

            # Run Kraken on PE reads
            kraken2 --db ${DBPATH} \
                    --confidence ${KRAKENCONF} \
                    --threads ${KRAKENTHREADS} \
                    --memory-mapping \
                    --report-minimizer-data \
                    --paired \
                    --use-names \
                    --minimum-hit-groups ${KRAKENMINHIT} \
                    --output ${KRAKEN_PE_OUTPUT} \
                    --report ${KRAKEN_PE_REPORT} \
                    ${PE_R1} ${PE_R2} &>> ${LOG_DIR}/${SAMPLE}.${OME}.kraken.pe.log
            gzip -f ${KRAKEN_PE_OUTPUT}

            # Run Kraken on SE reads
            kraken2 --db ${DBPATH} \
                    --confidence ${KRAKENCONF} \
                    --threads ${KRAKENTHREADS} \
                    --memory-mapping \
                    --report-minimizer-data \
                    --use-names \
                    --minimum-hit-groups ${KRAKENMINHIT} \
                    --output ${KRAKEN_SE_OUTPUT} \
                    --report ${KRAKEN_SE_REPORT} \
                    ${SE} &>> ${LOG_DIR}/${SAMPLE}.${OME}.kraken.se.log
            gzip -f ${KRAKEN_SE_OUTPUT}

            # Combine Kraken reports
            python ${CMBREPDIR} \
                -r ${KRAKEN_PE_REPORT} ${KRAKEN_SE_REPORT} \
                -o ${COMBINED_REPORT} \
                --only-combined \
                --no-headers \
                --has-minimizer &> ${LOG_DIR}/${SAMPLE}.${OME}.combine_kraken_reports.log && \
                rm -rf ${KRAKEN_PE_REPORT} ${KRAKEN_SE_REPORT}
        fi

        # Loop over all levels for Bracken
        for LEVEL in $LEVELS; do
            BRACKEN_OUTPUT="${BRACKEN_DIR}/${OME}.${LEVEL}.bracken"
            # Check if combined report exists and isn't empty and skip bracken if so
            if [ -s "${BRACKEN_OUTPUT}" ]; then
                echo "Bracken report for ${OME} of ${SAMPLE} ${LEVEL} already exists and is not empty. Skipping Bracken."
            else
                echo "Running Bracken for level $LEVEL on $OME of sample $SAMPLE"

                bracken -d ${DBPATH} \
                        -i ${COMBINED_REPORT} \
                        -o ${BRACKEN_OUTPUT} \
                        -r ${BRACKEN_KMER} \
                        -l ${LEVEL} \
                        -t ${KRAKENTHREADS} &>> ${LOG_DIR}/${SAMPLE}.${OME}.bracken.log

                # Check if Bracken output exists and isn't empty
                if [ ! -s "${BRACKEN_OUTPUT}" ]; then
                    echo "Error: Bracken output for ${OME} of ${SAMPLE} at level ${LEVEL} is empty."
                    echo "${SAMPLE} ${OME} ${LEVEL}" >> ${FAILED_FILES_LIST}
                fi
            fi
        done

        echo "Finished processing ${OME} for ${SAMPLE}"
    done
done < "${SAMPLES_LIST}"

echo "All samples and omes processed."
