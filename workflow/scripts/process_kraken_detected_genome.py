#!/usr/bin/env python3

import os
import gzip
import re
import glob
import sys
from pathlib import Path
import subprocess
import multiprocessing

import pyrodigal

def log(message, log_file):
    """Logs a message to both the log file and standard output."""
    with open(log_file, 'a') as f:
        f.write(message + '\n')
    print(message)

def is_gzip_valid(filepath):
    """Checks if a gzip file is valid."""
    try:
        with gzip.open(filepath, 'rb') as f:
            f.read(1)
        return True
    except gzip.BadGzipFile:
        return False

def download_file(url, dest_path, log_file):
    """Downloads a file using wget and logs the process."""
    try:
        # Using subprocess to call wget for efficient downloading
        cmd = f"wget --tries=5 --wait=10 --random-wait --retry-connrefused " \
              f"--no-verbose --retry-on-http-error=404,500,503 -O {dest_path} {url}"
        log(f"Downloading {url} to {dest_path} ...", log_file)
        ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if ret.returncode != 0:
            log(f"Failed to download {url}. Return code: {ret.returncode}", log_file)
            return False
        return True
    except Exception as e:
        log(f"Exception occurred while downloading {url}: {e}", log_file)
        return False

def run_gene_finder_on_sequence(args):
    """Runs Pyrodigal gene finder on a single sequence."""
    finder, sequence_id, sequence = args
    return (sequence_id, finder.find_genes(sequence))

def process_genome():
    """Main function to process the genome."""
    # Extract Snakemake inputs and outputs
    genome_id = snakemake.params.genome_id
    genome_accession = snakemake.params.genome_accession
    taxon_base_path = snakemake.params.taxon_base_path
    unsupported_taxa = snakemake.params.unsupported_taxa
    log_file = snakemake.log[0]
    assembly_summary = snakemake.input.assembly_summary
    kraken_taxa_genome = snakemake.input.kraken_taxa_genome
    process_check = snakemake.output.process_check
    threads = snakemake.threads

    log(f"{genome_id}: Processing genome/proteome ...", log_file)


    # Read assembly summary to find URL
    url = None
    alternate_accession = None
    with open(assembly_summary, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == genome_accession or fields[17] == genome_accession:
                url = fields[19]
                pot_alt_acc = fields[17]
                if pot_alt_acc.startswith('GCF') or pot_alt_acc.startswith('GCA'):
                    alternate_accession = fields[17]
                break

    if not url:
        log(f"{genome_id}: Unable to find URL in assembly summary. Exiting.", log_file)
        sys.exit(1)

    # Define the proteome fasta file
    proteome_fasta = os.path.join(taxon_base_path, f"{Path(url).name}_protein.faa.gz")

    # Determine if accession is GCA or GCF and prefer GCF
    accession_to_use = genome_accession
    url_to_use = url
    switched_accession = False
    if genome_accession.startswith("GCA"):
        log(f"{genome_id}: Original accession is GCA. Checking for corresponding GCF accession...", log_file)

        if alternate_accession and alternate_accession.startswith("GCF"):
            log(f"{genome_id}: Found corresponding GCF accession: {alternate_accession}. Updating URL and accession.", log_file)
            accession_to_use = alternate_accession
            # Update URL to GCF
            url_to_use = url_to_use.replace(genome_accession, accession_to_use).replace("GCA", "GCF")
            switched_accession = True
        else:
            log(f"{genome_id}: No corresponding GCF accession found. Proceeding with GCA.", log_file)

    # Now, proceed to download the required archives: genomic.fna.gz, protein.faa.gz, genomic.gff.gz
    archives = ['genomic.fna.gz', 'protein.faa.gz', 'genomic.gff.gz']
    for archive in archives:
        filename = f"{Path(url_to_use).name}_{archive}"
        filepath = os.path.join(taxon_base_path, filename)

        # Check if the archive exists and is valid
        if not os.path.exists(filepath) or not is_gzip_valid(filepath):
            # Remove invalid file
            if os.path.exists(filepath):
                os.remove(filepath)
            
            log(f"{genome_id}: {archive} is missing or invalid. Attempting to download.", log_file)
            success = download_file(f"{url_to_use}/{filename}", filepath, log_file)

            if not success or not is_gzip_valid(filepath):
                if switched_accession:
                    try_accession = accession_to_use
                    log(f"{genome_id}: Failed to download or validate {archive} from GCF. Attempting GCA.", log_file)

                    # Attempt to download from original GCA accession
                    alt_filename = f"{Path(url).name}_{archive}"
                    alt_filepath = os.path.join(taxon_base_path, alt_filename)

                    success = download_file(f"{url}/{alt_filename}", alt_filepath, log_file)
                elif alternate_accession:
                    try_accession = alternate_accession
                    log(f"{genome_id}: Failed to download or validate {archive}. Attempting alternate accession {alternate_accession}.", log_file)

                    alternate_url = url_to_use.replace(accession_to_use, alternate_accession)

                    # Check if alternate accession is GCF or GCA
                    if alternate_accession.startswith("GCF"):
                        alternat_url = alternate_url.replace("GCA", "GCF")
                    else:
                        alternate_url = alternate_url.replace("GCF", "GCA")

                    alt_filename = f"{Path(alternate_url).name}_{archive}"
                    alt_filepath = os.path.join(taxon_base_path, alt_filename)

                    success = download_file(f"{alternate_url}/{alt_filename}", alt_filepath, log_file)
                else:
                    log(f"{genome_id}: No alternate accession found. Skipping download.", log_file)
                    success = False
                    continue

                if success and is_gzip_valid(alt_filepath):
                    log(f"{genome_id}: Successfully downloaded and validated {archive} from accession {try_accession}.", log_file)
                    # If GCF was found and downloaded, delete GCA (not the genome though, we do that at the end)
                    if alt_filename.startswith("GCF") and genome_accession.startswith("GCA") and archive != "genomic.fna.gz":
                        # Remove GCA archive
                        gca_filename = f"{Path(url).name}_{archive}"
                        gca_filepath = os.path.join(taxon_base_path, gca_filename)
                        if os.path.exists(gca_filepath):
                            os.remove(gca_filepath)
                            log(f"{genome_id}: Removed GCA archive {gca_filepath}.", log_file)     

                    # If archive is proteome, update proteome_fasta
                    if archive == 'protein.faa.gz':
                        proteome_fasta = filepath             
                else:
                    log(f"{genome_id}: Failed to download or validate {archive} from alternate accession.", log_file)
                    # Remove invalid files
                    if os.path.exists(alt_filepath):
                        os.remove(alt_filepath)
            else:
                log(f"{genome_id}: Successfully downloaded and validated {archive}.", log_file)
                # If GCF was found and downloaded, delete GCA (not the genome though, we do that at the end)
                if filename.startswith("GCF") and genome_accession.startswith("GCA") and archive != "genomic.fna.gz":
                    # Remove GCA archive
                    gca_filename = f"{Path(url).name}_{archive}"
                    gca_filepath = os.path.join(taxon_base_path, gca_filename)
                    if os.path.exists(gca_filepath):
                        os.remove(gca_filepath)
                        log(f"{genome_id}: Removed GCA archive {gca_filepath}.", log_file)

                # If archive is proteome, update proteome_fasta
                if archive == 'protein.faa.gz':
                    proteome_fasta = filepath
        else:
            log(f"{genome_id}: {archive} is already valid.", log_file)

    # After downloading, validate proteome fasta
    log(f"{genome_id}: Checking proteome fasta ...", log_file)

    proteome_valid = False
    # Check if proteome fasta exists and isnt empty
    if os.path.exists(proteome_fasta) and os.path.getsize(proteome_fasta) > 0:
        try:
            # Using subprocess to call seqkit for validation
            cmd = f"seqkit stats --seq-type protein {proteome_fasta}"
            result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0 and 'invalid FASTA/Q format' not in result.stdout:
                proteome_valid = True
        except Exception as e:
            log(f"{genome_id}: Exception during proteome validation: {e}", log_file)

    if proteome_valid:
        log(f"{genome_id}: Valid proteome fasta. No gene prediction needed. Finishing up.", log_file)
    else:
        # Remove invalid proteome fasta
        if os.path.exists(proteome_fasta):
            os.remove(proteome_fasta)
            log(f"{genome_id}: Invalid proteome fasta. Deleting proteome file.", log_file)

        # Check if genome is a eukaryote
        log(f"{genome_id}: Checking if genome is a eukaryote ...", log_file)
        is_eukaryote = False
        with open(assembly_summary, 'r') as f:
            for line in f:
                if genome_id in line and re.search(unsupported_taxa, line):
                    is_eukaryote = True
                    break

        if is_eukaryote:
            log(f"{genome_id}: Eukaryote without valid proteome fasta. Finishing up.", log_file)
        else:
            # If proteome fasta is invalid or missing and genome is not eukaryote, run gene prediction
            log(f"{genome_id}: Running gene prediction with Pyrodigal ...", log_file)
            try:
                gene_finder = pyrodigal.GeneFinder(meta=True, mask=True, min_mask=0)
                # Read the input genomic FASTA file
                sequences = {}
                with gzip.open(kraken_taxa_genome, 'rt') as f:
                    current_seq_id = None
                    current_seq = []
                    for line in f:
                        if line.startswith(">"):
                            if current_seq_id is not None:
                                sequences[current_seq_id] = ''.join(current_seq).rstrip('N')
                            current_seq_id = line.strip().lstrip('>').replace(' ', '_')
                            current_seq = []
                        else:
                            current_seq.append(line.strip())
                    if current_seq_id is not None:
                        sequences[current_seq_id] = ''.join(current_seq).rstrip('N')

                # Run Pyrodigal gene prediction
                predicted_genes = {}
                with multiprocessing.Pool(processes=threads) as pool:
                    args = [(gene_finder, seq_id, seq) for seq_id, seq in sequences.items()]
                    results = pool.map(run_gene_finder_on_sequence, args)
                    for seq_id, genes in results:
                        predicted_genes[seq_id] = genes

                # Write the predicted proteins to the output file
                with gzip.open(proteome_fasta, 'wt') as f:
                    for seq_id, genes in predicted_genes.items():
                        genes.write_translations(f, sequence_id=seq_id, include_stop=False)

                # Optionally write GFF file
                gff_file = proteome_fasta.replace(".faa.gz", ".gff.gz")
                with gzip.open(gff_file, 'wt') as f:
                    for seq_id, genes in predicted_genes.items():
                        genes.write_gff(f, sequence_id=seq_id, header=False)

                log(f"{genome_id}: Gene prediction complete. Proteome fasta and GFF files written.", log_file)

            except Exception as e:
                log(f"{genome_id}: Exception during gene prediction: {e}", log_file)
                sys.exit(1)

    # Remove genome fastas for genomes accesions using glob
    for accession in [genome_accession, alternate_accession]:
        genome_file = glob.glob(f"{taxon_base_path}/{accession}*_genomic.fna.gz")
        if genome_file:
            os.remove(genome_file[0])
            log(f"{genome_id}: Removed {genome_file[0]}.", log_file)

    # Create the process check file
    Path(process_check).touch()
    log(f"{genome_id}: Processing complete.", log_file)

if __name__ == "__main__":
    # Snakemake provides the 'snakemake' object automatically when using the 'script:' directive
    process_genome()
