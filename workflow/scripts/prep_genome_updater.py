import os
import logging


### FUNCTIONS ###
def setup_logger(logfile=None, level=logging.DEBUG):
    """
    Setup a logger
    """

    # Formatting
    formatter = logging.Formatter('%(asctime)s - %(message)s', '%d/%m/%Y %I:%M:%S %p')

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Console handler
    ch = logging.StreamHandler()
    logger.addHandler(ch)

    # File handler
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return

def find_genome_updater():
    '''Function to find genome_updater.sh in the conda env path.'''
    conda_env_bin_dir = os.path.join(os.environ['CONDA_PREFIX'], 'bin')
    for root, dirs, files in os.walk(conda_env_bin_dir):
        for file in files:
            if file == 'genome_updater.sh':
                return os.path.join(root, file)
    return None
###

### VARS ###
# Define the original and replacement functions
original_function = '''filter_taxids_ncbi() # parameter: ${1} assembly_summary file, ${2} ncbi_tax file - return number of lines
{
    # Keep only selected taxid lineage, removing at the end duplicated entries from duplicates on taxids
    tmp_lineage=$(tmp_file "lineage.tmp")
    for tx in ${taxids//,/ }; do
        txids_lin=$(grep "[^0-9]${tx}[^0-9]" "${2}" | cut -f 1) #get only taxids in the lineage section
        echolog " - $(count_lines "${txids_lin}") children taxids in the lineage of ${tx}" "0"
        echo "${txids_lin}" >> "${tmp_lineage}" 
    done
    lineage_taxids=$(sort ${tmp_lineage} | uniq | tr '\\n' ',')${taxids} # put lineage back into the taxids variable with the provided taxids
    rm "${tmp_lineage}"

    # Join with assembly_summary based on taxid field 6
    join -1 6 -2 1 <(sort -k 6,6 "${1}") <(echo "${lineage_taxids//,/$'\\n'}" | sort -k 1,1) -t$'\\t' -o ${join_as_fields1} | sort | uniq > "${1}_taxids"
    mv "${1}_taxids" "${1}"
    count_lines_file "${1}"
}
'''

replacement_function = '''filter_taxids_ncbi() { # parameter: ${1} assembly_summary file, ${2} ncbi_tax file - return number of lines
    tmp_taxids=$(tmp_file "taxids.tmp")
    
    # Write the provided taxids directly to the temporary file
    for tx in ${taxids//,/ }; do
        echo "${tx}" >> "${tmp_taxids}"
    done

    # Join with assembly_summary based on taxid field 6
    join -1 6 -2 1 <(sort -k 6,6 "${1}") <(sort -k 1,1 "${tmp_taxids}") -t$'\\t' -o ${join_as_fields1} | sort | uniq > "${1}_taxids"
    mv "${1}_taxids" "${1}"
    rm "${tmp_taxids}"
    
    count_lines_file "${1}"
}
'''
###


# Setup logger
log = snakemake.log[0]
setup_logger(log, level=logging.INFO)

# Find `genome_updater.sh` in conda environment dir
genome_updater_path = find_genome_updater()
if not genome_updater_path:
    raise FileNotFoundError("genome_updater.sh not found in the conda environment")

# Read the original script
logging.info(f"Reading genome_updater.sh from {genome_updater_path}")
with open(genome_updater_path, 'r') as file:
    script_content = file.read()

# Replace the first instance of the original function with the replacement function
logging.info("Replacing the original function with the replacement function")
updated_script_content = script_content.replace(original_function, replacement_function, 1)

# Write the modified script to a temporary path
logging.info("Writing the modified script to a temporary path")
tmp_path = genome_updater_path.replace('genome_updater.sh', 'genome_updater.sh.tmp')
with open(tmp_path, 'w') as file:
    file.write(updated_script_content)

# Move the modified script to the original path
logging.info("Moving the modified script to the original path")
os.rename(tmp_path, genome_updater_path)

# Create the output file to indicate completion
output_path = snakemake.output[0]
logging.info(f"Creating output file at {output_path}")
with open(output_path, 'w') as file:
    file.write("genome_updater ready")

# Make file executable
logging.info("Making the modified script executable")
os.chmod(genome_updater_path, 0o755)

logging.info("Done.")
