# smORF module of IMP3. From the smORF prediction pipeline of Luis F. Delgado and Julia Ortis Sunyer
# (gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper; Davin, Ortis Sunyer, Delgado et al.,
# Nat Commun 2026, doi:10.1038/s41467-026-72762-5). GPL-3.0 (LICENSE.GPL-3.0 in this folder).
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import glob

def get_info(filein):
    samples = dict()
    with open(filein, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if parts[0] == "regexp":
                group = parts[1]
                pattern = parts[2]
                # Use glob to expand the pattern
                matched_files = glob.glob(pattern)
                for filepath in matched_files:
                    # Extract sample name from the path (customize as needed)
                    # Example: sample name is the parent directory of 'Assembly'
                    path_parts = filepath.split(os.sep)
                    try:
                        idx = path_parts.index("Assembly")
                        sample_name = path_parts[idx - 1]
                    except ValueError:
                        # Fallback: use filename or another method
                        sample_name = os.path.basename(filepath).split(".")[0]
                    samples[sample_name] = group
            else:
                # Standard line: sample_name, group, path
                sample_name = parts[0]
                group = parts[1]
                samples[sample_name] = group
    return samples

def extract_data_from_faa(path):
    results = []
    for file in glob.glob(os.path.join(path, "*.faa")):
        with open(file, 'r') as f:
            sample = os.path.basename(file).split(".")[0]
            smorf = None
            length = 0
            for line in f:
                if line.startswith('>'):
                    if smorf is not None:
                        results.append([sample, smorf, length])
                    handle = line[1:].split()[0]
                    smorf = handle.split("::")[1] if "::" in handle else handle
                    length = 0
                else:
                    length += len(line.strip())
            if smorf is not None:
                results.append([sample, smorf, length])
    df = pd.DataFrame(results, columns=["Sample", "smORF", "Length"])
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot smORF sequence length distribution with groups.')
    parser.add_argument('-i', type=str, required=True, help='Path to directory with .faa files')
    parser.add_argument('--group_info', type=str, required=True, help='Path to group mapping file')
    parser.add_argument('-o', type=str, required=True, help='Output file path for plot (e.g., output.png)')
    args = parser.parse_args()

    # Get group mapping and process data
    samples_to_groups = get_info(args.group_info)
    df = extract_data_from_faa(args.i)
    
    # Add group information and clean data
    df["Group"] = df["Sample"].map(samples_to_groups)
    df = df.dropna(subset=["Group"])

    # Create plot
    plt.figure(figsize=(14, 6))
    min_len = df["Length"].min()
    max_len = df["Length"].max()
    
    sns.histplot(
        data=df,
        x="Length",
        hue="Group",
        multiple="dodge",
        discrete=True,
        palette="viridis",
        edgecolor="black"
    )

    plt.title("smORF Length Distribution by Group")
    plt.xlabel("smORF Length (AA)")
    plt.ylabel("Number of smORFs")
    plt.xticks(ticks=range(min_len, max_len + 1), rotation=90)
    plt.xlim(min_len - 0.5, max_len + 0.5)
    plt.tight_layout()

    # Save and show plot
    plt.savefig(args.o, dpi=300)
#    plt.show()
