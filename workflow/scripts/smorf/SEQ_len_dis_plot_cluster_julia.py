import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import re

def glob_to_regex(glob_pattern):
    """
    Convert a glob pattern (with *) to a regex pattern.
    Only supports '*' wildcard.
    """
    # Escape all regex special chars except '*'
    regex = re.escape(glob_pattern)
    regex = regex.replace(r'\*', '.*')
    match = re.search(r'/([^/]*\*[^/]*)/', glob_pattern)
    if match:
        sample_pattern = match.group(1)
        regex = sample_pattern.replace('*', '.*')
        return '^' + regex + '$'
    else:
        return '.*'

def get_group_patterns_and_samples(filein):
    """
    Returns both regex patterns and direct sample mappings.
    Supports both formats:
    - regexp,group_name,pattern
    - sample_name,group_name,path (optional)
    """
    group_patterns = []
    sample_groups = {}
    
    with open(filein, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
                
            if parts[0] == "regexp":
                # Format: regexp,group_name,pattern
                group = parts[1]
                pattern = parts[2]
                regex = glob_to_regex(pattern)
                group_patterns.append((group, re.compile(regex)))
            else:
                # Format: sample_name,group_name[,path]
                sample_name = parts[0]
                group_name = parts[1]
                sample_groups[sample_name] = group_name
    
    return group_patterns, sample_groups

def sample_to_group(sample, group_patterns, sample_groups):
    """
    Returns the group for a sample, checking direct mapping first, then regex patterns.
    """
    if sample is None:
        return None
    
    # First check direct sample mapping
    if sample in sample_groups:
        return sample_groups[sample]
    
    # Then check regex patterns
    for group, regex in group_patterns:
        if regex.match(sample):
            return group
    
    return None

def extract_data_from_fasta(fasta_path):
    results = []
    with open(fasta_path, 'r') as f:
        current_sample = None
        current_smorf = None
        current_length = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_smorf is not None:
                    results.append([current_sample, current_smorf, current_length])
                header = line[1:]
                if "::" in header:
                    sample_part = header.split("::", 1)[0]
                    current_sample = sample_part.strip()
                    current_smorf = header.split()[0]
                else:
                    current_sample = None
                    current_smorf = header.split()[0]
                current_length = 0
            else:
                current_length += len(line)
        if current_smorf is not None:
            results.append([current_sample, current_smorf, current_length])
    return pd.DataFrame(results, columns=["Sample", "smORF", "Length"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot smORF sequence length distribution with groups.')
    parser.add_argument('-i', type=str, required=True, help='Path to input FASTA file')
    parser.add_argument('--group_info', type=str, required=True, help='Path to group mapping file')
    parser.add_argument('-o', type=str, required=True, help='Output file path for plot')
    args = parser.parse_args()

    # Get both group regex patterns and direct sample mappings
    group_patterns, sample_groups = get_group_patterns_and_samples(args.group_info)
    print(f"Found {len(group_patterns)} regex patterns and {len(sample_groups)} direct sample mappings")
    
    # Extract data from single FASTA file
    df = extract_data_from_fasta(args.i)
    print(f"Extracted {len(df)} sequences from FASTA")
    
    # Debug: Check sample names before grouping
    print("Sample names found:", df["Sample"].unique()[:10])
    print("Direct mappings:", sample_groups)
    print("Regex patterns:", [(g, p.pattern) for g, p in group_patterns])
    
    # Map sample names to groups using both methods
    df["Group"] = df["Sample"].apply(lambda s: sample_to_group(s, group_patterns, sample_groups))
    
    # Debug: Check group assignment
    print("Groups assigned:", df["Group"].value_counts())
    
    # Filter out rows without group assignment
    df_filtered = df.dropna(subset=["Group"])
    print(f"After filtering: {len(df_filtered)} sequences with groups")
    
    # Check if we have any data left
    if len(df_filtered) == 0:
        print("ERROR: No sequences matched any group patterns or direct mappings!")
        print("Available sample names:", df["Sample"].unique())
        print("Direct mappings:", sample_groups)
        print("Regex patterns:", [(g, p.pattern) for g, p in group_patterns])
        exit(1)
    
    # Rest of the plotting code remains the same...
    # Check if we have valid length data
    if df_filtered["Length"].isna().all():
        print("ERROR: All length values are NaN!")
        exit(1)
    
    # Calculate range
    min_len = int(df_filtered["Length"].min())
    max_len = int(df_filtered["Length"].max())
    print(f"Length range: {min_len} to {max_len}")
    
    # Create plot
    plt.figure(figsize=(14, 6))
    
    sns.histplot(
        data=df_filtered,
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
    
    # Only set custom ticks if range is reasonable
    if max_len - min_len <= 50:
        plt.xticks(ticks=range(min_len, max_len + 1), rotation=90)
    
    plt.xlim(min_len - 0.5, max_len + 0.5)
    plt.tight_layout()

    # Save plot
    plt.savefig(args.o, dpi=300)
    print(f"Plot saved to {args.o}")

