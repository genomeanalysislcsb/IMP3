# smORF module of IMP3. From the smORF prediction pipeline of Luis F. Delgado and Julia Ortis Sunyer
# (gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper; Davin, Ortis Sunyer, Delgado et al.,
# Nat Commun 2026, doi:10.1038/s41467-026-72762-5). GPL-3.0 (LICENSE.GPL-3.0 in this folder).
import argparse
import os
import re
import fnmatch
from collections import defaultdict
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

'''
def get_info(filein):
    samples = defaultdict(list)
    with open(filein, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            group = parts[1]
            if line.startswith("regexp"):
                try:
                    pattern = next(p for p in parts[2].split("/") if "*" in p)
                except StopIteration:
                    raise ValueError("No wildcard pattern found in regexp line.")
                main_dir, rest_dir = parts[2].split(pattern, 1)
                main_dir = main_dir.rstrip("/")
                rest_dir = rest_dir.lstrip("/")
                cpattern = re.compile(fnmatch.translate(pattern))
                for f in os.listdir(main_dir):
                    if cpattern.fullmatch(f) and os.path.exists(os.path.join(main_dir, f, rest_dir)):
                        samples[f].append(group)
            else:
                sample_name = parts[0]
                samples[sample_name].append(group)
    return samples
'''

def get_info(filein):
    samples = defaultdict(list)
    with open(filein, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            group = parts[1]
            if line.startswith("regexp"):
                try:
                    # Find the pattern with wildcard
                    pattern = next(p for p in parts[2].split("/") if "*" in p)
                except StopIteration:
                    raise ValueError("No wildcard pattern found in regexp line.")
                
                # Split the path at the pattern
                path_parts = parts[2].split("/")
                pattern_index = next(i for i, part in enumerate(path_parts) if "*" in part)
                
                main_dir = "/".join(path_parts[:pattern_index])
                rest_dir = "/".join(path_parts[pattern_index + 1:])
                
                # Check if main_dir exists
                if not os.path.exists(main_dir):
                    print(f"Warning: Directory {main_dir} does not exist, skipping pattern {pattern}")
                    continue
                
                cpattern = re.compile(fnmatch.translate(pattern))
                
                try:
                    for f in os.listdir(main_dir):
                        if cpattern.fullmatch(f):
                            full_path = os.path.join(main_dir, f, rest_dir)
                            if os.path.exists(full_path):
                                samples[f].append(group)
                            else:
                                print(f"Warning: Path {full_path} does not exist")
                except OSError as e:
                    print(f"Error reading directory {main_dir}: {e}")
                    continue
            else:
                sample_name = parts[0]
                samples[sample_name].append(group)
    return samples

def parse_cluster_file(cluster_path):
    representative_to_samples = defaultdict(set)
    with open(cluster_path, 'r') as file:
        for line in file:
            representative, member = line.strip().split('\t')
            rep_sample = representative.split('::', 1)[0]
            mem_sample = member.split('::', 1)[0]
            representative_to_samples[representative].add(rep_sample)
            representative_to_samples[representative].add(mem_sample)
    return representative_to_samples


def process_single_representative(rep, samples, sample_to_group):
    group_counts = defaultdict(int)
    for sample in samples:
        for group in sample_to_group.get(sample, []):
            group_counts[group] += 1
    result = {"Representative": rep}
    result.update(group_counts)
    return result


def count_group_presence_parallel(representative_to_samples, sample_to_group, threads):
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(process_single_representative, rep, samples, sample_to_group)
            for rep, samples in representative_to_samples.items()
        ]
        for future in futures:
            results.append(future.result())
    return results


def main():
    parser = argparse.ArgumentParser(description="Count unique samples per group for each representative protein.")
    parser.add_argument("-i", required=True, help="Input cluster.tsv file")
    parser.add_argument("-g", required=True, help="Sample-to-group mapping file")
    parser.add_argument("-o", required=True, help="Output TSV file")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of CPUs to use (default: 4)")
    args = parser.parse_args()

    sample_to_group = get_info(args.g)
    representative_to_samples = parse_cluster_file(args.i)
    results = count_group_presence_parallel(representative_to_samples, sample_to_group, args.threads)

    df = pd.DataFrame(results)
    group_columns = [col for col in df.columns if col != "Representative"]
    df = df[["Representative"] + sorted(group_columns, key=str.lower)]
    df.to_csv(args.o, sep="\t", index=False)

    print(f"Results saved to {args.o}")


if __name__ == "__main__":
    main()
