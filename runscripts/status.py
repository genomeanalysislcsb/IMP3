#!/usr/bin/env python

import subprocess
import sys
import time

RETRY_COUNT = 5  # Number of retries
RETRY_DELAY = 10  # Delay (in seconds) between retries

def get_sacct_output(jobid, cluster):
    # Verify clujster is a single string with no spaces and jobid is a number
    if (cluster and not cluster.replace("-", "").replace("_", "").isalnum()) or not jobid.isdigit():
        print(f"WARNING: Unexpected format in job info: {jobid} {cluster}", file=sys.stderr)

    # with Slurm multi-cluster the job id is only unique per cluster; on a single cluster there is no cluster part
    cmd_sacct = ["sacct", "--jobs", jobid, "--format", "State", "--noheader"]
    if cluster:
        cmd_sacct[1:1] = ["--clusters", cluster]
    try:
        output_sacct = subprocess.check_output(cmd_sacct).decode('utf-8').strip()
    except subprocess.CalledProcessError:
        print(f"Error: Failed to run sacct command: {' '.join(cmd_sacct)}", file=sys.stderr)
        print(f"Error message: {sys.exc_info()}", file=sys.stderr)
        return None
    return output_sacct

def main():
    if len(sys.argv) != 2:
        print(f"Error: Expected single argument with job info, but got: {sys.argv[1]}")
        sys.exit(1)

    # Extract job info from the provided argument
    jobinfo = sys.argv[1].strip()
    if ';' in jobinfo:
        jobinfo_parts = jobinfo.split(';')
        if len(jobinfo_parts) == 2:
            jobid, cluster = jobinfo_parts
        elif len(jobinfo_parts) == 1:
            jobid = jobinfo_parts[0]
            cluster = ''
        else:
            print(f"Error: Unexpected format in job info: {jobinfo}")
            sys.exit(1)
    elif jobinfo.startswith("Submitted batch job"):
        jobinfo_parts = jobinfo.split(' ')
        if len(jobinfo_parts) == 4:
            jobid = jobinfo_parts[3]
            cluster = ''
        elif len(jobinfo_parts) == 7:
            jobid = jobinfo_parts[3]
            cluster = jobinfo_parts[6]
        else:
            print(f"Error: Unexpected format in job info: {jobinfo}")
            sys.exit(1)
    else:
        print(f"Error: Expected job info, but got: {jobinfo}")
        sys.exit(1)

    # Retry logic for sacct command
    for _ in range(RETRY_COUNT):
        output_status = get_sacct_output(jobid, cluster)
        if output_status:
            break
        time.sleep(RETRY_DELAY)
    else:
        # All retries exhausted and still no output
        print(f"No output from sacct for job {jobid} on cluster {cluster} after {RETRY_COUNT} retries.", file=sys.stderr)
        print(f"Trying with seff command for job {jobid} on cluster {cluster}.", file=sys.stderr)
        # Try with `seff` as final measure
        cmd_seff = ["seff", jobid]
        try:
            output_seff = subprocess.check_output(cmd_seff).decode('utf-8').strip()
            # Get line 4 for status
            output_status = output_seff.splitlines()[3].split()[1]
        except subprocess.CalledProcessError:
            print(f"Error: Failed to run seff command: {' '.join(cmd_seff)}", file=sys.stderr)
            print(f"Error message: {sys.exc_info()}", file=sys.stderr)
            print("unknown")
            sys.exit(0)

    # Process the output
    output = output_status.splitlines()[0].split()[0]

    # Determine job status based on the output
    running_status = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
    if "COMPLETED" in output:
        print("success")
    elif any(r in output for r in running_status):
        print("running")
    else:
        print("failed")

if __name__ == '__main__':
    main()
