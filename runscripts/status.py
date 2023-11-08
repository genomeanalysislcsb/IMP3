#!/usr/bin/env python

import subprocess
import sys
import time

RETRY_COUNT = 3  # Number of retries
RETRY_DELAY = 5  # Delay (in seconds) between retries

def get_sacct_output(jobid, cluster):
    cmd_sacct = ["sacct", "--allclusters", "--clusters", cluster, "--jobs", jobid, "--format", "State", "--noheader"]
    try:
        output_sacct = subprocess.check_output(cmd_sacct).decode('utf-8').strip()
    except subprocess.CalledProcessError:
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
        output_sacct = get_sacct_output(jobid, cluster)
        if output_sacct:
            break
        time.sleep(RETRY_DELAY)
    else:
        # All retries exhausted and still no output
        print(f"No output from sacct for job {jobid} on cluster {cluster} after {RETRY_COUNT} retries.", file=sys.stderr)
        print("unknown")
        sys.exit(0)

    # Process the output
    output = output_sacct.splitlines()[0].split()[0]

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
