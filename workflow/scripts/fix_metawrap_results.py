#!/usr/bin/env python3

import os
import sys
import glob

mw_stats = sys.argv[1]
out_file = sys.argv[2]

if os.stat(mw_stats).st_size > 0:
    stats_data = {}
    got_grid_header = False
    got_gtdbtk_header = False
    # Get stats
    print('Getting MetaWrap bin stats.')
    with open(mw_stats,'r') as mws:
        # Get mw stats header
        header = next(mws).strip('\n \t')
        for line in mws:
            line = line.strip('\n \t').split('\t')
            bin = line[0]
            stats_data[bin] = line
    print("Reading GRiD results.")
    # Get all GTBtk stat files
    grid_out = glob.glob("Binning/MetaWrap/bins/*/grid/reads.GRiD.txt")
    if grid_out:
        for bin_grid_data in grid_out:
            bin = bin_grid_data.split('/')[-3]
            with open(bin_grid_data,'r') as f:
                if not got_grid_header:
                    header = header + '\t' + '\t'.join(next(f).strip(' \t\n').split('\t')[1:])
                    got_grid_header = True
                else:
                    next(f)
                stats_data[bin].extend(next(f).strip(' \t\n').split('\t')[1:])

    print("Reading GTDBtk results.")
    with open(out_file,'w') as of:
        # Get all GTBtk stat files
        gtdbtk_out = glob.glob("Binning/MetaWrap/bins/*/GTDB/classify/gtdbtk.*.summary.tsv")
        if gtdbtk_out:
            for bin_gtdbtk_data in gtdbtk_out:
                bin = bin_gtdbtk_data.split('/')[-4]
                with open(bin_gtdbtk_data,'r') as f:
                    if not got_gtdbtk_header:
                        header = header + '\t' + '\t'.join(next(f).strip(' \t').split('\t')[1:])
                        of.write(header)
                        got_grid_header = True
                    else:
                        next(f)
                    of.write('\t'.join(stats_data[bin] + next(f).split('\t')[1:]))
        else:
            of.write(header)
            for bin in stats_data.keys():
                of.write('\t'.join(stats_data[bin] + ['\n']))
else:
    print("No output from MetaWrap!")
    os.mknod(output)


