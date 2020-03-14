#!/bin/bash -l

ass=mgmt
input0="Assembly/$ass.assembly.merged.fa"

output0=Analysis/amr
output1="Analysis/amr.done"

__conda_setup="$('/work/projects/ecosystem_biology/local_tools/miniconda/miniconda2/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"

conda activate abricate

mkdir -p $output0
for db in resfinder card argannot ecoh vfdb plasmidfinder ecoli_vf ncbi
do
    echo "# ABRicate $db start"
    abricate --db $db $input0 >Analysis/amr/abricate.$db.tsv
    echo "# ABricate $db done"
done

touch $output1

conda deactivate

