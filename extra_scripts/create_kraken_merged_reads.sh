

for i in sample; do
paste -d 'x' \
      /mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/imp3/M1.1-V1/Preprocessing/mt.r1.preprocessed.fq \
      /mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/imp3/M1.1-V1/Preprocessing/mt.r2.preprocessed.fq | less -RS