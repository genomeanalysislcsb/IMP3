#!/bin/bash

for i in `ls ../conda/*/etc/conda/activate.d/activate-r-base.sh`
do
	cp $i $i.org
	sed 's/^R CMD javareconf/#R CMD javareconf/' $i.org >$i
done
