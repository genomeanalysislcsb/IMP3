#!/bin/bash

cd $IMPROOT/conda

for i in `ls */etc/conda/activate.d/activate-r-base.sh`
do
	cp $i $i.org
	sed 's/^R CMD javareconf/#R CMD javareconf/' $i.org >$i
done
