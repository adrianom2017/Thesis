#!/bin/sh
#
bootstrap=1000
adjust=0
echo "NEM for network inference, with bootstrap = ${bootstrap}"
#
for i in {1,3}
do
echo "Network ${i}"
bsub -W 05:00 -n 1 -R "rusage[mem=2000]" "Rscript --vanilla ./DREAM5_NEM.R $i $bootstrap $adjust > output/outputNEM${adjust}${i}"
echo "NEM submitted for network ${i}"
done
mv lsf.* lsf/
