#! /bin/bash
for TAXA in "mammals" "birds" "squamates" "fish" "amphibians"
do
    (Rscript --no-save --no-restore --verbose runBayouTest.R > ../output/runs/outfile.$TAXA.Rout 2>&1 $TAXA) &

done
