#! /bin/bash
NGEN="100000"
for TAXA in "mammals" "birds" "squamates" "fish" "amphibians"
do
    (Rscript --no-save --no-restore --verbose runBetaBayou.R > ../output/runs/log.$TAXA.Rout 2>&1 $TAXA $NGEN) &

done
