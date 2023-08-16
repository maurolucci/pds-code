#!/usr/bin/zsh
TIMEOUT=7200
NRUN=5
mkdir -p results/solver/{ours,bozeman-smith,jovanovic}{,-r}{,-z}/sol
for OTHER_ARGS in '' '-z'; do
    echo $OTHER_ARGS
    time ./experiment -s forts -u -r --early-stop -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/ours-r$OTHER_ARGS.csv" --fort-stats results/solver/ours-r$OTHER_ARGS-fort.csv -w results/solver/ours-r$OTHER_ARGS/sol $(cat cases.txt)
    date
    time ./experiment -s forts -u --early-stop -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/ours$OTHER_ARGS.csv" --fort-stats results/solver/ours$OTHER_ARGS-fort.csv -w results/solver/ours$OTHER_ARGS/sol $(cat cases.txt)
    date
    time ./experiment -s bozeman2 --fort-init=3 -u -r -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/bozeman-smith-r$OTHER_ARGS.csv" --fort-stats results/solver/bozeman-smith-r$OTHER_ARGS-fort.csv -w results/solver/bozeman-smith-r$OTHER_ARGS/sol $(cat cases.txt)
    date
    time ./experiment -s bozeman2 --fort-init=3 -u -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/bozeman-smith$OTHER_ARGS.csv" --fort-stats results/solver/bozeman-smith$OTHER_ARGS-fort.csv -w results/solver/bozeman-smith$OTHER_ARGS/sol $(cat cases.txt)
    date
    time ./experiment -s gurobi -u -r -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/jovanovic-r$OTHER_ARGS.csv" --fort-stats results/solver/jovanovic-r$OTHER_ARGS-fort.csv -w results/solver/jovanovic-r$OTHER_ARGS/sol $(cat cases.txt)
    date
    time ./experiment -s gurobi -u -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/jovanovic$OTHER_ARGS.csv" --fort-stats results/solver/jovanovic$OTHER_ARGS-fort.csv -w results/solver/jovanovic$OTHER_ARGS/sol $(cat cases.txt)
    date
done
