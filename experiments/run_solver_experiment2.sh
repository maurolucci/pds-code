#!/usr/bin/zsh
TIMEOUT=3600
NRUN=5
mkdir -p results/solver/{blasius,bozeman-smith,jovanovic,brimkov,brimkov-geq,brimkov-alt}{,-z}/sol
for OTHER_ARGS in '' '-z'; do
    echo $OTHER_ARGS
    if [ "$OTHER_ARGS" = "-z" ]; then   
        CASES="cases.txt" 
    fi
    if [ "$OTHER_ARGS" = "" ]; then   
        CASES="cases2000.txt" 
    fi
    time ./experiment -s forts -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/blasius$OTHER_ARGS.csv" --fort-stats results/solver/blasius$OTHER_ARGS-fort.csv -w results/solver/blasius$OTHER_ARGS/sol $(cat $CASES)
    date
    time ./experiment -s bozeman2 --fort-init=3 -u -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/bozeman-smith$OTHER_ARGS.csv" --fort-stats results/solver/bozeman-smith$OTHER_ARGS-fort.csv -w results/solver/bozeman-smith$OTHER_ARGS/sol $(cat $CASES)
    date
    time ./experiment -s jovanovic -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/jovanovic$OTHER_ARGS.csv" --fort-stats results/solver/jovanovic$OTHER_ARGS-fort.csv -w results/solver/jovanovic$OTHER_ARGS/sol $(cat $CASES)
    date
    time ./experiment -s brimkov -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/brimkov$OTHER_ARGS.csv" --fort-stats results/solver/brimkov$OTHER_ARGS-fort.csv -w results/solver/brimkov$OTHER_ARGS/sol $(cat $CASES)
    date
    time ./experiment -s brimkov-geq -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/brimkov-geq$OTHER_ARGS.csv" --fort-stats results/solver/brimkov-geq$OTHER_ARGS-fort.csv -w results/solver/brimkov-geq$OTHER_ARGS/sol $(cat $CASES)
    date
    time ./experiment -s brimkov-alt -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/brimkov-alt$OTHER_ARGS.csv" --fort-stats results/solver/brimkov-alt$OTHER_ARGS-fort.csv -w results/solver/brimkov-alt$OTHER_ARGS/sol $(cat $CASES)
    date
done
