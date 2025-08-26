#!/usr/bin/zsh
TIMEOUT=3600
NRUN=5
for OTHER_ARGS in ''; do
    echo $OTHER_ARGS
    if [ "$OTHER_ARGS" = "" ]; then   
        CASES="cases.txt" 
    fi
    time ./experiment -s forts -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/blasius$OTHER_ARGS.csv" --forts-stats results/solver/blasius$OTHER_ARGS-fort.csv -w results/solver/blasius$OTHER_ARGS/sol $(cat $CASES)
    date
    time ./experiment -s efpss -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/efps$OTHER_ARGS.csv" --forts-stats results/solver/efps$OTHER_ARGS-fort.csv -w results/solver/efps$OTHER_ARGS/sol $(cat $CASES)
    date
done
