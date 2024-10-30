#!/usr/bin/zsh
TIMEOUT=3600
NRUN=5
for OTHER_ARGS in '-z'; do
    echo $OTHER_ARGS
    for SOLVER in 'paths' 'paths-forts'; do
        echo $SOLVER
        if [ "$SOLVER" = "paths" ]; then   
            CASES="cases-z-7000.txt" 
        fi
        if [ "$SOLVER" = "paths-forts" ]; then   
            CASES="cases-z.txt" 
        fi
        for INIT in '' '11' '12' '13' '21' '22' '31' '32' '33' '1221' '1321'; do
            time ./experiment -s $SOLVER$INIT -n $NRUN $OTHER_ARGS --timeout=$TIMEOUT -o "results/solver/$SOLVER$INIT$OTHER_ARGS.csv" --fort-stats results/solver/$SOLVER$INIT$OTHER_ARGS-fort.csv -w results/solver/$SOLVER$INIT$OTHER_ARGS/sol $(cat $CASES)
            date
        done
    done
done
