#!/usr/bin/zsh
TIMEOUT=36000
NRUNS=5
for INSTANCE in inputs/{case9241pegase.graphml,psd-{Texas,Western,Eastern,USA}.graphml}; do
    ./experiment -r -s forts -n$NRUNS -v -t $TIMEOUT -b 1 --fort-init=3 --early-stop=2 --write-bounds=results/bounds-ours --fort-stats=results/forts-ours -w results/bounds-ours/sol -o results/bounds-ours/$INSTANCE.csv -- $INSTANCE
done
for INSTANCE in inputs/{case9241pegase.graphml,psd-{Texas,Western,Eastern,USA}.graphml}; do
    ./experiment -r -s bozeman2 -n$NRUNS -v -t $TIMEOUT -b 1 --fort-init=3 --write-bounds=results/bounds-smith -o results/bounds-smith/$INSTANCE.csv --fort-stats=results/forts-smith -- $INSTANCE 
done
for INSTANCE in inputs/{case9241pegase.graphml,psd-{Texas,Western,Eastern,USA}.graphml}; do
    ./experiment -r -s jovanovic -n$NRUNS -v -t $TIMEOUT -b 1 --early-stop=2 --write-bounds=results/bounds-milp -o results/bounds-milp/$INSTANCE.csv -- $INSTANCE
done

