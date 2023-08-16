#!/usr/bin/zsh
for MODE in all domination no-necessary no-domination simple none; do
#for MODE in no-domination simple none; do
time ./experiment -z -s forts --early-stop -b 2 -u -r $MODE -n 5 --timeout=5400 -o results/all-zi/red_${MODE}.csv -w results/all-zi/sol $(cat cases.txt)
    date
    time ./experiment -s forts --early-stop -b 2 -u -r $MODE -n 5 --timeout=5400 -o results/normal/red_${MODE}.csv -w results/normal/sol $(cat cases.txt)
    date
done
