#!/bin/bash

for i in `seq 1 $1`;
do
    nohup ./estimation_stdev input.txt out_plus$i.txt $i &
done

for i in `seq 1 $1`;
do
    nohup ./estimation_stdev input.txt out_minus$i.txt -$i &
done
