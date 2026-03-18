#!/bin/bash

squeue_result=$(squeue --me -o '%.18i %.9f %.20j %.4t %.9M %.4D %.8a %.12R' --qos regular_1)
echo "$squeue_result" | tail -n +2 | grep " R " | wc -l | tr -d '\n' &&
    echo -n " / " &&
    echo $(echo "$squeue_result" | wc -l)-1 | bc
echo "$squeue_result" | tail -n $(echo $LINES-10 | bc)
