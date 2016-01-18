#! /bin/bash

echo "$1" > adjacency_matrix.txt

for (( i=0; i<$1; i++ ))
do
    for (( j=0; j<$1; j++ ))
    do
        echo -n "0 " >> adjacency_matrix.txt
    done
    echo >> adjacency_matrix.txt
done

