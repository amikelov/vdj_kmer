#!/bin/bash

mkdir -p $1/kmer_format
for i in $(ls $1)
    do paste <(paste <( tail -n +2 $i| cut -f33) <( tail -n +2 $i| cut -f 6)) <(tail -n +2 $i | cut -f3) |\
    sed 's/\*00\((.....\))//g' |sed 's/\*00\((...\))//g'|sed 's/\*00\((....\))//g'|sed 's/\*00\((..\))//g' |sed "s/^/${i}\t/" \
    > $1/$i.txt
done
