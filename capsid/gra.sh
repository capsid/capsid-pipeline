#!/bin/bash
# Ivan Borozan 
# Mar-13-2013
# Calculate genome relative abundence using 'stats' from readscan

xeno=$1
human=$2
out=$3

#cat $xeno | sort -t$'\t' +0 -1 -T $temp > $xeno.sorted
#cat $human | sort -t$'\t' +0 -1 -T $temp | uniq > $human.sorted

cat $human.sorted | uniq > $human.sorted.uniq

join -v1 -t$'\t' -1 1 -2 1 $xeno.sorted  $human.sorted.uniq  > $out/pathogen.sam

# rm $xeno.sorted $human.sorted $human.sorted.uniq

gzip $out/pathogen.sam

echo $(date +%T)" Calculating genome relative abundence..."

readscan.pl stats -R $READSCAN_PATHOGEN_REF -T $READSCAN_TAXON $out/pathogen.sam.gz > $out/pathogen.gra.txt

echo $(date +%T)" Genome relative abundence finished..."
