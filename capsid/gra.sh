#!/bin/bash
# Calculate genome relative abundence using 'stats' from readscan

xeno=$1
human=$2
out=$3

directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#cat $xeno | sort -t$'\t' +0 -1 -T $temp > $xeno.sorted
#cat $human | sort -t$'\t' +0 -1 -T $temp | uniq > $human.sorted

cat $human | uniq > $human.uniq

join -v1 -t$'\t' -1 1 -2 1 $xeno  $human.uniq  > $out/pathogen.sam

# rm $xeno.sorted $human.sorted $human.sorted.uniq

gzip -f $out/pathogen.sam

echo $(date +%T)" Calculating genome relative abundance..."

if [[ ! $READSCAN_PATHOGEN_REF ]]; then
	echo "Missing environment variable: READSCAN_PATHOGEN_REF"
	exit 1
fi

if [[ ! $READSCAN_TAXON ]]; then
	echo "Missing environment variable: READSCAN_TAXON"
	exit 1
fi

# We have a dependency here, and we don't have an obvious way to deploy it. So let's 
# bundle it into a local lib and add that to @INC on the command line. 

echo $(date +%T)" perl -I "${directory}/lib" "${directory}/readscan.pl" stats -R $READSCAN_PATHOGEN_REF -T $READSCAN_TAXON $out/pathogen.sam.gz"
perl -I "${directory}/lib" "${directory}/readscan.pl" stats --data -R $READSCAN_PATHOGEN_REF -T $READSCAN_TAXON $out/pathogen.sam.gz > $out/pathogen.gra.txt

echo $(date +%T)" Genome relative abundance finished..."
