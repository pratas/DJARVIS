#!/bin/bash
#
# DB.fa top.csv NEAN.fq
#
READS_FILE="/home/pratas/DATABASE/NEANDERTHAL/FASTQ/NEAN.fq";
DB_FILE="/home/pratas/DATABASE/DB/DB.fa";
TOP_FILE="top_nean.csv";
THRESHOLD="5.0";
PID=$(echo $$); # CURRENT PROCESS
KESTREL_ARGS=" -v -n 10 -F -m 14:1:1:0/0 -g 0.9 -t 0.7 ";
#
#
#
cat $TOP_FILE \
| awk '{ if($3 >= "'$THRESHOLD'" && $3 <= 100) print $1"\t"$2"\t"$3"\t"$4; }' \
| awk '{ print $4;}' \
| tr '|' '\t' \
| awk '{ print $2;}' \
> GIS-$PID;
#
#
#
cat GIS-$PID | while read line
  do
  echo "Running $line";
  namex=`echo $line | tr ' ' '_'`;
  ./goose-extractreadbypattern $line < $DB_FILE > $namex.fa;
  ./KESTREL $KESTREL_ARGS -o filt-$namex.fq $namex.fa $READS_FILE;
  ./bowtie-build $namex.fa index-$namex;
  ./bowtie -a -v 3 --sam index-$namex filt-$namex.fq | ./samtools view -bh > aligned-$namex.bam;
  mapDamage -i aligned-$namex.bam -r $namex.fa;
  done
#
#
#
