#!/bin/bash
#
#
# $1 FASTA filename containning the reference
# $2 FASTQ reads
#
#
GET_FALCON=0;
GET_GULL=0;
GET_GOOSE=0;
GET_KESTREL=0;
GET_MAPDAMAGE=0;
GET_SAMTOOLS=0;
GET_BOWTIE=0;
#
#
RUN_FALCON=0;
RUN_GULL=0;
RUN_KESTREL=0;
RUN_BOWTIE=1;
RUN_MAPDAMAGE=1;
#
#
#==============================================================================
# GET FALCON
#
if [[ "$GET_FALCON" -eq "1" ]]; then
  rm -fr falcon FALCON FALCON-* *.pl
  git clone https://github.com/pratas/falcon.git
  cd falcon/src/
  cmake .
  make
  cp FALCON ../../
  cp FALCON-FILTER ../../
  cp FALCON-EYE ../../
  cd ../../
  cp falcon/scripts/*.pl .
fi
#==============================================================================
# GET GULL
#
if [[ "$GET_GULL" -eq "1" ]]; then
  rm -fr GULL/ GULL-map GULL-visual
  git clone https://github.com/pratas/GULL.git
  cd GULL/src/
  cmake .
  make
  cp GULL-map ../../
  cp GULL-visual ../../
  cd ../../
fi
#==============================================================================
# GET GOOSE
#
if [[ "$GET_GOOSE" -eq "1" ]]; then
  rm -fr goose/ goose-*
  git clone https://github.com/pratas/goose.git
  cd goose/src/
  make
  cp goose-* ../../
  cd ../../
fi
#==============================================================================
# GET KESTREL
#
if [[ "$GET_KESTREL" -eq "1" ]]; then
  rm -rf kestrel/ KESTREL
  git clone https://github.com/pratas/kestrel.git
  cd kestrel/src/
  cmake .
  make
  cp KESTREL ../../
  cd ../../
fi
#==============================================================================
# GET MAPDAMAGE 2
#
if [[ "$GET_MAPDAMAGE" -eq "1" ]]; then
  # 1. install mapDamage
  git clone https://github.com/ginolhac/mapDamage.gi
  cd mapDamage
  sudo python setup.py install
  # 2. install pysam
  # 2.1 cynthon
  mkdir cynthon
  cd cynthon
  wget https://pypi.python.org/packages/cb/0c/9d64e6ed68e76eb7d8a1ca8ec0d42ba0cac31dae10c6fbe6b3a5385b83a7/Cython-0.25.2-cp26-cp26mu-manylinux1_x86_64.whl#md5=81f24ce7637de392299d7cb9b12ec8cd
  pip install Cython --install-option="--no-cython-compile"
  cd ..
  # 2.2 pysam
  rm -f pysam-0.7.5.tar.gz
  wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/pysam/pysam-0.7.5.tar.gz
  tar -vxf pysam-0.7.5.tar.gz
  pysam-0.7.5/
  python setup.py build
  sudo python setup.py install
  cd ..
  # 3. Install R
  # https://cloud.r-project.org/
  # 4. Install the R packages:  
  R
  install.packages("inline") #select 39
  install.packages("gam")
  install.packages("Rcpp")
  install.packages("ggplot2")
  quit("yes")
  # 5.
  sudo apt-get install libgsl0-dev
  # 6. 
  R
  install.packages("RcppGSL") #select 39
  quit("yes")
fi
#==============================================================================
# GET SAMTOOLS 1.3.1
#
if [[ "$GET_SAMTOOLS" -eq "1" ]]; then
  rm -f samtools-1.3.1.tar.bz2
  wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
  tar -xvf samtools-1.3.1.tar.bz2
  cd samtools-1.3.1/
  ./configure --without-curses
  make
  cp samtools ../
  cd ..
  rm -fr samtools-1.3.1.*
fi
#==============================================================================
# GET BOWTIE
#
if [[ "$GET_BOWTIE" -eq "1" ]]; then
  sudo apt-get install libtbb-dev
  rm -fr bowtie/
  git clone https://github.com/BenLangmead/bowtie.git
  cd bowtie/
  make
  cp bowtie ../
  cp bowtie-build ../
  cd ..
fi
#==============================================================================
# BUILD SAMPLE
if [[ "$BUILD_SAMPLE" -eq "1" ]]; then
  rm -f NEAN.fq;
  for((x=25 ; x<=56 ; ++x)); # ONLY UNMAPPED DATA
    do
    ./samtools bam2fq HN-C$x.bam >> NEAN.fq;
    rm -f HN-C$x;
    done
fi
#==============================================================================
###############################################################################
###############################################################################
############################################################################### 
#==============================================================================
# RUN FALCON
if [[ "$RUN_FALCON" -eq "1" ]]; then
  (time ./FALCON -v -n 8 -t 2000 -F -Z -m 20:100:1:5/10 -c 200 -y complexity.nean NEAN DB.fa ) &> REPORT-FALCON ;
  (time ./FALCON-FILTER -v -F -sl 0.001 -du 20000000 -t 0.5 -o positions.nean complexity.nean ) &> REPORT-FALCON-FILTER ;
  (time ./FALCON-EYE -v -F  -e 500 -s 4 -sl 4.15 -o draw.map positions.nean ) &> REPORT-FALCON-EYE ;
fi
#==============================================================================
# RUN GULL
if [[ "$RUN_GULL" -eq "1" ]]; then
  cat top.csv | awk '{ if($3 > 0.001) print $1"\t"$2"\t"$3"\t"$4; }' \
  | awk '{ print $4;}' | tr '|' '\t' | awk '{ print $2;}' > GIS;
  idx=0;
  cat GIS | while read line
    do
    namex=`echo $line | tr ' ' '_'`;
    if [[ "$idx" -eq "0" ]]; then
      printf "%s" "$namex" > FNAMES.fil;
      else
      printf ":%s" "$namex" >> FNAMES.fil;
    fi
    ./goose-extractreadbypattern $line < DB.fa > $namex;
    ((idx++));
    done
  ./GULL-map -v -m 20:100:1:5/10 -c 30 -n 8 -x MATRIX.csv `cat FNAMES.fil`
  ./GULL-visual -v -w 25 -a 8 -x HEATMAP.svg MATRIX.csv
fi
#==============================================================================
###############################################################################
###############################################################################
###############################################################################                                                                                 
#==============================================================================
# RUN KESTREL
if [[ "$RUN_KESTREL" -eq "1" ]]; then
  (time ./KESTREL -v -n 10 -F -m 13:200:1:4/10 -m 6:1:0:0/0 -g 0.9 -o filtered.fq -t 0.7 $1 $2 ) &> REPORT;
  # IF ON -> ALIGN NOT USING $1, BUT RATHER filtered.fq
fi

#==============================================================================
# RUN ALIGNMENT
#
if [[ "$RUN_BOWTIE" -eq "1" ]]; then
  ./bowtie-build $1 index_file
  ./bowtie -a -v 3 --sam index_file $2 | ./samtools view -bh > aligned.bam
fi
#==============================================================================
# RUN MAPDAMAGE
#
if [[ "$RUN_MAPDAMAGE" -eq "1" ]]; then
  mapDamage -i aligned.bam -r $1
fi 
