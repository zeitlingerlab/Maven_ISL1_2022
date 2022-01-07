#!/bin/bash

#Melanie Weilert
#January 2020
#For .bam files with large amount of reads, split by chromosome, deduplicate, and remerge.
#Without splitting tasks, R cannot process such a large sample and returns corrupted results.
#Note: This is meant only to be run on single-end aligned BAM files.

#Set defaults
discard=""
gpu=0
memfrac_gpu=.8
nworkers=16

helpmessage="Usage:\n\
bash determine_nexus_barcodes.sh -m [model_dir] -o [output_prefix] -t [tasks] -r [regions] -c [chrom_sizes] -b [BSgenome]\n\
-h  Display this help message.\n\
-m  Model directory.\n\
-o  Output path/prefix for bigwig files.\n\
-t  Tasks to merge, separated by a comma. \n\
-r  Regions to generate predictions/contribution across.\n\
-c  Chromosome sizes in a .txt file.\n\
-b  BSgenome object for bigwig merging in an R script. \n\
-n  Number of workers to use while merging bigwig files [default = 16] \n\
-g  GPU number [default = 0]\n\
-x  Fraction of GPU memory to use [default = .8]\n"

#Option parsing
while getopts "m:o:t:r:c:b:n:g:x:h" opt; do
  case ${opt} in
    m )
      model_dir=$OPTARG
      ;;
    o )
      output_path=$OPTARG
      ;;
    t )
      tasks=$OPTARG
      ;;
    r )
      regions=$OPTARG
      ;;
    c )
      chrom_sizes=$OPTARG
      ;;
    b )
      bsgenome=$OPTARG
      ;;
    n )
      nworkers=$OPTARG
      ;;
    g )
      gpu=$OPTARG
      ;;
    x )
      memfrac_gpu=$OPTARG
      ;;
    h )
      echo -e $helpmessage
      exit 1
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument. Try -h" 1>&2
      exit 1
      ;;
  esac
done

# shift $((OPTIND -1)) #parse through options
if [ -z "$*" ]; then echo "No args, please type: bash run_export_bw_and_merge.sh -h"; exit 0; fi #create error message.

#Run the Export Large BW file
python /n/projects/mw2098/shared_code/bpnet/bpnet_export_large_bw.py --gpu $gpu --memfrac-gpu $memfrac_gpu \
--contrib-method deeplift -r $regions -m $model_dir -o $output_path -c $chrom_sizes

for task in $(echo $tasks | sed "s/,/ /g")
do

  echo "Merging bigwigs for following task: ${task}..."
  Rscript /n/projects/mw2098/shared_code/rscripts/merge_bigwigs.r --files=$output_path\_chr*.$task\.contrib.counts.bw -o $output_path\_$task.contrib.counts.bw -g $bsgenome -c $nworkers
  Rscript /n/projects/mw2098/shared_code/rscripts/merge_bigwigs.r --files=$output_path\_chr*.$task\.contrib.profile.bw -o $output_path\_$task.contrib.profile.bw -g $bsgenome -c $nworkers
  Rscript /n/projects/mw2098/shared_code/rscripts/merge_bigwigs.r --files=$output_path\_chr*.$task\.preds.pos.bw -o $output_path\_$task.preds.pos.bw -g $bsgenome -c $nworkers
  Rscript /n/projects/mw2098/shared_code/rscripts/merge_bigwigs.r --files=$output_path\_chr*.$task\.preds.neg.bw -o $output_path\_$task.preds.neg.bw -g $bsgenome -c $nworkers

done

echo "Done! Cleaning up intermediate files..."

rm $output_path\_chr*

echo "Done for real!"
