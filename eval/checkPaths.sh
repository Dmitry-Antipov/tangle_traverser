#!/bin/bash
module load seqtk
module load winnowmap

if [ "$#" -lt 2 ]; then
    echo "Usage: $0  <gap_fasta> <ont reads fastq>"
    echo "All relevant reads should be in single fastq file"
    exit 1
fi

base_path=$(dirname $(readlink -e $0))


gap_seq=$1
reads=$2

# compress the reads

prefix=`echo $reads |sed s/.fastq//g`
compressed_reads=${prefix}_compress.fasta.gz
if [ ! -e $compressed_reads ];then
   echo "Compressing $reads"
   cat $reads  | seqtk hpc |pigz -c - >${compressed_reads}
fi


name=`echo $gap_seq |sed s/.fasta//g`
alignment_paf=$name.paf
if [ ! -e $alignment_paf ]; then
winnowmap -c -t 16 -x map-ont -H $gap_seq $compressed_reads > $alignment_paf
echo " winnowmapped"
fi

# Process alignment_paf to keep only the read with largest alignment length for each query
echo "Filtering alignments to keep only the best read per query..."
awk '
   BEGIN { FS="\t"; OFS="\t" }
   {
      if (!max[$1] || $10 > max[$1]) {
         max[$1] = $10
         best[$1] = $0
      }
   }
   END {
      for (query in best) {
         print best[query]
      }
   }
' "$alignment_paf" > "${name}_best_alignments.paf"


# Sum up the total matches in largest alignment for each read.
total_alignment_length=$(awk '{sum+=$10} END {print sum}' "${name}_best_alignments.paf")
echo "Total alignment length across all best alignments: $total_alignment_length"
