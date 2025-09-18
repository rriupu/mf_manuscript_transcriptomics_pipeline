DIR=$(realpath $(dirname $0))

raw_TSS_path=${DIR}/TSS_annotations/raw_annotations.bed
chrom_sizes=${DIR}/hg38.chromsizes

Rscript ${DIR}/get_entrez_ids.R \
    -i $raw_TSS_path \
    -o ${DIR}/TSS_annotations/mRNA_entrez_annotations.bed

chrs=($(cut -f 1 ${DIR}/hg38.chromsizes))
for chr in ${chrs[@]}; do

    grep -e "$chr\t" ${DIR}/TSS_annotations/mRNA_entrez_annotations.bed \
    >> ${DIR}/TSS_annotations/mRNA_entrez_annotations_standard_chrs.bed

done

unique_gene_ids=($(cut -f 4 ${DIR}/TSS_annotations/mRNA_entrez_annotations_standard_chrs.bed | sort -k1,1 | uniq))

for gene in ${unique_gene_ids[@]}; do

    grep -e "\t${gene}\t" ${DIR}/TSS_annotations/mRNA_entrez_annotations_standard_chrs.bed \
    | awk -v OFS='\t' '{print $4"::"$6"::"$1, $2, $3}' \
    | bedtools sort \
    | bedtools merge \
    | awk -v OFS='\t' '{split($1, a, "::"); print a[3], $2, $3, a[1], 0, a[2]}' \
    | bedtools slop -l 300 -r 100 -s -i stdin -g $chrom_sizes \
    >> ${DIR}/TSS_annotations/mRNA_entrez_annotations_processed.bed

done

bedtools sort -i ${DIR}/TSS_annotations/mRNA_entrez_annotations_processed.bed \
> ${DIR}/TSS_annotations/mRNA_entrez_annotations_processed_sorted.bed