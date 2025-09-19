DIR=$(realpath $(dirname $0))

raw_TSS_path=$1
chrom_sizes=$2
output_dir=$3

Rscript ${DIR}/get_entrez_ids.R \
    -i $raw_TSS_path \
    -o ${output_dir}/mRNA_entrez_annotations.bed

chrs=($(cut -f 1 $chrom_sizes))
for chr in ${chrs[@]}; do

    grep -P "$chr\t" ${output_dir}/mRNA_entrez_annotations.bed \
    >> ${output_dir}/mRNA_entrez_annotations_standard_chrs.bed

done

unique_gene_ids=($(cut -f 4 ${output_dir}/mRNA_entrez_annotations_standard_chrs.bed | sort -k1,1 | uniq))

for gene in ${unique_gene_ids[@]}; do

    # grep -P "\t${gene}\t" ${output_dir}/mRNA_entrez_annotations_standard_chrs.bed \
    awk -v OFS='\t' -v gene=${gene} '$4 == gene {print $4"::"$6"::"$1, $2, $3}' ${output_dir}/mRNA_entrez_annotations_standard_chrs.bed \
    | bedtools sort \
    | bedtools merge \
    | awk -v OFS='\t' '{split($1, a, "::"); print a[3], $2, $3, a[1], 0, a[2]}' \
    | bedtools slop -l 300 -r 100 -s -i stdin -g $chrom_sizes \
    >> ${output_dir}/mRNA_entrez_annotations_processed.bed

done

bedtools sort -i ${output_dir}/mRNA_entrez_annotations_processed.bed \
> ${output_dir}/mRNA_entrez_annotations_processed_sorted.bed