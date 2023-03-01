#!/bin/bash

## Adjust_Make_Genetic_Map
## This function adjusts outputs from the original MakeReference so that they're compatible with CookHLA's MakeGeneticMap expects as inputs
# This only works for GRCh38

if [ "$1" == "--help" ]; then
  echo "--path_Ref: path to reference"
  echo "--help: this help text"
  exit 0
fi


POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --path_Ref)
      pathRef="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

#need to wrangle the beagle inputs a bit, since our MakeReference gives outputs that are different from what CookHLA's MakeGeneticMap expects as inputs

# first need to rename the variants since makereference does weird things with indels, and renames the chromosome wrong, and only keep the variants in the target file
#bcftools query -f '%ID\n' ${pathTarget}.vcf.gz > ${pathTarget}_list_variants.txt

#bcftools view --header-only ${pathRef}_tmp.bgl.phased.vcf.gz  > ${pathRef}.bgl.phased_header.txt
#tabix -p vcf ${pathRef}_tmp.bgl.phased.vcf.gz
#bcftools query -f '%ID\n' ${pathRef}_tmp.bgl.phased.vcf.gz > ${pathRef}.bgl.phased_list.variants.txt
#sed 's|^chr6:[0-9]*:||g' ${pathRef}.bgl.phased_list.variants.txt | \
#  sed 's|:[a-zA-Z]*$||g' | \
#  sed 's|^HLA.*|A|g' > ${pathRef}.bgl.phased_ref_alleles.txt
#sed 's|^chr6:[0-9]*:||g' ${pathRef}.bgl.phased_list.variants.txt | \
#  sed 's|^[a-zA-Z]*:||g' | \
#  sed 's|^HLA.*|T|g' > ${pathRef}.bgl.phased_alt_alleles.txt

#bcftools view --no-header ${pathRef}_tmp.bgl.phased.vcf.gz | \
#  awk 'FNR==NR{a[NR]=$1;next}{$4=a[FNR]}1' ${pathRef}.bgl.phased_ref_alleles.txt - | \
#  awk 'FNR==NR{a[NR]=$1;next}{$5=a[FNR]}1' ${pathRef}.bgl.phased_alt_alleles.txt - | \
#  tr [:blank:] \\t | \
#  cat ${pathRef}.bgl.phased_header.txt - | \
#  sed 's|^6|chr6|g' | \
#  sed 's|##contig=<ID=6>|##contig=<ID=chr6>|g' | \
#  bgzip > ${pathRef}.bgl.phased.vcf.gz

#make plink files
#plink \
#  --vcf ${pathRef}.bgl.phased.vcf.gz \
#  --double-id \
#  --make-bed \
#  --threads 5 \
#  --memory 30000 \
#  --freq \
#  --out ${pathRef}

#make FRQ.frq file
#mv ${pathRef}.frq ${pathRef}.FRQ.frq

#now the vcf.gz file needs to be brought back to an old beagle format
zcat ${pathRef}.bgl.phased.vcf.gz | java -jar dependency/vcf2beagle.jar "." ${pathRef}
sleep 60s
gzip -d -f ${pathRef}.bgl.gz
mv ${pathRef}.bgl ${pathRef}.bgl.phased_tmp

#now wrangle the header of the bgl.phased_tmp file
head -n 1 ${pathRef}.bgl.phased_tmp | \
  sed 's|^I|P|g' | \
  sed 's|id|pedigree|g' > ${pathRef}header_bgl_file.txt
head -n 1 ${pathRef}.bgl.phased_tmp >> ${pathRef}header_bgl_file.txt

n_col=$(awk '{print NF}' ${pathRef}header_bgl_file.txt | sort -nu | tail -n 1)
n_zero=$(($n_col-2))
printf "change1 change2" > ${pathRef}zeros_file_tmp.txt

awk "BEGIN{for(c=0;c<${n_zero};c++) print "0"}" > ${pathRef}_tmp.txt
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' ${pathRef}_tmp.txt | paste ${pathRef}zeros_file_tmp.txt - > ${pathRef}zeros_file.txt

sed 's|change1|fID|g' ${pathRef}zeros_file.txt |  \
  sed 's|change2|father|g' >> ${pathRef}header_bgl_file.txt
sed 's|change1|mID|g' ${pathRef}zeros_file.txt |  \
  sed 's|change2|mother|g' >> ${pathRef}header_bgl_file.txt
sed 's|change1|C|g' ${pathRef}zeros_file.txt |  \
  sed 's|change2|gender|g' >> ${pathRef}header_bgl_file.txt

#now put it all together
tail -n +2 ${pathRef}.bgl.phased_tmp | cat ${pathRef}header_bgl_file.txt - > ${pathRef}.bgl.phased

#make markers file
#awk '{print $2, $4, $6, $5}' ${pathRef}.bim | tr [:blank:] \\t > ${pathRef}.markers

#and remove temporary files.
rm ${pathRef}.bgl.phased_tmp
rm ${pathRef}header_bgl_file.txt
rm ${pathRef}zeros_file.txt
rm ${pathRef}zeros_file_tmp.txt
rm ${pathRef}_tmp.txt
#rm ${pathRef}.bgl.phased_ref_alleles.txt
#rm ${pathRef}.bgl.phased_alt_alleles.txt
#rm ${pathRef}.bgl.phased_header.txt
#rm ${pathRef}.bgl.phased_list.variants.txt