#!/bin/bash

## Adjust_Make_Genetic_Map
## This function calls the MakeGeneticMap function, but allows the use of beagle output from our MakeReference, since it gives outputs that are different from what CookHLA's MakeGeneticMap expects as inputs
#This only works for GRCh38

if [ "$1" == "--help" ]; then
  echo "--path_Tapas: path to HLA-TAPAS folder"
  echo "--path_Ref: path to reference from MakeGeneticMap (plink files, bgl.phased.vcf.gz file, FRQ.frq file, and .markers file."
  echo "--n_rounds: number of MCMC mach1 iteration"
  echo "--path_out: path to output folder"
  echo "--prefix_out: prefix of files written in the output folder"
  echo "--path_target: path to plink files of target samples (i.e. those to be eventually imputed for HLA alleles)"
  echo "--help: this help text"
  exit 0
fi


POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --path_Tapas)
      pathTapas="$2"
      shift # past argument
      shift # past value
      ;;
    --path_Ref)
      pathRef="$2"
      shift # past argument
      shift # past value
      ;;
    --n_rounds)
      rounds="$2"
      shift # past argument
      shift
      ;;
    --path_out)
      pathOut="$2"
      shift # past argument
      shift
      ;;
    --prefix_out)
      prefixOut="$2"
      shift # past argument
      shift
      ;;
    --path_target)
      pathTarget="$2"
      shift # past argument
      shift
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

cd ${pathTapas}

#need to wrangle the beagle inputs a bit, since our MakeReference gives outputs that are different from what CookHLA's MakeGeneticMap expects as inputs

# first need to rename the variants since makereference does weird things with indels, and renames the chromosome wrong, and only keep the variants in the target file
#bcftools query -f '%ID\n' ${pathTarget}.vcf.gz > ${pathTarget}_list_variants.txt

bcftools view --header-only ${pathRef}.bgl.phased.vcf.gz  > ${pathRef}.bgl.phased_header.txt
tabix -p vcf ${pathRef}.bgl.phased.vcf.gz
bcftools query -f '%ID\n' ${pathRef}.bgl.phased.vcf.gz > ${pathRef}.bgl.phased_list.variants.txt
sed 's|^chr6:[0-9]*:||g' ${pathRef}.bgl.phased_list.variants.txt | \
  sed 's|:[a-zA-Z]*$||g' | \
  sed 's|^HLA.*|A|g' > ${pathRef}.bgl.phased_ref_alleles.txt
sed 's|^chr6:[0-9]*:||g' ${pathRef}.bgl.phased_list.variants.txt | \
  sed 's|^[a-zA-Z]*:||g' | \
  sed 's|^HLA.*|T|g' > ${pathRef}.bgl.phased_alt_alleles.txt

bcftools view --no-header ${pathRef}.bgl.phased.vcf.gz | \
  awk 'FNR==NR{a[NR]=$1;next}{$4=a[FNR]}1' ${pathRef}.bgl.phased_ref_alleles.txt - | \
  awk 'FNR==NR{a[NR]=$1;next}{$5=a[FNR]}1' ${pathRef}.bgl.phased_alt_alleles.txt - | \
  tr [:blank:] \\t | \
  cat ${pathRef}.bgl.phased_header.txt - | \
  sed 's|^6|chr6|g' | \
  sed 's|##contig=<ID=6>|##contig=<ID=chr6>|g' | \
  bgzip > ${pathRef}_fixed.bgl.phased.vcf.gz
  #bcftools view -i "ID=@${pathTarget}_list_variants.txt" -Oz > ${pathRef}_fixed.bgl.phased.vcf.gz

#make plink files
plink \
  --vcf ${pathRef}_fixed.bgl.phased.vcf.gz \
  --double-id \
  --make-bed \
  --threads 5 \
  --memory 30000 \
  --freq \
  --out ${pathRef}_fixed

#make FRQ.frq file
mv ${pathRef}_fixed.frq ${pathRef}_fixed.FRQ.frq

#make markers file
awk '{print $2, $4, $6, $5}' ${pathRef}_fixed.bim | tr [:blank:] \\t  > ${pathRef}_fixed.markers


#now the vcf.gz file needs to be brought back to an old beagle format
zcat ${pathRef}_fixed.bgl.phased.vcf.gz | java -jar dependency/vcf2beagle.jar "." ${pathRef}_fixed
gzip -d -f ${pathRef}_fixed.bgl.gz
mv ${pathRef}_fixed.bgl ${pathRef}.bgl.phased_tmp

#now wrangle the header of the bgl.phased_tmp file
head -n 1 ${pathRef}.bgl.phased_tmp | \
  sed 's|^I|P|g' | \
  sed 's|id|pedigree|g' > ${pathOut}header_bgl_file.txt
head -n 1 ${pathRef}.bgl.phased_tmp >> ${pathOut}header_bgl_file.txt

n_col=$(awk '{print NF}' ${pathOut}header_bgl_file.txt | sort -nu | tail -n 1)
n_zero=$(($n_col-2))
printf "change1 change2" > ${pathOut}zeros_file_tmp.txt

awk "BEGIN{for(c=0;c<${n_zero};c++) print "0"}" > ${pathRef}tmp.txt
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
}' ${pathRef}tmp.txt | paste ${pathOut}zeros_file_tmp.txt - > ${pathOut}zeros_file.txt

rm ${pathOut}zeros_file_tmp.txt

sed 's|change1|fID|g' ${pathOut}zeros_file.txt |  \
  sed 's|change2|father|g' >> ${pathOut}header_bgl_file.txt
sed 's|change1|mID|g' ${pathOut}zeros_file.txt |  \
  sed 's|change2|mother|g' >> ${pathOut}header_bgl_file.txt
sed 's|change1|C|g' ${pathOut}zeros_file.txt |  \
  sed 's|change2|gender|g' >> ${pathOut}header_bgl_file.txt

#now put it all together and remove temporary files.
tail -n +2 ${pathRef}.bgl.phased_tmp | cat ${pathOut}header_bgl_file.txt - > ${pathRef}_fixed.bgl.phased

#make markers file
awk '{print $2, $4, $6, $5}' ${pathRef}_fixed.bim | tr [:blank:] \\t  > ${pathRef}_fixed.markers

#remove intermediate files
rm ${pathRef}.bgl.phased_tmp
rm ${pathOut}header_bgl_file.txt
rm ${pathOut}zeros_file.txt
rm ${pathOut}zeros_file_tmp.txt
rm ${pathRef}.bgl.phased_ref_alleles.txt
rm ${pathRef}.bgl.phased_alt_alleles.txt
rm ${pathRef}.bgl.phased_header.txt
rm ${pathRef}.bgl.phased_list.variants.txt


#now call genetic map
python -m MakeGeneticMap \
    -i ${pathTarget} \
    -hg 38 \
    -ref ${pathRef}_fixed \
    --mach_mcmc_rounds ${rounds} \
    -o ${pathOut}${prefixOut}_fixed































