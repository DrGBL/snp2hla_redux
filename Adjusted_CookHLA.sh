#!/bin/bash

## Adjust_Make_Genetic_Map
## This function calls the MakeGeneticMap function, but allows the use of beagle output from our MakeReference, since it gives outputs that are different from what CookHLA's MakeGeneticMap expects as inputs
#This only works for GRCh38

if [ "$1" == "--help" ]; then
  echo "--path_Tapas: path to HLA-TAPAS folder"
  echo "--path_Ref: path to reference from MakeGeneticMap (plink files, bgl.phased.vcf.gz file, FRQ.frq file, and .markers file."
  echo "--n_rounds: number of MCMC mach iteration"
  echo "--iter: number of MCMC beagle iterations"
  echo "--burnin: number of beagle burnins"
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
    --burnin)
      burnin="$2"
      shift # past argument
      shift
      ;;
    --n_iter)
      n_iter="$2"
      shift # past argument
      shift
      ;;
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

#first need to wrangle the beagle inputs a bit, since our MakeReference gives outputs that are different from what CookHLA's MakeGeneticMap expects as inputs
#specifically the vcf.gz file needs to be brought back to an old beagle format
zcat ${pathRef}.bgl.phased.vcf.gz | java -jar dependency/vcf2beagle.jar "." ${pathRef}
gzip -d ${pathRef}.bgl.gz
mv ${pathRef}.bgl ${pathRef}.bgl.phased_tmp

#now wrangle the header of the bgl.phased_tmp file
head -n 1 ${pathRef}.bgl.phased_tmp | \
  sed 's|^I|P|g' | \
  sed 's|id|pedigree|g' > ${pathOut}header_bgl_file.txt
head -n 1 ${pathRef}.bgl.phased_tmp >> ${pathOut}header_bgl_file.txt

n_col=$(awk '{print NF}' ${pathOut}header_bgl_file.txt | sort -nu | tail -n 1)
n_zero=$(($n_col-2))
printf "change1 change2" > ${pathOut}zeros_file.txt
for (( c=1; c<=${n_zero}; c++ ))
do
  echo "0" | paste -d " " ${pathOut}zeros_file.txt - > ${pathOut}zeros_file_tmp.txt
  cp ${pathOut}zeros_file_tmp.txt ${pathOut}zeros_file.txt
done

sed 's|change1|fID|g' ${pathOut}zeros_file.txt |  \
  sed 's|change2|father|g' >> ${pathOut}header_bgl_file.txt
sed 's|change1|mID|g' ${pathOut}zeros_file.txt |  \
  sed 's|change2|mother|g' >> ${pathOut}header_bgl_file.txt
sed 's|change1|C|g' ${pathOut}zeros_file.txt |  \
  sed 's|change2|gender|g' >> ${pathOut}header_bgl_file.txt

#now put it all together and remove temporary files.
tail -n +2 ${pathRef}.bgl.phased_tmp | cat ${pathOut}header_bgl_file.txt - > ${pathRef}.bgl.phased
rm ${pathRef}.bgl.phased_tmp
rm ${pathOut}header_bgl_file.txt
rm ${pathOut}zeros_file.txt
rm ${pathOut}zeros_file_tmp.txt

#lastly make sure that the allele order is the same in the plink files, so rearrange freq table and plink files
#plink \
#  --vcf ${pathRef}.bgl.phased.vcf.gz \
#  --make-bed \
#  --threads 5 \
#  --memory 15000 \
#  --freq \
#  --double-id \
#  --out ${pathRef}

#mv ${pathRef}.frq ${pathRef}.FRQ.frq

#now CookHLA without genetic map (will be done automatically)
python CookHLA.py \
    -i ${pathTarget} \
    -hg 38 \
    -o ${pathOut}${prefixOut} \
    -ref ${pathRef} \
    --iter ${iter} \
    --burnin ${burnin} \
    -mem 10g \
    --nthreads 10




























