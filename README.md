# SNP2HLA Redux
This is a modification, simplification, and overhaul of **HLA-TAPAS** which handles HLA reference panel construction (*MakeReference*) and HLA imputation (*SNP2HLA*). It also handles *Cook-HLA* and building genetic maps (*MakeGeneticMap*), but these features are in development.

This tool does not handle association studies, for which we recommend the use of proper GWAS software to better adjust for relatedness (e.g. REGENIE, SAIGE).

Modifications:

(1) The native use of GRCh38 instead of successive liftovers to hg19 and hg18. Note that at this point it *requires* that your data be in GRCh38. There is no plan to allow for GRCh37/hg19 or other builds in the future.

(2) More genes.

(3) Using BEAGLE v5.4 for phasing and imputation.

(4) More options to control the QC parameters of *MakeReference* and *SNP2HLA*

There was no change in the algorithm otherwise. Specifically, the actual imputation is only performed by Beagle, and the rest of the software is used as a wrapper for Beagle.

## Credit and citation
You will find that much of the code and documentation was either inspired by or directly pulled from **HLA-TAPAS**.

Hence, please cite [the following paper](https://www.nature.com/articles/s41588-021-00935-7) if you use this work.

Luo, Y., Kanai, M., Choi, W. et al. A high-resolution HLA reference panel capturing global population diversity enables multi-ancestry fine-mapping in HIV host response. Nat Genet 53, 1504–1516 (2021). https://doi.org/10.1038/s41588-021-00935-7

## Requirements & Dependencies

To install the software, simply clone this git in your environment.

This software was tested in a CentOS environment.

Python system requires next settings.
- python>=3.7
- pandas>=1.0.3

R statistical programming language requires next settings.
- R>=3.6
- argparse
- stringr
- purrr
- dplyr
- multidplyr
- tidyr
- data.table
- parallel
- rcompanion

You will also need the softwares in the 'dependency/' folder i.e. Plink1.9, software (mach) from the Gonçalo Abecasis lab, and code/software from the Brian Browning lab (including BEAGLE 5.4). All copyrights to their respective authors.

Once downloaded, you might need to change the file permission for PLINK.
```
$ chmod +x dependency/plink
```

## Usage

### Building reference panel

All code blocks below should be ran from the snp2hla_redux root folder.

```
python -m MakeReference \
  --variants reference_variant_panel_plink_prefix \
  --chped hla_calls.ped \
  --genes A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1 \ #order is important, you can add more genes, see help command output
  --hg 38 \   #only 38 is allowed
  --mind 0.3 \
  --hardy 0.00000005 \
  --maf 0.00000005 \
  --miss 0.05 \
  --hla_maf 0.00000005 \
  --out full_path_out_reference \
  --mem 50G \
  --burnin 20 \
  --iter 100 \
  --nthreads 20 \
  --phasing \
  --window 10 \
  --overlap 1.8 \
  --tmp /tmp
```

Here *reference_variant_panel_plink_prefix* is the prefix of the plink files of the reference sample and hla_calls.ped is a ped file with the HLA alleles. For example, for an individual with IID and FID 1000000, for which we have three HLA genes (A, B, and C), their row would be:

```
1000000 1000000 0       0       0       0       A*01:01 A*02:02 B*01:01 B*02:02 C*01:01 C*02:02
```

Full details of the other options can be read using `python -m MakeReference --help`.


### Imputation

```
python -m SNP2HLA \
  --target target_variants_plink_prefix \
  --reference path_out_reference \
  --out full_path_out_imputed \
  --nthreads 20 \
  --burnin 10 \
  --iter 100 \
  --window 40 \
  --overlap 2 \
  --maf 0.00005 \
  --mem 30g
```

Here, the *path_out_reference* is the same as from *MakeReference* above. Again, full details can be obtained with `python -m SNP2HLA --help`.

### A few notes on the output files

The output is a vcf with the reference panel variants and the HLA alleles, all imputed. Each imputed variant and alleles will include an R2 value. This should not be interpreted as a dosage. Please refer to the Beagle documentation for this.

Also note that the major assumption made by SNP2HLA (and hence by SNP2HLA_reudx) is that each HLA allele is made into a biallelic SNP, and imputed as such. Hence, *it does not use the fact that HLA genes are multiallelic, and will sometimes impute more than 2 possible alleles for a gene, for a given sample*. Users of this software should be aware of this for their QC.
