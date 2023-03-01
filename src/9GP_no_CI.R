library(stringr)

arg <- commandArgs(trailingOnly = TRUE)
args <- arg[1]
args2 <- arg[2]
args3 <- arg[3]
args4 <- arg[4]
args5 <- arg[5]
args6 <- arg[6]
args7 <- arg[7]
args8 <- arg[8]
args9 <- arg[9]

#args <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon2.0.5.raw_imputation_out.vcf.HLA_A"
#args2 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon2.1.5.raw_imputation_out.vcf.HLA_A"
#args3 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon2.1.raw_imputation_out.vcf.HLA_A"
#args4 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon3.0.5.raw_imputation_out.vcf.HLA_A"
#args5 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon3.1.5.raw_imputation_out.vcf.HLA_A"
#args6 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon3.1.raw_imputation_out.vcf.HLA_A"
#args7 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon4.0.5.raw_imputation_out.vcf.HLA_A"
#args8 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon4.1.5.raw_imputation_out.vcf.HLA_A"
#args9 <- "/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/bangladesh_ukb_genetic_map.MHC.QC.exon4.1.raw_imputation_out.vcf.HLA_A"


gene <- arg[10]
HLA_EXON1 <- as.matrix(read.table(args,colClasses = 'character'))
HLA_EXON2 <- as.matrix(read.table(args2,colClasses = 'character'))
HLA_EXON3 <- as.matrix(read.table(args3,colClasses = 'character'))
HLA_EXON4 <- as.matrix(read.table(args4,colClasses = 'character'))
HLA_EXON5 <- as.matrix(read.table(args5,colClasses = 'character'))
HLA_EXON6 <- as.matrix(read.table(args6,colClasses = 'character'))
HLA_EXON7 <- as.matrix(read.table(args7,colClasses = 'character'))
HLA_EXON8 <- as.matrix(read.table(args8,colClasses = 'character'))
HLA_EXON9 <- as.matrix(read.table(args9,colClasses = 'character'))

a <- max(nrow(HLA_EXON1),nrow(HLA_EXON2),nrow(HLA_EXON3),nrow(HLA_EXON4),nrow(HLA_EXON5),nrow(HLA_EXON6),nrow(HLA_EXON7),nrow(HLA_EXON8),nrow(HLA_EXON9))

# Exception handling for absent HLA genes.
if (a <= 1) {
  print(paste0("No HLA_", gene, " in this study."))
  quit(save="no")
}

new_HLA_vcf<-matrix(0,a,ncol(HLA_EXON1))
new_HLA_vcf[,1:9]<-HLA_EXON1[1:a,1:9]

FID<-HLA_EXON1[1,10:ncol(HLA_EXON1)]

for(i in (2:a))
{
  for(j in (10:ncol(HLA_EXON1)))
  {
    #in the previous code GP position 1 was for homozygous for the given HLA allele. For some reason this flipped in GRCh38 with UKB and Beagle 5.4. So it's position 3 now
    #i.e. I use "strsplit((strsplit(HLA_EXON1[i,j],":")[[1]][3]),",")[[1]][3]" instead of "strsplit((strsplit(HLA_EXON1[i,j],":")[[1]][3]),",")[[1]][1]" below
    EXON1 <- as.numeric(strsplit((strsplit(HLA_EXON1[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON1[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON2 <- as.numeric(strsplit((strsplit(HLA_EXON2[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON2[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON3 <- as.numeric(strsplit((strsplit(HLA_EXON3[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON3[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON4 <- as.numeric(strsplit((strsplit(HLA_EXON4[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON4[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON5 <- as.numeric(strsplit((strsplit(HLA_EXON5[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON5[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON6 <- as.numeric(strsplit((strsplit(HLA_EXON6[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON6[i,j],":")[[1]][3]),",")[[1]][2])/2
    if(i <=nrow(HLA_EXON7)){
      EXON7 <- as.numeric(strsplit((strsplit(HLA_EXON7[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON7[i,j],":")[[1]][3]),",")[[1]][2])/2
      EXON8 <- as.numeric(strsplit((strsplit(HLA_EXON8[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON8[i,j],":")[[1]][3]),",")[[1]][2])/2
      EXON9 <- as.numeric(strsplit((strsplit(HLA_EXON9[i,j],":")[[1]][3]),",")[[1]][3])+as.numeric(strsplit((strsplit(HLA_EXON9[i,j],":")[[1]][3]),",")[[1]][2])/2
      new_HLA_vcf[i,j] <- (max(EXON1,EXON2,EXON3)+max(EXON4,EXON5,EXON6)+max(EXON7,EXON8,EXON9))/3
    }
    else{
      new_HLA_vcf[i,j] <- (max(EXON1,EXON2,EXON3)+max(EXON4,EXON5,EXON6))/2
    }
  }
}

# Normalization
for(i in (10:ncol(HLA_EXON1))){
  sum_col <- sum(as.numeric(new_HLA_vcf[,i]))
  for(j in (2:a)){
    new_HLA_vcf[j,i] <- as.numeric(new_HLA_vcf[j,i])/sum_col
  }
}

# [,4] <- First allele, [,5] <- second allele 
# [,6] <- First pp,     [,7] <- second pp     
# [,8] <- Confidence                          


HLA<-matrix(0,ncol(HLA_EXON1)-9,8)
for(i in (10:ncol(new_HLA_vcf))){
  HLA[i-9,1] <- FID[i-9]
  HLA[i-9,2] <- FID[i-9]
  HLA[i-9,3] <- gene
  HLA[i-9,4] <- new_HLA_vcf[which(new_HLA_vcf[,i]==max(new_HLA_vcf[,i])),3][1]
  HLA[i-9,6] <- new_HLA_vcf[which(new_HLA_vcf[,i]==max(new_HLA_vcf[,i])),i][1]
  t <- 0
  max_pp <- 0
  for(j in (2:(nrow(new_HLA_vcf))))
  {
    if (new_HLA_vcf[j,3] != HLA[i-9,4] && new_HLA_vcf[j,i] > as.numeric(HLA[i-9,6])/2 && max_pp < new_HLA_vcf[j,i]){
      HLA[i-9,5] <- new_HLA_vcf[j,3]
      HLA[i-9,7] <- new_HLA_vcf[j,i]
      HLA[i-9,8] <- as.numeric(HLA[i-9,6])+as.numeric(HLA[i-9,7])
      t <- 1
      max_pp <- new_HLA_vcf[j,i]
    }
  }
  if(t==0){
    HLA[i-9,5] <- HLA[i-9,4]
    HLA[i-9,7] <- HLA[i-9,6]
    HLA[i-9,8] <- HLA[i-9,6]
  }
}

for(i in (1:nrow(HLA))){
  #commented out the original code which assumed that alleles were coded as HLA_A_0101_exon2
  #this code assumes that alleles are coded as HLA_A*01:01_exon2
  #i1 <- strsplit(HLA[i,4],"_")[[1]][3]
  #i2 <- strsplit(HLA[i,5],"_")[[1]][3]
  #i3 <- paste0(strsplit(i1,"")[[1]][1],strsplit(i1,"")[[1]][2])
  #i4 <- paste0(strsplit(i2,"")[[1]][1],strsplit(i2,"")[[1]][2])
  #HLA[i,4] <- paste0(i3,",",i4)
  #HLA[i,5] <- paste0(i1,",",i2)
  i1<-gsub(":", "", gsub("\\*", "", paste0(str_extract(HLA[i,4], "\\*[0-9]+:"), ",", str_extract(HLA[i,5], "\\*[0-9]+:"))))
  i2<-gsub("\\*", "", paste0(str_extract(HLA[i,4], "\\*[0-9]+:[0-9]+"),",",str_extract(HLA[i,5], "\\*[0-9]+:[0-9]+")))
  HLA[i,4] <- i1
  HLA[i,5] <- i2
}

write.table(HLA,paste0(args,".alleles"),quote=F, col.names=F, row.names=F)
