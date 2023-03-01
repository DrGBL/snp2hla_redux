
args<- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(stringr)

input_beagle <- read.table(args[1],stringsAsFactors = F)
input_marker <- read.table(args[2],stringsAsFactors = F)
ref_beagle <- read.table(args[3],stringsAsFactors = F)
ref_marker <- read.table(args[4],stringsAsFactors = F)


#input_beagle <- read.table("/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/topmed_munged_for_hydra_ukb_hla_imputation.COPY.subset.bgl.phased",stringsAsFactors = F)
#input_marker <- read.table("/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/topmed_munged_for_hydra_ukb_hla_imputation.COPY.subset.markers",stringsAsFactors = F)
#ref_beagle <- read.table("/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/beagle5_typed_all_ukb_topmed_reference_two_fields_fixed.subset.bgl.phased",stringsAsFactors = F)
#ref_marker <- read.table("/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/beagle5_typed_all_ukb_topmed_reference_two_fields_fixed.subset.markers",stringsAsFactors = F)



#only retrieve markers in both
input_beagle<-input_beagle %>%
  filter(V2 %in% ref_beagle$V2)
ref_beagle<-ref_beagle %>%
  filter(V2 %in% input_beagle$V2)

input_marker<-input_marker %>%
  filter(V1 %in% ref_marker$V1) #%>%
  #arrange(V1) %>%
  #distinct() %>%
  #group_by(V1) %>%
  #filter(row_number()==1)
ref_marker<-ref_marker %>%
  filter(V1 %in% input_marker$V1) #%>%
  #arrange(V1) %>%
  #distinct() #%>%
  #group_by(V1) %>%
  #filter(row_number()==1)
  
#test<-read.table("/project/richards/guillaume.butler-laporte/hla_tapas/my_tapas/bangladesh_imputation/test_files/beagle5_typed_all_ukb_topmed_reference_two_fields_fixed.subset.markers",stringsAsFactors = F) %>%
#  filter(V1 %in% input_marker$V1) %>%
#  arrange(V1) %>%
#  distinct() %>%
#  group_by(V1) %>%
#  filter(n()>1)

var_to_remove<-c()
for(i in 1:nrow(ref_marker)){
  #first need to rearrange those god-awful A and T for indels
  #input_ref<-gsub(":", "", str_extract(input_marker[i,1], ":[ACGT]*:"))
  #input_alt<-gsub(":", "", str_extract(ref_marker[i,1], ":[ACGT]*$"))
  #if(nchar(input_ref)>1 | nchar(input_alt)>1){
  #  if(input_marker[i,3]=="A"){
  #    input_marker[i,3]<-"0"
  #    input_marker[i,4]<-"1"
  #  } else {
  #    input_marker[i,4]<-"1"
  #    input_marker[i,3]<-"0"
  #  }
  #  input_beagle[i+5,which(input_beagle[i+5,]=="A")]<-"0"
  #  input_beagle[i+5,which(input_beagle[i+5,]=="T")]<-"1"
  #}
  
  a1<-input_marker[i,3]
  a2<-input_marker[i,4]
  b1<-ref_marker[i,3]
  b2<-ref_marker[i,4]
  
#  if(b1==0 | b1==0){
#    ref<-gsub(":", "", str_extract(ref_marker[i,1], ":[ACGT]*:"))
#    alt<-gsub(":", "", str_extract(ref_marker[i,1], ":[ACGT]*$"))
#    if(ref_marker[i,3]=="0"){
#      ref_marker[i,3]<-ref
#      ref_marker[i,4]<-alt
#    } else {
#      ref_marker[i,4]<-ref
#      ref_marker[i,3]<-alt
#    }
#    ref_beagle[i+5,which(ref_beagle[i+5,]=="0")]<-ref
#    ref_beagle[i+5,which(ref_beagle[i+5,]=="1")]<-alt
#    
#    b1<-ref_marker[i,3]
#    b2<-ref_marker[i,4]
#  }
  
  if( !( (a1==b1 & a2==b2) | (a1==b2 & a2==b1) ) ){
    #var_to_remove<-c(var_to_remove,input_marker[i,1])
    next
    stop("Alleles are not the same for certain variants in the input and reference panels")
  }
  
  if(a1 == b2 & a2==b1) {
    ref_marker[i,3]<-b2
    ref_marker[i,4]<-b1
  }
}

#now do the same as the bgl2GC_trick_bgl R script from Cook HLA
input_beagle<-as.matrix(input_beagle)
input_marker<-as.matrix(input_marker)
ref_beagle<-as.matrix(ref_beagle)
ref_marker<-as.matrix(ref_marker)

first_order_chr <-c("G")
second_order_chr <-c("C")
thrid_order_chr <-c("G")

startSNP_row<-6

# for the input first
beagle_row <- (nrow(input_beagle))
beagle_col <- (ncol(input_beagle))
number_snp <- (nrow(input_marker))


#checking not including snp

making_frist_line=beagle_row-number_snp

print("if this process is stop,please sort the marker and bgl by position(position is adjusted by redefineBPv1BH.py)")

for(t in (startSNP_row:beagle_row))
{
  stopifnot(input_beagle[t,2]==input_marker[t-making_frist_line,1])
  
  
}  

print("good snp order (marker and bgl_data) ")

for(g in (startSNP_row:beagle_row))
{
  snpline <- input_beagle[g,-(1:2)] 
  is.char1 = (snpline == input_marker[g-making_frist_line,3])
  is.char2 = (snpline == input_marker[g-making_frist_line,4])
  is.char3 = !is.char1 & !is.char2
  snpline[is.char1] <- first_order_chr
  snpline[is.char2] <- second_order_chr
  snpline[is.char3] <- thrid_order_chr
  input_beagle[g,-(1:2)] <- snpline 
  
  
  
}

#conversing GC overlap markers

for(k in (1:nrow(input_marker)))
{
  
  input_marker[,3]<-first_order_chr
  
  input_marker[,4]<-second_order_chr 
  
}

save_bgl<-args[5]
save_marker<-args[6]

write.table(input_beagle,save_bgl,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(input_marker,save_marker,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)




# for the ref now
beagle_row <- (nrow(ref_beagle))
beagle_col <- (ncol(ref_beagle))
number_snp <- (nrow(ref_marker))


#checking not including snp

making_frist_line=beagle_row-number_snp

print("if this process is stop,please sort the marker and bgl by position(position is adjusted by redefineBPv1BH.py)")

for(t in (startSNP_row:beagle_row))
{
  stopifnot(ref_beagle[t,2]==ref_marker[t-making_frist_line,1])
  
  
}  

print("good snp order (marker and bgl_data) ")

for(g in (startSNP_row:beagle_row))
{
  snpline <- ref_beagle[g,-(1:2)] 
  is.char1 = (snpline == ref_marker[g-making_frist_line,3])
  is.char2 = (snpline == ref_marker[g-making_frist_line,4])
  is.char3 = !is.char1 & !is.char2
  snpline[is.char1] <- first_order_chr
  snpline[is.char2] <- second_order_chr
  snpline[is.char3] <- thrid_order_chr
  ref_beagle[g,-(1:2)] <- snpline 
  
  
  
}

#conversing GC overlap markers

for(k in (1:nrow(ref_marker)))
{
  
  ref_marker[,3]<-first_order_chr
  
  ref_marker[,4]<-second_order_chr 
  
}

save_bgl<-args[7]
save_marker<-args[8]

write.table(ref_beagle,save_bgl,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(ref_marker,save_marker,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)


































