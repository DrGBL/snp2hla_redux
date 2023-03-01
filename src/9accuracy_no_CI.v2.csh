set INPUT=$1
set INPUT2=$2
set INPUT3=$3
set INPUT4=$4
set INPUT5=$5
set INPUT6=$6
set INPUT7=$7
set INPUT8=$8
set INPUT9=$9

set out=$10

#Reminder: ["A", "B", "C", "E", "F", "G", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5"]

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT2.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT3.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT4.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT5.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT6.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT7.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT8.HLA_DRB5

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_E
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_F
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_G
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DMA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DMB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DOA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DOB
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DPB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DRA
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DRB3
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DRB4
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT9.HLA_DRB5

grep HLA_A $INPUT >> $INPUT.HLA_A
grep HLA_B $INPUT >> $INPUT.HLA_B
grep HLA_C $INPUT >> $INPUT.HLA_C
grep HLA_E $INPUT >> $INPUT.HLA_E
grep HLA_F $INPUT >> $INPUT.HLA_F
grep HLA_G $INPUT >> $INPUT.HLA_G
grep HLA_DMA $INPUT >> $INPUT.HLA_DMA
grep HLA_DMB $INPUT >> $INPUT.HLA_DMB
grep HLA_DOA $INPUT >> $INPUT.HLA_DOA
grep HLA_DOB $INPUT >> $INPUT.HLA_DOB
grep HLA_DPA1 $INPUT >> $INPUT.HLA_DPA1
grep HLA_DPB1 $INPUT >> $INPUT.HLA_DPB1
grep HLA_DQA1 $INPUT >> $INPUT.HLA_DQA1
grep HLA_DQB1 $INPUT >> $INPUT.HLA_DQB1
grep HLA_DRA $INPUT >> $INPUT.HLA_DRA
grep HLA_DRB1 $INPUT >> $INPUT.HLA_DRB1
grep HLA_DRB3 $INPUT >> $INPUT.HLA_DRB3
grep HLA_DRB4 $INPUT >> $INPUT.HLA_DRB4
grep HLA_DRB5 $INPUT >> $INPUT.HLA_DRB5

grep HLA_A $INPUT2 >> $INPUT2.HLA_A
grep HLA_B $INPUT2 >> $INPUT2.HLA_B
grep HLA_C $INPUT2 >> $INPUT2.HLA_C
grep HLA_E $INPUT2 >> $INPUT2.HLA_E
grep HLA_F $INPUT2 >> $INPUT2.HLA_F
grep HLA_G $INPUT2 >> $INPUT2.HLA_G
grep HLA_DMA $INPUT2 >> $INPUT2.HLA_DMA
grep HLA_DMB $INPUT2 >> $INPUT2.HLA_DMB
grep HLA_DOA $INPUT2 >> $INPUT2.HLA_DOA
grep HLA_DOB $INPUT2 >> $INPUT2.HLA_DOB
grep HLA_DPA1 $INPUT2 >> $INPUT2.HLA_DPA1
grep HLA_DPB1 $INPUT2 >> $INPUT2.HLA_DPB1
grep HLA_DQA1 $INPUT2 >> $INPUT2.HLA_DQA1
grep HLA_DQB1 $INPUT2 >> $INPUT2.HLA_DQB1
grep HLA_DRA $INPUT2 >> $INPUT2.HLA_DRA
grep HLA_DRB1 $INPUT2 >> $INPUT2.HLA_DRB1
grep HLA_DRB3 $INPUT2 >> $INPUT2.HLA_DRB3
grep HLA_DRB4 $INPUT2 >> $INPUT2.HLA_DRB4
grep HLA_DRB5 $INPUT2 >> $INPUT2.HLA_DRB5

grep HLA_A $INPUT3 >> $INPUT3.HLA_A
grep HLA_B $INPUT3 >> $INPUT3.HLA_B
grep HLA_C $INPUT3 >> $INPUT3.HLA_C
grep HLA_E $INPUT3 >> $INPUT3.HLA_E
grep HLA_F $INPUT3 >> $INPUT3.HLA_F
grep HLA_G $INPUT3 >> $INPUT3.HLA_G
grep HLA_DMA $INPUT3 >> $INPUT3.HLA_DMA
grep HLA_DMB $INPUT3 >> $INPUT3.HLA_DMB
grep HLA_DOA $INPUT3 >> $INPUT3.HLA_DOA
grep HLA_DOB $INPUT3 >> $INPUT3.HLA_DOB
grep HLA_DPA1 $INPUT3 >> $INPUT3.HLA_DPA1
grep HLA_DPB1 $INPUT3 >> $INPUT3.HLA_DPB1
grep HLA_DQA1 $INPUT3 >> $INPUT3.HLA_DQA1
grep HLA_DQB1 $INPUT3 >> $INPUT3.HLA_DQB1
grep HLA_DRA $INPUT3 >> $INPUT3.HLA_DRA
grep HLA_DRB1 $INPUT3 >> $INPUT3.HLA_DRB1
grep HLA_DRB3 $INPUT3 >> $INPUT3.HLA_DRB3
grep HLA_DRB4 $INPUT3 >> $INPUT3.HLA_DRB4
grep HLA_DRB5 $INPUT3 >> $INPUT3.HLA_DRB5

grep HLA_A $INPUT4 >> $INPUT4.HLA_A
grep HLA_B $INPUT4 >> $INPUT4.HLA_B
grep HLA_C $INPUT4 >> $INPUT4.HLA_C
grep HLA_E $INPUT4 >> $INPUT4.HLA_E
grep HLA_F $INPUT4 >> $INPUT4.HLA_F
grep HLA_G $INPUT4 >> $INPUT4.HLA_G
grep HLA_DMA $INPUT4 >> $INPUT4.HLA_DMA
grep HLA_DMB $INPUT4 >> $INPUT4.HLA_DMB
grep HLA_DOA $INPUT4 >> $INPUT4.HLA_DOA
grep HLA_DOB $INPUT4 >> $INPUT4.HLA_DOB
grep HLA_DPA1 $INPUT4 >> $INPUT4.HLA_DPA1
grep HLA_DPB1 $INPUT4 >> $INPUT4.HLA_DPB1
grep HLA_DQA1 $INPUT4 >> $INPUT4.HLA_DQA1
grep HLA_DQB1 $INPUT4 >> $INPUT4.HLA_DQB1
grep HLA_DRA $INPUT4 >> $INPUT4.HLA_DRA
grep HLA_DRB1 $INPUT4 >> $INPUT4.HLA_DRB1
grep HLA_DRB3 $INPUT4 >> $INPUT4.HLA_DRB3
grep HLA_DRB4 $INPUT4 >> $INPUT4.HLA_DRB4
grep HLA_DRB5 $INPUT4 >> $INPUT4.HLA_DRB5

grep HLA_A $INPUT5 >> $INPUT5.HLA_A
grep HLA_B $INPUT5 >> $INPUT5.HLA_B
grep HLA_C $INPUT5 >> $INPUT5.HLA_C
grep HLA_E $INPUT5 >> $INPUT5.HLA_E
grep HLA_F $INPUT5 >> $INPUT5.HLA_F
grep HLA_G $INPUT5 >> $INPUT5.HLA_G
grep HLA_DMA $INPUT5 >> $INPUT5.HLA_DMA
grep HLA_DMB $INPUT5 >> $INPUT5.HLA_DMB
grep HLA_DOA $INPUT5 >> $INPUT5.HLA_DOA
grep HLA_DOB $INPUT5 >> $INPUT5.HLA_DOB
grep HLA_DPA1 $INPUT5 >> $INPUT5.HLA_DPA1
grep HLA_DPB1 $INPUT5 >> $INPUT5.HLA_DPB1
grep HLA_DQA1 $INPUT5 >> $INPUT5.HLA_DQA1
grep HLA_DQB1 $INPUT5 >> $INPUT5.HLA_DQB1
grep HLA_DRA $INPUT5 >> $INPUT5.HLA_DRA
grep HLA_DRB1 $INPUT5 >> $INPUT5.HLA_DRB1
grep HLA_DRB3 $INPUT5 >> $INPUT5.HLA_DRB3
grep HLA_DRB4 $INPUT5 >> $INPUT5.HLA_DRB4
grep HLA_DRB5 $INPUT5 >> $INPUT5.HLA_DRB5

grep HLA_A $INPUT6 >> $INPUT6.HLA_A
grep HLA_B $INPUT6 >> $INPUT6.HLA_B
grep HLA_C $INPUT6 >> $INPUT6.HLA_C
grep HLA_E $INPUT6 >> $INPUT6.HLA_E
grep HLA_F $INPUT6 >> $INPUT6.HLA_F
grep HLA_G $INPUT6 >> $INPUT6.HLA_G
grep HLA_DMA $INPUT6 >> $INPUT6.HLA_DMA
grep HLA_DMB $INPUT6 >> $INPUT6.HLA_DMB
grep HLA_DOA $INPUT6 >> $INPUT6.HLA_DOA
grep HLA_DOB $INPUT6 >> $INPUT6.HLA_DOB
grep HLA_DPA1 $INPUT6 >> $INPUT6.HLA_DPA1
grep HLA_DPB1 $INPUT6 >> $INPUT6.HLA_DPB1
grep HLA_DQA1 $INPUT6 >> $INPUT6.HLA_DQA1
grep HLA_DQB1 $INPUT6 >> $INPUT6.HLA_DQB1
grep HLA_DRA $INPUT6 >> $INPUT6.HLA_DRA
grep HLA_DRB1 $INPUT6 >> $INPUT6.HLA_DRB1
grep HLA_DRB3 $INPUT6 >> $INPUT6.HLA_DRB3
grep HLA_DRB4 $INPUT6 >> $INPUT6.HLA_DRB4
grep HLA_DRB5 $INPUT6 >> $INPUT6.HLA_DRB5

grep HLA_A $INPUT7 >> $INPUT7.HLA_A
grep HLA_B $INPUT7 >> $INPUT7.HLA_B
grep HLA_C $INPUT7 >> $INPUT7.HLA_C
grep HLA_E $INPUT7 >> $INPUT7.HLA_E
grep HLA_F $INPUT7 >> $INPUT7.HLA_F
grep HLA_G $INPUT7 >> $INPUT7.HLA_G
grep HLA_DMA $INPUT7 >> $INPUT7.HLA_DMA
grep HLA_DMB $INPUT7 >> $INPUT7.HLA_DMB
grep HLA_DOA $INPUT7 >> $INPUT7.HLA_DOA
grep HLA_DOB $INPUT7 >> $INPUT7.HLA_DOB
grep HLA_DPA1 $INPUT7 >> $INPUT7.HLA_DPA1
grep HLA_DPB1 $INPUT7 >> $INPUT7.HLA_DPB1
grep HLA_DQA1 $INPUT7 >> $INPUT7.HLA_DQA1
grep HLA_DQB1 $INPUT7 >> $INPUT7.HLA_DQB1
grep HLA_DRA $INPUT7 >> $INPUT7.HLA_DRA
grep HLA_DRB1 $INPUT7 >> $INPUT7.HLA_DRB1
grep HLA_DRB3 $INPUT7 >> $INPUT7.HLA_DRB3
grep HLA_DRB4 $INPUT7 >> $INPUT7.HLA_DRB4
grep HLA_DRB5 $INPUT7 >> $INPUT7.HLA_DRB5

grep HLA_A $INPUT8 >> $INPUT8.HLA_A
grep HLA_B $INPUT8 >> $INPUT8.HLA_B
grep HLA_C $INPUT8 >> $INPUT8.HLA_C
grep HLA_E $INPUT8 >> $INPUT8.HLA_E
grep HLA_F $INPUT8 >> $INPUT8.HLA_F
grep HLA_G $INPUT8 >> $INPUT8.HLA_G
grep HLA_DMA $INPUT8 >> $INPUT8.HLA_DMA
grep HLA_DMB $INPUT8 >> $INPUT8.HLA_DMB
grep HLA_DOA $INPUT8 >> $INPUT8.HLA_DOA
grep HLA_DOB $INPUT8 >> $INPUT8.HLA_DOB
grep HLA_DPA1 $INPUT8 >> $INPUT8.HLA_DPA1
grep HLA_DPB1 $INPUT8 >> $INPUT8.HLA_DPB1
grep HLA_DQA1 $INPUT8 >> $INPUT8.HLA_DQA1
grep HLA_DQB1 $INPUT8 >> $INPUT8.HLA_DQB1
grep HLA_DRA $INPUT8 >> $INPUT8.HLA_DRA
grep HLA_DRB1 $INPUT8 >> $INPUT8.HLA_DRB1
grep HLA_DRB3 $INPUT8 >> $INPUT8.HLA_DRB3
grep HLA_DRB4 $INPUT8 >> $INPUT8.HLA_DRB4
grep HLA_DRB5 $INPUT8 >> $INPUT8.HLA_DRB5

grep HLA_A $INPUT9 >> $INPUT9.HLA_A
grep HLA_B $INPUT9 >> $INPUT9.HLA_B
grep HLA_C $INPUT9 >> $INPUT9.HLA_C
grep HLA_E $INPUT9 >> $INPUT9.HLA_E
grep HLA_F $INPUT9 >> $INPUT9.HLA_F
grep HLA_G $INPUT9 >> $INPUT9.HLA_G
grep HLA_DMA $INPUT9 >> $INPUT9.HLA_DMA
grep HLA_DMB $INPUT9 >> $INPUT9.HLA_DMB
grep HLA_DOA $INPUT9 >> $INPUT9.HLA_DOA
grep HLA_DOB $INPUT9 >> $INPUT9.HLA_DOB
grep HLA_DPA1 $INPUT9 >> $INPUT9.HLA_DPA1
grep HLA_DPB1 $INPUT9 >> $INPUT9.HLA_DPB1
grep HLA_DQA1 $INPUT9 >> $INPUT9.HLA_DQA1
grep HLA_DQB1 $INPUT9 >> $INPUT9.HLA_DQB1
grep HLA_DRA $INPUT9 >> $INPUT9.HLA_DRA
grep HLA_DRB1 $INPUT9 >> $INPUT9.HLA_DRB1
grep HLA_DRB3 $INPUT9 >> $INPUT9.HLA_DRB3
grep HLA_DRB4 $INPUT9 >> $INPUT9.HLA_DRB4
grep HLA_DRB5 $INPUT9 >> $INPUT9.HLA_DRB5


# echo "Main Calling."
Rscript src/9GP_no_CI.R $INPUT.HLA_A $INPUT2.HLA_A $INPUT3.HLA_A $INPUT4.HLA_A $INPUT5.HLA_A $INPUT6.HLA_A $INPUT7.HLA_A $INPUT8.HLA_A $INPUT9.HLA_A A
Rscript src/9GP_no_CI.R $INPUT.HLA_B $INPUT2.HLA_B $INPUT3.HLA_B $INPUT4.HLA_B $INPUT5.HLA_B $INPUT6.HLA_B $INPUT7.HLA_B $INPUT8.HLA_B $INPUT9.HLA_B B
Rscript src/9GP_no_CI.R $INPUT.HLA_C $INPUT2.HLA_C $INPUT3.HLA_C $INPUT4.HLA_C $INPUT5.HLA_C $INPUT6.HLA_C $INPUT7.HLA_C $INPUT8.HLA_C $INPUT9.HLA_C C
Rscript src/9GP_no_CI.R $INPUT.HLA_E $INPUT2.HLA_E $INPUT3.HLA_E $INPUT4.HLA_E $INPUT5.HLA_E $INPUT6.HLA_E $INPUT7.HLA_E $INPUT8.HLA_E $INPUT9.HLA_E E
Rscript src/9GP_no_CI.R $INPUT.HLA_F $INPUT2.HLA_F $INPUT3.HLA_F $INPUT4.HLA_F $INPUT5.HLA_F $INPUT6.HLA_F $INPUT7.HLA_F $INPUT8.HLA_F $INPUT9.HLA_F F
Rscript src/9GP_no_CI.R $INPUT.HLA_G $INPUT2.HLA_G $INPUT3.HLA_G $INPUT4.HLA_G $INPUT5.HLA_G $INPUT6.HLA_G $INPUT7.HLA_G $INPUT8.HLA_G $INPUT9.HLA_G G
Rscript src/9GP_no_CI.R $INPUT.HLA_DMA $INPUT2.HLA_DMA $INPUT3.HLA_DMA $INPUT4.HLA_DMA $INPUT5.HLA_DMA $INPUT6.HLA_DMA $INPUT7.HLA_DMA $INPUT8.HLA_DMA $INPUT9.HLA_DMA DMA
Rscript src/9GP_no_CI.R $INPUT.HLA_DMB $INPUT2.HLA_DMB $INPUT3.HLA_DMB $INPUT4.HLA_DMB $INPUT5.HLA_DMB $INPUT6.HLA_DMB $INPUT7.HLA_DMB $INPUT8.HLA_DMB $INPUT9.HLA_DMB DMB
Rscript src/9GP_no_CI.R $INPUT.HLA_DOA $INPUT2.HLA_DOA $INPUT3.HLA_DOA $INPUT4.HLA_DOA $INPUT5.HLA_DOA $INPUT6.HLA_DOA $INPUT7.HLA_DOA $INPUT8.HLA_DOA $INPUT9.HLA_DOA DOA
Rscript src/9GP_no_CI.R $INPUT.HLA_DOB $INPUT2.HLA_DOB $INPUT3.HLA_DOB $INPUT4.HLA_DOB $INPUT5.HLA_DOB $INPUT6.HLA_DOB $INPUT7.HLA_DOB $INPUT8.HLA_DOB $INPUT9.HLA_DOB DOB
Rscript src/9GP_no_CI.R $INPUT.HLA_DPA1 $INPUT2.HLA_DPA1 $INPUT3.HLA_DPA1 $INPUT4.HLA_DPA1 $INPUT5.HLA_DPA1 $INPUT6.HLA_DPA1 $INPUT7.HLA_DPA1 $INPUT8.HLA_DPA1 $INPUT9.HLA_DPA1 DPA1
Rscript src/9GP_no_CI.R $INPUT.HLA_DPB1 $INPUT2.HLA_DPB1 $INPUT3.HLA_DPB1 $INPUT4.HLA_DPB1 $INPUT5.HLA_DPB1 $INPUT6.HLA_DPB1 $INPUT7.HLA_DPB1 $INPUT8.HLA_DPB1 $INPUT9.HLA_DPB1 DPB1
Rscript src/9GP_no_CI.R $INPUT.HLA_DQA1 $INPUT2.HLA_DQA1 $INPUT3.HLA_DQA1 $INPUT4.HLA_DQA1 $INPUT5.HLA_DQA1 $INPUT6.HLA_DQA1 $INPUT7.HLA_DQA1 $INPUT8.HLA_DQA1 $INPUT9.HLA_DQA1 DQA1
Rscript src/9GP_no_CI.R $INPUT.HLA_DQB1 $INPUT2.HLA_DQB1 $INPUT3.HLA_DQB1 $INPUT4.HLA_DQB1 $INPUT5.HLA_DQB1 $INPUT6.HLA_DQB1 $INPUT7.HLA_DQB1 $INPUT8.HLA_DQB1 $INPUT9.HLA_DQB1 DQB1
Rscript src/9GP_no_CI.R $INPUT.HLA_DRA $INPUT2.HLA_DRA $INPUT3.HLA_DRA $INPUT4.HLA_DRA $INPUT5.HLA_DRA $INPUT6.HLA_DRA $INPUT7.HLA_DRA $INPUT8.HLA_DRA $INPUT9.HLA_DRA DRA
Rscript src/9GP_no_CI.R $INPUT.HLA_DRB1 $INPUT2.HLA_DRB1 $INPUT3.HLA_DRB1 $INPUT4.HLA_DRB1 $INPUT5.HLA_DRB1 $INPUT6.HLA_DRB1 $INPUT7.HLA_DRB1 $INPUT8.HLA_DRB1 $INPUT9.HLA_DRB1 DRB1
Rscript src/9GP_no_CI.R $INPUT.HLA_DRB3 $INPUT2.HLA_DRB3 $INPUT3.HLA_DRB3 $INPUT4.HLA_DRB3 $INPUT5.HLA_DRB3 $INPUT6.HLA_DRB3 $INPUT7.HLA_DRB3 $INPUT8.HLA_DRB3 $INPUT9.HLA_DRB3 DRB3
Rscript src/9GP_no_CI.R $INPUT.HLA_DRB4 $INPUT2.HLA_DRB4 $INPUT3.HLA_DRB4 $INPUT4.HLA_DRB4 $INPUT5.HLA_DRB4 $INPUT6.HLA_DRB4 $INPUT7.HLA_DRB4 $INPUT8.HLA_DRB4 $INPUT9.HLA_DRB4 DRB4
Rscript src/9GP_no_CI.R $INPUT.HLA_DRB5 $INPUT2.HLA_DRB5 $INPUT3.HLA_DRB5 $INPUT4.HLA_DRB5 $INPUT5.HLA_DRB5 $INPUT6.HLA_DRB5 $INPUT7.HLA_DRB5 $INPUT8.HLA_DRB5 $INPUT9.HLA_DRB5 DRB5

### Merge the outputs from '9GP_no.R'
# cat *.alleles > $INPUT.9.merged

#  cat ${INPUT}.HLA_A.alleles \
#      ${INPUT}.HLA_B.alleles \
#      ${INPUT}.HLA_C.alleles \
#      ${INPUT}.HLA_DRB1.alleles \
#      ${INPUT}.HLA_DQA1.alleles \
#      ${INPUT}.HLA_DQB1.alleles \
#      ${INPUT}.HLA_DPA1.alleles \
#      ${INPUT}.HLA_DPB1.alleles > ${out}.alleles 
