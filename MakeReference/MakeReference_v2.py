# -*- coding: utf-8 -*-

import os, sys, re
import subprocess

#from MakeReference.src.HLAtoSequences import HLAtoSequences
#from MakeReference.src.encodeVariants import encodeVariants
from MakeReference.src.encodeHLA import encodeHLA

from MakeReference.src.ATtrick import ATtrick
from MakeReference.src.redefineBPv1BH import redefineBP

########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def MakeReference_v2(_CHPED, _OUT, _hg, _genes="A,B,C,E,F,G,H,J,K,L,V,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB4,DRB5,MICA,MICB", _variants=None,
                     _java_heap_mem='2g', _java_stack_mem='1g', _java_tmp_folder="/tmp",
                     _p_dependency="dependency/", f_save_intermediates=True, f_phasing=False, _nthreads=1,
                     _burnin=3, _iter=12, _map="null",
                     _mind=0.3, _hardy=0.000001, _maf=0.01, _miss=0.05,
                     _hla_maf=0.0001, _window=2.8, _overlap=1.25):


    ########## < Core Variables > ##########

    ### [1] Major Path Variables
    if os.path.exists(_p_dependency):
        p_dependency = _p_dependency
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Folder for dependent software('{}') can't be found. Please check '--dependency' argument again.".format(_p_dependency))
        sys.exit()

    # dependent software.
    _p_plink = os.path.join(p_dependency, "plink")
    _p_beagle = os.path.join(p_dependency, "beagle.jar")
    _p_linkage2beagle = os.path.join(p_dependency, "linkage2beagle.jar")
    _p_beagle2vcf = os.path.join(p_dependency, "beagle2vcf.jar")

    if not os.path.exists(_p_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(_p_dependency))
        sys.exit()

    if not os.path.exists(_p_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'Beagle 5.4' (https://faculty.washington.edu/browning/beagle/b4_1.html#download) in '{0}'\n".format(_p_dependency))
        sys.exit()

    if not os.path.exists(_p_linkage2beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'linkage2beagle.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) in '{0}'\n".format(_p_dependency))
        sys.exit()

    if not os.path.exists(_p_beagle2vcf):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'beagle2vcf.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) in '{0}'\n".format(_p_dependency))
        sys.exit()



    ### [2] Memory representation check.
    p_Mb = re.compile(r'\d+m')
    p_Gb = re.compile(r'\d+[gG]')

    if not (bool(p_Mb.match(_java_heap_mem)) or bool(p_Gb.match(_java_heap_mem))):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given Java memory value('{}') has bizzare representation. Please check '--mem' argument again.".format(_java_heap_mem))
        sys.exit()


# 
# ### [3] Dictionary Files
# 
# _dictionary_AA_seq = _dictionary_AA + ".txt" # From now on, official extension of HLA sequence information dictionary is ".txt". (2018. 9. 25.)
# _dictionary_AA_map = _dictionary_AA + ".map"
# 
# _dictionary_SNPS_seq = _dictionary_SNPS + ".txt"
# _dictionary_SNPS_map = _dictionary_SNPS + ".map"
# 
# 
# if not os.path.exists(_dictionary_AA_map):
#   print(std_ERROR_MAIN_PROCESS_NAME + "AA dictionary map file can't be found('{0}'). Please check '--dict-AA' argument again.\n".format(_dictionary_AA_map))
# sys.exit()
# 
# if not os.path.exists(_dictionary_AA_seq):
#   print(std_ERROR_MAIN_PROCESS_NAME + "AA dictionary seq file can't be found('{0}'). Please check '--dict-AA' argument again.\n".format(_dictionary_AA_seq))
# sys.exit()
# 
# if not os.path.exists(_dictionary_SNPS_map):
#   print(std_ERROR_MAIN_PROCESS_NAME + "SNPS dictionary map file can't be found('{0}'). Please check '--dict-SNPS' argument again.\n".format(_dictionary_SNPS_map))
# sys.exit()
# 
# if not os.path.exists(_dictionary_SNPS_seq):
#   print(std_ERROR_MAIN_PROCESS_NAME + "SNPS dictionary seq file can't be found('{0}'). Please check '--dict-SNPS' argument again.\n".format(_dictionary_SNPS_seq))
# sys.exit()



    ### Intermediate path.
    OUTPUT = _OUT if not _OUT.endswith('/') else _OUT.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)



    ########## < Core Variables 2 > ##########

    ### Flag for plain SNP markers.
    f_plain_SNP = bool(_variants)

    # Input 1 : HLA type data
    HLA_DATA = _CHPED

    if not os.path.exists(_CHPED):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA type can't be found('{}'). Please check '--chped' argument again.".format(_CHPED))
        sys.exit()

    # Input 2 : Plain SNP data
    if f_plain_SNP:

        if not os.path.exists(_variants+'.bed'):
            print(std_ERROR_MAIN_PROCESS_NAME + "One of SNP data can't be found('{}'). Please check '--variants' argument again.".format(_variants+'.bed'))
            sys.exit()
        if not os.path.exists(_variants+'.bim'):
            print(std_ERROR_MAIN_PROCESS_NAME + "One of SNP data can't be found('{}'). Please check '--variants' argument again.".format(_variants+'.bim'))
            sys.exit()
        if not os.path.exists(_variants+'.fam'):
            print(std_ERROR_MAIN_PROCESS_NAME + "One of SNP data can't be found('{}'). Please check '--variants' argument again.".format(_variants+'.fam'))
            sys.exit()

        SNP_DATA = _variants
        SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, os.path.basename(_variants))

    else:
        print(std_WARNING_MAIN_PROCESS_NAME + "SNP data wasn't given. MakeReference_v2 will generate a reference panel only with HLA type information.")


    # Output prefix
    #AA_CODED = OUTPUT + '.AA.CODED'
    HLA_CODED = OUTPUT + ".HLA"
    #SNPS_CODED = OUTPUT + '.SNPS.CODED'

    plink = ' '.join([_p_plink, "--noweb", "--silent"])
    beagle = ' '.join(["java", "-Xmx{}".format(_java_heap_mem), "-Xss{}".format(_java_stack_mem), "-Djava.io.tmpdir={}".format(_java_tmp_folder), "-jar", _p_beagle])
    linkage2beagle = ' '.join(["java", "-Xmx{}".format(_java_heap_mem), "-Djava.io.tmpdir={}".format(_java_tmp_folder), "-jar", _p_linkage2beagle])
    beagle2vcf = ' '.join(["java", "-Xmx{}".format(_java_heap_mem), "-Djava.io.tmpdir={}".format(_java_tmp_folder), "-jar", _p_beagle2vcf])




    ########## <Flags for Code Block> ##########

    #ENCODE_AA = 1
    ENCODE_HLA = 1
    #ENCODE_SNPS = 1

    EXTRACT_FOUNDERS = 1
    MERGE = 1
    #QC = 1

    PREPARE = 1
    PHASE = 1
    CLEANUP = 1

    # (2019. 01. 10.) Last three code blocks won't be implemented



    ########## <Making Reference Panel> ##########

    print(std_MAIN_PROCESS_NAME + "Making Reference Panel for \"{0}\"".format(OUTPUT))
    index = 1


    if ENCODE_HLA:

        print("[{}] Encoding HLA alleles.".format(index))

        ### (1) Encoded HLA ( *.HLA.{ped,map} ) ###
        encodeHLA(HLA_DATA, OUTPUT + ".HLA", _genes, _hg)

        ### (2) Final Encoded Outputs ( *.HLA.{bed,bim,fam,nosex,log} ) ###
        command = ' '.join([plink, "--file", OUTPUT + '.HLA', "--make-bed", "--out", OUTPUT + '.HLA'])
        # print(command)
        os.system(command)

        index += 1


        """
        Generaed outputs :
            (1) Encoded HLA ( *.HLA.{ped,map} )
            (2) Final Encoded Outputs ( *.HLA.{bed,bim,fam}, *.HLA.{nosex,log} )
            
        Final outputs :
            - *.HLA.{bed,bim,fam,factors,nosex,log}
            
        Outputs to remove :
            - *.HLA.{ped,map}
        """

        if not f_save_intermediates:
            os.system("rm " + (OUTPUT + ".HLA.ped"))
            os.system("rm " + (OUTPUT + ".HLA.map"))



        if EXTRACT_FOUNDERS:

            print("[{}] Extracting founders.".format(index))

            """
            if ($EXTRACT_FOUNDERS) then
                echo "[$i] Extracting founders."; @ i++
                plink --bfile $SNP_DATA --filter-founders --mind 0.3 --alleleACGT --make-bed --out $SNP_DATA.FOUNDERS
            
                # Initial QC on Reference SNP panel
                plink --bfile $SNP_DATA.FOUNDERS --hardy        --out $SNP_DATA.FOUNDERS.hardy  # ?? 92?? ?? position?? HWE test? ??
                plink --bfile $SNP_DATA.FOUNDERS --freq         --out $SNP_DATA.FOUNDERS.freq   # ?? --freq ??? allele frequency????? ???. 
                plink --bfile $SNP_DATA.FOUNDERS --missing      --out $SNP_DATA.FOUNDERS.missing
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe      | awk ' $9 < 0.000001 { print $2 }' | sort -u > remove.snps.hardy 
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.freq.frq       | awk ' $5 < 0.01 { print $2 } '             > remove.snps.freq
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.missing.lmiss  | awk ' $5 > 0.05 { print $2 } '              > remove.snps.missing
                cat remove.snps.*                                            | sort -u                                     > all.remove.snps
            
                plink --bfile $SNP_DATA.FOUNDERS --allow-no-sex --exclude all.remove.snps --make-bed --out $SNP_DATA.FOUNDERS.QC
            
                # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file
            
                plink --bfile $OUTPUT.HLA --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.HLA.FOUNDERS
                plink --bfile $OUTPUT.SNPS.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.SNPS.FOUNDERS
                plink --bfile $OUTPUT.AA.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.AA.FOUNDERS
            
                rm remove.snps.*
            endif
            """

            ### (1) --filter-founders to variants data ( *.FOUNDERS )
            command = ' '.join([plink, "--bfile", SNP_DATA, "--filter-founders", "--mind", _mind, "--alleleACGT", "--make-bed", "--out", SNP_DATA2+'.FOUNDERS'])
            # print(command)
            os.system(command)


            ### (2) QC ( *.FOUNDERS.{hardy,freq,missing )
            # Initial QC on Reference SNP panel
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--hardy", "--out", SNP_DATA2+'.FOUNDERS.hardy'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--freq", "--out", SNP_DATA2+'.FOUNDERS.freq'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--missing", "--out", SNP_DATA2+'.FOUNDERS.missing'])
            # print(command)
            os.system(command)


            ### (3) Stuffs to remove ( remove.snps.hardy, remove.snps.freq, remove.snps.missing, all.remove.snps )
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.hardy.hwe', "|", "awk", "' $9 <", _hardy, "{ print $2 }'", "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.hardy")])
            # print(command)
            os.system(command)
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.freq.frq', "|", "awk", "' $5 <", _maf, "{ print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.freq")])
            # print(command)
            os.system(command)
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.missing.lmiss', "|", "awk", "' $5 >", _miss, "{ print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.missing")])
            # print(command)
            os.system(command)
            command = ' '.join(["cat", os.path.join(INTERMEDIATE_PATH, "remove.snps.*"), "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "all.remove.snps")])
            # print(command)
            os.system(command)


            ### (4) Filtering out Quality-controled FOUNDERS ( *.FOUNDERS.QC )
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--allow-no-sex", "--exclude", os.path.join(INTERMEDIATE_PATH, "all.remove.snps"), "--make-bed", "--out", SNP_DATA2+'.FOUNDERS.QC'])
            # print(command)
            os.system(command)

            # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

            ### (5) --filter-founders to HLA information ( *.{HLA,AA,SNPS}.FOUNDERS )
            command = ' '.join([plink, "--bfile", OUTPUT+'.HLA', "--filter-founders", "--maf", _hla_maf, "--make-bed", "--out", OUTPUT+'.HLA.FOUNDERS'])
            # print(command)
            os.system(command)
            #command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.SNPS.FOUNDERS'])
            # print(command)
            #os.system(command)
            #command = ' '.join([plink, "--bfile", OUTPUT+'.AA.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.AA.FOUNDERS'])
            # print(command)
            #os.system(command)


            """
            Generated outputs :
                (1) --filter-founders to variants data ( *.FOUNDERS )
                (2) QC ( *.FOUNDERS.{hardy,freq,missing )
                (3) Stuffs to remove ( remove.snps.hardy, remove.snps.freq, remove.snps.missing, all.remove.snps )
                (4) Filtering out Quality-controled FOUNDERS ( *.FOUNDERS.QC )
                (5) --filter-founders to HLA information ( *.{HLA,AA,SNPS}.FOUNDERS )
                
            Final outputs :
                - *.FOUNDERS.QC
                - *.{HLA,AA,SNPS}.FOUNDERS
                
            Outputs to remove :
                - *.FOUNDERS
                - *.FOUNDERS.{hardy,freq,missing}
                - remove.snps.hardy, remove.snps.freq, remove.snps.missing, all.remove.snps
                
            """


            if not f_save_intermediates:

                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.bed"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.bim"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.fam"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.log"))
                if os.path.exists(SNP_DATA2 + ".FOUNDERS.nosex"):
                    os.system("rm " + (SNP_DATA2 + ".FOUNDERS.nosex"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.hardy.*"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.freq.*"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.missing.*"))
                os.system("rm " + os.path.join(INTERMEDIATE_PATH, "remove.snps.*"))
                os.system("rm " + os.path.join(INTERMEDIATE_PATH, "all.remove.snps"))


            index += 1


        if MERGE:

            print("[{}] Merging SNP, HLA, and amino acid datasets.".format(index))

            """
            echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
            echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
            echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
            echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
            plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
            rm $OUTPUT.HLA.???
            rm $OUTPUT.AA.CODED.???
            rm $OUTPUT.SNPS.CODED.???
            rm merge_list
    
            """

            ### (1) Stuffs to merge ( merge_list )
            TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")


            command = ' '.join(["echo", OUTPUT + '.HLA.FOUNDERS.bed', OUTPUT + '.HLA.FOUNDERS.bim', OUTPUT + '.HLA.FOUNDERS.fam', ">", TMP_merged_list])
            # print(command)
            os.system(command)

            #command = ' '.join(["echo", OUTPUT + '.AA.FOUNDERS.bed', OUTPUT + '.AA.FOUNDERS.bim', OUTPUT + '.AA.FOUNDERS.fam', ">>", TMP_merged_list])
            # print(command)
            #os.system(command)

            #command = ' '.join(["echo", OUTPUT + '.SNPS.FOUNDERS.bed', OUTPUT + '.SNPS.FOUNDERS.bim', OUTPUT + '.SNPS.FOUNDERS.fam', ">>", TMP_merged_list])
            # print(command)
            #os.system(command)


            ### (2) Merging the above stuffs ( *.MERGED.FOUNDERS )
            command = ' '.join(
                [plink, "--bfile", SNP_DATA2 + '.FOUNDERS.QC', "--merge-list", TMP_merged_list, "--make-bed", "--out",
                 OUTPUT + '.MERGED.FOUNDERS'])
            # print(command)
            os.system(command)


            """
            Generated Outputs :
                (1) Stuffs to merge ( merge_list )
                (2) Merging the above stuffs ( *.MERGED.FOUNDERS )
            
            Final outputs : 
                - *.MERGED.FOUNDERS
            
            Outputs to remove :
                - *.{AA,HLA,SNPS}.FOUNDERS.{bed,bim,fam}
                - merge_list
            
            
            """

            if not f_save_intermediates:

                os.system("rm " + (OUTPUT + ".HLA.*"))
                #os.system("rm " + (OUTPUT + ".AA.*"))
                #os.system("rm " + (OUTPUT + ".SNPS.*"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.QC.*"))
                os.system("rm " + TMP_merged_list)



            index += 1
        
        
        
        """
            awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order
        
            plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --make-bed --out $OUTPUT
        
            # Calculate allele frequencies
            plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
            rm $SNP_DATA.FOUNDERS.*
            rm $OUTPUT.MERGED.FOUNDERS.*
            rm $OUTPUT.*.FOUNDERS.???
            rm allele.order
            rm all.remove.snps
        
        """
        
        TMP_allele_order = OUTPUT + ".refallele"
        
        
        command = ' '.join([plink, "--bfile", OUTPUT + '.MERGED.FOUNDERS', "--freq", "--out", OUTPUT + '.MERGED.FOUNDERS.FRQ'])
            # print(command)
        os.system(command)
        
        command = ' '.join(
            ["awk", '\'{if (NR > 1){if (($3 == "a" && $4 == "p") || ($4 == "a" && $3 == "p")){print $2 "\tp"}}}\'',
             OUTPUT + '.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
        # print(command)
        os.system(command)
        
        command = ' '.join(
            [plink, "--bfile", OUTPUT + '.MERGED.FOUNDERS', "--a1-allele", TMP_allele_order, "--make-bed", "--out", OUTPUT])
        # print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--freq", "--out", OUTPUT + '.FRQ',
                            "--a1-allele", TMP_allele_order])
        # print(command)
        os.system(command)
        
        """
            Generated Outputs :
                (1) Frequency file to use in filtering out some snps ( *.MERGED.FOUNDERS.FRQ )
                (2) List up SNPs which have extreme allele frequency. ( all.remove.snps )
                (3) Making reference allele ( allele.order )
                (4) Filtering out SNPs which have extreme allele frequency(The reference panel for association test.) ( *.{bed,bim,fam} )
                (5) Allele frequency info. of output reference panel ( *.FRQ )
        
            Final outputs : 
                - *.{bed,bim,fam}
                - *.FRQ
        
            Outputs to remove :
                - *.MERGED.FOUNDERS.FRQ
                - all.remove.snps
                - allele.order
                
        """
        
        if not f_save_intermediates:
        
            os.system("rm " + (OUTPUT + ".MERGED.FOUNDERS.*"))
            os.system("rm " + (OUTPUT + ".FRQ.log"))
            os.system("rm " + (OUTPUT + ".log"))
            os.system("rm " + TMP_allele_order)
        
            if os.path.exists(OUTPUT + ".FRQ.nosex"):
                os.system("rm " + (OUTPUT + ".FRQ.nosex"))
            if os.path.exists(OUTPUT + ".nosex"):
                os.system("rm " + (OUTPUT + ".nosex"))
        
        index += 1

    
        
        
        
        if f_phasing:


            if PREPARE:

                print("[{}] Preparing files for Beagle.".format(index))

                """
                [Source from Buhm Han.]

                awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
                plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
                awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
                cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped

                echo "[$i] Converting to beagle format.";  @ i++
                linkage2beagle pedigree=$OUTPUT.nopheno.ped data=$OUTPUT.dat beagle=$OUTPUT.bgl standard=true > $OUTPUT.bgl.log


                [Source from Yang.]

                awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
                plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
                plink --bfile $OUTPUT --recode --transpose --out $OUTPUT
                # awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
                # cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped

                echo "[$i] Converting to beagle format.";  @ i++
                beagle2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnout $OUTPUT.vcf

                I will make this code block based on source given by Yang. for now.

                """

                # ATtrick : 'P' to 'T', 'A' to 'A'
                [bim_ATtrick, a1_allele_ATtrick] = ATtrick(OUTPUT + '.bim', OUTPUT)


                # command = ' '.join(["awk", '\'{print $2" "$4" "$5" "$6}\'', OUTPUT + '.bim', ">", OUTPUT + '.markers'])
                command = ' '.join(["awk", '\'{print $2" "$4" "$6" "$5}\'', bim_ATtrick, ">", OUTPUT + '.ATtrick.markers'])
                # print(command)
                os.system(command)
                """
                Plink works by setting ALT allele as a1-allele, which is the 5th column of bim file.
                However, VCF file sets a2-allele, which is the 6th column of plink bim file, as ALT allele.

                beagle2vcf converts 4th column of *.markers file to ALT allele of newly generated vcf file.
                That's why I reordered ($5, $6) to ($6, $5).
                """

                # Manipulate duplicated Base poistion
                redefined_markers = redefineBP(OUTPUT + '.ATtrick.markers', OUTPUT + '.markers')

                # Applying the above manipulated base position information.
                # command = [
                #     "paste <(cut -f1-3 %s) <(awk '{print $2}' %s) <(cut -f5- %s) > %s" % (bim_ATtrick, redefined_markers, bim_ATtrick, OUTPUT+'.ATtrick.redefined.bim')
                # ]
                # print(command)
                # # os.system(command)
                # subprocess.call(command)

                command = ' '.join(
                    [plink, "--bed", OUTPUT+'.bed', "--bim", bim_ATtrick, "--fam", OUTPUT+'.fam',
                     "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT+'.ATtrick',
                     "--a1-allele", a1_allele_ATtrick])
                # print(command)
                os.system(command)

                command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT + '.ATtrick.map', ">", OUTPUT + '.ATtrick.dat'])
                # print(command)
                os.system(command)

                command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT + '.ATtrick.ped', ">", OUTPUT + '.ATtrick.nopheno.ped'])
                # print(command)
                os.system(command)

                index += 1

                print("[{}] Converting PLINK to BEAGLE format.".format(index))

                command = ' '.join([linkage2beagle, "pedigree=" + OUTPUT + '.ATtrick.nopheno.ped', "data=" + OUTPUT + '.ATtrick.dat',
                                    "beagle=" + OUTPUT + '.ATtrick.bgl', "standard=true", ">", OUTPUT + '.ATtrick.bgl.log'])
                # print(command)
                os.system(command)

                index += 1

                # for Beagle 4.1
                print("[{}] Converting BEAGLE to VCF format.".format(index))

                command = ' '.join([beagle2vcf, '6', redefined_markers, OUTPUT + '.ATtrick.bgl', '0', '>', OUTPUT+'.bgl.vcf'])
                # print(command)
                os.system(command)

                index += 1


                if not f_save_intermediates:
                    # os.system("rm {}".format())
                    os.system("rm {}".format(OUTPUT + '.ATtrick.bgl.log'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.bgl'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.log'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.dat'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.nopheno.ped'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.map'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.ped'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.markers'))
                    os.system("rm {}".format(bim_ATtrick))
                    os.system("rm {}".format(a1_allele_ATtrick))



            if PHASE:

                #print("[{}] Phasing reference using Beagle4.1.".format(index))
                print("[{}] Phasing reference using Beagle 5.4.".format(index))

                '''
                # Beagle v3.x.x
                beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log
    
                # Beagle v4.1
                beagle gt=OUTPUT+'.ATtrick.redefined.vcf' out=OUTPUT+'.ATtrick.redefined.phased' nthreads=1 impute=false niterations=10 lowmem=true >> $OUTPUT.bgl.log
    
                '''

                # command = ' '.join([beagle, "gt={}".format(OUTPUT+'.bglv4.bgl.vcf'),
                #                     "nthreads=1", "impute=false",
                #                     "niterations=10", "lowmem=true", "out={}".format(OUTPUT+'.bglv4.bgl.phased'),
                #                     '>', OUTPUT+'.bglv4.bgl.phased.vcf.log'])

                #command = ' '.join([beagle, "gt={}".format(OUTPUT+'.bgl.vcf'),
                #                   "impute=false", "nthreads={}".format(_nthreads),
                #                    "niterations=5", "lowmem=true", "out={}".format(OUTPUT+'.bgl.phased')])
                
                if _map=="null":                
                   
                  command = ' '.join([beagle, "gt={}".format(OUTPUT+'.bgl.vcf'),
                                      "impute=true", "nthreads={}".format(_nthreads),
                                      "burnin={}".format(_burnin), "iterations={}".format(_iter),
                                      "out={}".format(OUTPUT+'.bgl.phased'), "window={}".format(_window), "overlap={}".format(_overlap)])
                
                else:
                
                  command = ' '.join([beagle, "gt={}".format(OUTPUT+'.bgl.vcf'),
                                      "impute=true", "nthreads={}".format(_nthreads),
                                      "burnin={}".format(_burnin), "iterations={}".format(_iter),
                                      "out={}".format(OUTPUT+'.bgl.phased'), "map={}".format(_map), "window={}".format(_window), "overlap={}".format(_overlap)])
                                    
                # print(command)

                try:
                    #os.system(command)
                    f_log = open(OUTPUT+'.bgl.phased.vcf.log', 'w')
                    subprocess.run(re.split(r'\s+',command), check=True, stdout=f_log, stderr=f_log)

                except subprocess.CalledProcessError:
                    # fail.
                    print(std_ERROR_MAIN_PROCESS_NAME + "Phasing failed. See log file('{}').".format(OUTPUT+'.bgl.phased.vcf.log'))
                    sys.exit()
                else:
                    # succeed.
                    f_log.close()
                    if not f_save_intermediates:
                        # os.system("rm {}".format())
                        os.system("rm {}".format(OUTPUT + '.bgl.vcf'))

                        # remove redundant log file.
                        os.system("rm {}".format(OUTPUT+'.bgl.phased.log'))

                index += 1



        if CLEANUP:

            print("[{}] Removing unnecessary files, and munging output for rest of pipeline".format(index))

            '''
            rm $OUTPUT.nopheno.ped
            rm $OUTPUT.bgl.gprobs
            rm $OUTPUT.bgl.r2
            rm $OUTPUT.bgl
            rm $OUTPUT.ped
            rm $OUTPUT.map
            rm $OUTPUT.dat
            rm $OUTPUT.phasing.log
            '''

            #rm_tlist = ('.nopheno.ped', '.bgl.gprobs', '.bgl.r2', '.bgl', '.ped', '.map', '.dat')

            #for i in rm_tlist:
            #    print("rm " + OUTPUT + i)
            #    os.system("rm " + OUTPUT + i)

            index += 1
            
            #command = ''.join(['mv ', OUTPUT, '.bgl.phased.vcf.gz ', OUTPUT, '_tmp.bgl.phased.vcf.gz'])
            #os.system(command)
            
            command = ''.join(["bash MakeReference/src/Adjust_output.sh --path_Ref ", OUTPUT])
            os.system(command)
            
            #command = ''.join(['rm ', OUTPUT, '_tmp.bgl.phased.vcf.gz'])
            #os.system(command)

        print("[{}] Making reference panel for HLA and is Done!".format(index))

        __return__ = OUTPUT

    return _OUT


