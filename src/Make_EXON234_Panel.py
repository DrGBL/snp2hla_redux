#-*- coding: utf-8 -*-

import os, sys, re
from src.redefineBPv1BH import redefineBP
from src.BGL2SortBGl import BGL2SortBGL_WS
from math import floor

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "E", "F", "G", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5"]
isClassI = {"A": True, "B": True, "C": True, "E": True, "F": True, "G": True,
            "DMA": False, "DMB": False, "DOA": False,
            "DPA1": False, "DPB1": False, "DQA1": False, "DQB1": False,
            "DRA": False, "DRB1": False, "DRB3": False, "DRB4": False, "DRB5": False}

# Patterns to use.
# p_HLA_2field = re.compile(r'^HLA_(\w+)_\d{4}')
# p = re.compile(r'^([A-Za-z0-9_-]+)\s+(\w+)') # Frist two columns (ex. 'P pedigree' or 'rs969931 29602876', ... )
# p = re.compile(r'^(\S+)\s+(\S+)') # Frist two columns (ex. 'P pedigree' or 'rs969931 29602876', ... )


def Make_EXON234_Panel(infile, outfile, BEAGLE2LINKAGE, PLINK, __save_intermediates=False):


    # REF_base = os.path.basename(outfile)
    # OUTPUT_dir = os.path.dirname(outfile)
    # outfile = os.path.join(OUTPUT_dir, REF_base)



    ### STEP1_Collect_SNP_HLA_DATA

    # # In STEP1, New *.markers file will be used just next step.
    # command = "grep rs {} > {}".format(infile + ".markers", outfile+".STEP1_SNP.markers")
    # # print(command)
    # os.system(command)
    #
    # command = "grep \'HLA_[A-Z]_[0-9][0-9][0-9][0-9]\' {} > {}".format(infile + ".markers", outfile+".STEP1_class1_4dit.markers")
    # # print(command)
    # os.system(command)
    #
    # command = "grep \'HLA_[A-Z][A-Z][A-Z][0-9]_[0-9][0-9][0-9][0-9]\' {} > {}".format(infile + ".markers", outfile+".STEP1_class2_4dit.markers")
    # # print(command)
    # os.system(command)
    #
    # command = 'cat {} {} {} > {}'.format(outfile+".STEP1_SNP.markers", outfile+".STEP1_class1_4dit.markers",
    #                                      outfile+".STEP1_class2_4dit.markers", outfile+".STEP1_SNP_4dit.markers")
    # # print(command)
    # os.system(command)
    #
    #
    # # Remove
    # if not __save_intermediates:
    #     os.system('rm {}'.format(outfile+".STEP1_SNP.markers"))
    #     os.system('rm {}'.format(outfile+".STEP1_class1_4dit.markers"))
    #     os.system('rm {}'.format(outfile+".STEP1_class2_4dit.markers"))


    #p_MkRef_ToExclude = re.compile(r'^(AA_|SNP_|INS_|HLA_\w+_\d{2}$)')
    p_MkRef_ToExclude = re.compile(r'^(AA_|SNP_|INS_|HLA_\w+_\d{2}$)')

    with open(infile + ".markers", 'r') as f_in_markers, open(outfile+".STEP1_SNP_4dit.markers", 'w') as f_out_markers:
        for line in f_in_markers:
            l = line.split()

            m = p_MkRef_ToExclude.match(l[0])

            if not m:
                # To save
                f_out_markers.write(line)





    ### STEP2_EXON234_MARKERS

    [outbgl, outmarker] = HLA2EXON234(outfile+".STEP1_SNP_4dit.markers",
                                      infile + ".bgl.phased", outfile+".STEP2_exon234.bgl.phased",
                                      infile + ".markers", outfile+".STEP2_exon234.markers")

    # Remove
    if not __save_intermediates:
        os.system('rm {}'.format(outfile+".STEP1_SNP_4dit.markers"))


    ### STEP3_SORT

    # Dispersing genomic positions of given marker file (*.markers)
    refiend_outmarker = redefineBP(outmarker, outfile+".STEP3_refined.markers")
    # print(refiend_outmarker)


    # Sorting the dispersed marker file.
    command = 'sort -gk 2 {} > {}'.format(refiend_outmarker, outfile+'.markers')
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outmarker))
            os.system('rm {}'.format(refiend_outmarker))


    # Sorting beagle file to the oreder of the above sorted markers file
    sorted_outbgl = BGL2SortBGL_WS(outfile+'.markers', outbgl, outfile + ".bgl.phased")
    # print(sorted_outbgl)
    if not os.path.exists(sorted_outbgl):
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate '{}'.".format(sorted_outbgl))
        sys.exit()
    else:
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outbgl))



    ### STEP_4_Make_plink_file

    command = 'cat {} | {} {}'.format(sorted_outbgl, BEAGLE2LINKAGE, outfile + ".STEP4_tmp") # *.ped, *.dat (cf. 'java -jar' is included in 'BEAGLE2LINKAGE'.)
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".STEP4_tmp.dat")) # *.dat file is unnecessary.


    command = 'cut -d \' \' -f-5 {} > {}'.format(outfile + ".STEP4_tmp.ped", outfile + ".STEP4_tmp.ped.left") # ['FID', 'IID', 'PID', 'MID', 'Sex']
    # print(command)
    os.system(command)

    command = 'cut -d \' \' -f6- {} > {}'.format(outfile + ".STEP4_tmp.ped", outfile + ".STEP4_tmp.ped.right") # genotype information part.
    # print(command)
    os.system(command)


    command = 'paste -d \' -9 \' {} /dev/null /dev/null /dev/null {} > {}'.format(outfile + ".STEP4_tmp.ped.left", outfile + ".STEP4_tmp.ped.right", outfile + ".ped")
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".STEP4_tmp.ped"))
            os.system('rm {}'.format(outfile + ".STEP4_tmp.ped.left"))
            os.system('rm {}'.format(outfile + ".STEP4_tmp.ped.right"))


    # (1) rsid, (2) bp, (3) allele1
    os.system(' '.join(["cut -d \' \' -f1", outfile + ".markers", ">", outfile + ".STEP4_map.rsid"]))

    os.system(' '.join(["cut -d \' \' -f2", outfile + ".markers", ">", outfile + ".STEP4_map.bp"]))

    os.system(' '.join(["cut -d \' \' -f3", outfile + ".markers", ">", outfile + ".STEP4_map.allele1"]))


    os.system(' '.join(
        ["paste -d \'6  0 \'", "/dev/null", "/dev/null", outfile + ".STEP4_map.rsid", "/dev/null", "/dev/null",
         outfile + ".STEP4_map.bp", ">", outfile + ".map"]))

    # os.system(' '.join(
    #     ["paste -d \'   \'", outfile + ".STEP4_map.rsid", outfile + ".STEP4_map.bp", ">", outfile + ".refallele"]))

    os.system(' '.join(
        ["paste -d \' \'", outfile + ".STEP4_map.rsid", outfile + ".STEP4_map.allele1", ">", outfile + ".refallele"]))

    """
    (2019. 07. 09.)
    To make '*.refallele' file, I think right part is supposed to be 'outfile + ".STEP4_map.allele1"' not 'outfile + ".STEP4_map.bp"'
    """


    # bed, bim, fam files.
    command = ' '.join([PLINK, '--ped {} --map {} --make-bed --reference-allele {} --out {}'.format(
        outfile + ".ped",
        outfile + ".map",
        outfile + ".refallele",
        outfile
    )])
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".STEP4_map.rsid"))
            os.system('rm {}'.format(outfile + ".STEP4_map.bp"))
            os.system('rm {}'.format(outfile + ".STEP4_map.allele1"))
            os.system('rm {}'.format(outfile + ".ped"))
            os.system('rm {}'.format(outfile + ".map"))
            os.system('rm {}'.format(outfile + ".log"))
            os.system('rm {}'.format(outfile + ".refallele"))


    # Allele Frequency file(*.frq)
    command = ' '.join([PLINK, '--bfile {} --keep-allele-order --freq --out {}'.format(outfile, outfile + ".FRQ")])
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".FRQ.log"))


    return outfile



def HLA2EXON234(markerchoicefile, inbgl, outbgl, inmarker, outmarker):

    """
    Originally, this function was 'HLA2EXON234.py' source file.

    """


    selectMarkers = {}

    # did hla-a from here https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000206503;r=6:29941260-29945884;t=ENST00000376809
    #also did hla-b
    HLA_EXON2_POSITION = [['HLA_A', str(floor((29942757+29943026)/2))], ['HLA_B', str(floor((31356957+31356688)/2))], ['HLA_C', str(floor((31271868+31271599)/2))], 
                          ['HLA_E', str(floor((30489726+30489995)/2))], ['HLA_F', str(floor((29723658+29723927)/2))], ['HLA_G', str(floor((29828047+29828316)/2))],
                          ['HLA_DMA', str(floor((32950803+32950519)/2))], ['HLA_DMB', str(floor((32938965+32938684)/2))], ['HLA_DOA', str(floor((33008261+33008013)/2))], ['HLA_DOB', str(floor((32815313+32815044)/2))],
                          ['HLA_DPA1', str(floor((33069886+33069641)/2))], ['HLA_DPB1', str(floor((33080672+33080935)/2))], ['HLA_DQA1', str(floor((32641310+32641558)/2))], ['HLA_DQB1', str(floor((32665067+32664798)/2))],
                          ['HLA_DRA', str(floor((32442448+32442693)/2))],
                          ['HLA_DRB1', str(floor((32584378+32584109)/2))], ['HLA_DRB3', str(floor((32455037+32454768)/2))], ['HLA_DRB4', str(floor((32547835+32547566)/2))], ['HLA_DRB5', str(floor((32522174+32521905)/2))]]
    
    HLA_EXON3_POSITION = [['HLA_A', str(floor((29943268+29943543)/2))], ['HLA_B', str(floor((31356442+31356167)/2))], ['HLA_C', str(floor((31271348+31271073)/2))], 
                          ['HLA_E', str(floor((30490240+30490515)/2))], ['HLA_F', str(floor((29724173+29724448)/2))], ['HLA_G', str(floor((29828543+29828818)/2))],
                          ['HLA_DMA', str(floor((32949889+32949611)/2))], ['HLA_DMB', str(floor((32937456+32937172)/2))], ['HLA_DOA', str(floor((33007592+33007311)/2))], ['HLA_DOB', str(floor((32814601+32814320)/2))],
                          ['HLA_DPA1', str(floor((33069300+33069019)/2))], ['HLA_DPB1', str(floor((33084950+33085231)/2))], ['HLA_DQA1', str(floor((32641972+32642253)/2))], ['HLA_DQB1', str(floor((32662248+32661967)/2))],
                          ['HLA_DRA', str(floor((32443185+32443466)/2))],
                          ['HLA_DRB1', str(floor((32581838+32581557)/2))], ['HLA_DRB3', str(floor((32452465+32452184)/2))], ['HLA_DRB4', str(floor((32544828+32544547)/2))], ['HLA_DRB5', str(floor((32519651+32519370)/2))]]
    
    HLA_EXON4_POSITION = [['HLA_A', str(floor((29944122+29944397)/2))], ['HLA_B', str(floor((31355592+31355317)/2))], ['HLA_C', str(floor((31270485+31270210)/2))],
                          ['HLA_E', str(floor((30491137+30491412)/2))], ['HLA_F', str(floor((29725031+29725306)/2))], ['HLA_G', str(floor((29829418+29829693)/2))]]

    # ["A", "B", "C", "E", "F", "G", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5"]

    with open(markerchoicefile) as fin:
        for l in fin:
            c = l.split()
            selectMarkers[c[0]] = True

    with open(inbgl) as pf, open(outbgl, 'w') as of:
        for l in pf:
            c = l.split()
            header = c[:2]
            data = c[2:]
            if header[0] == "M" and header[1] not in selectMarkers:
                # 이 블럭에서 2-digit HLA allele들이 걸려나가줌
                continue

            if "HLA" in l:

                # Exon 2
                header[1] = header[1] + '_exon2'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')

                # Exon 3
                header[1] = (header[1].replace("_exon2", ""))
                header[1] = header[1] + '_exon3'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')

                # Exon 4
                header[1] = (header[1].replace("_exon3", ""))

                if "DRB1" in l:
                    continue
                if "DRB3" in l:
                    continue
                if "DRB4" in l:
                    continue
                if "DRB5" in l:
                    continue
                if "DMA" in l:
                    continue
                if "DMB" in l:
                    continue
                if "DRA" in l:
                    continue
                if "DOA" in l:
                    continue
                if "DOB" in l:
                    continue
                if "DQA1" in l:
                    continue
                if "DQB1" in l:
                    continue
                if "DPA1" in l:
                    continue
                if "DPB1" in l:
                    continue

                header[1] = header[1] + '_exon4'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')
                header[1] = (header[1].replace("_exon4", ""))

            if "HLA" in l:
                continue

            # "HLA_A_0101" 이런애들은 바로 위 if구문에서 마무리가 됨.
            # 여기는 사실상 "M rs1234" 이런 rs_id를 가지는 SNP들이 여기서 처리됨.
            newdata = []
            for j in range(len(data)):
                newdata.append(data[j])
            of.write(' '.join(header + newdata) + '\n')



    with open(inmarker) as mf, open(outmarker, 'w') as of:
        for l in mf:
            c = l.split()
            rsid_bp = c[:2]
            allele = c[2:]

            if rsid_bp[0] not in selectMarkers:
                continue

            if "HLA" in l:
                rsid_bp[0] = rsid_bp[0] + '_exon2'
                for i in range(19):
                    if HLA_EXON2_POSITION[i][0] in rsid_bp[0]:
                        rsid_bp[1] = HLA_EXON2_POSITION[i][1]
                        #of.write(' '.join((map(str, rsid_bp), map(str, allele))) + '\n')
                        of.write(' '.join(rsid_bp + allele) + '\n')
                        rsid_bp[0] = (rsid_bp[0].replace("_exon2", ""))
                    if HLA_EXON2_POSITION[i][0] in rsid_bp[0]:
                        continue

                rsid_bp[0] = rsid_bp[0] + '_exon3'
                for i in range(19):
                    if HLA_EXON3_POSITION[i][0] in rsid_bp[0]:
                        rsid_bp[1] = HLA_EXON3_POSITION[i][1]
                        #of.write(' '.join((map(str, rsid_bp), map(str, allele))) + '\n')
                        of.write(' '.join(rsid_bp + allele) + '\n')
                        rsid_bp[0] = (rsid_bp[0].replace("_exon3", ""))
                    if HLA_EXON3_POSITION[i][0] in rsid_bp[0]:
                        continue

                rsid_bp[0] = rsid_bp[0] + '_exon4'

                for i in range(6):
                    if HLA_EXON4_POSITION[i][0] in rsid_bp[0]:
                        rsid_bp[1] = HLA_EXON4_POSITION[i][1]

                        #of.write(' '.join((map(str, rsid_bp), map(str, allele))) + '\n')
                        of.write(' '.join(rsid_bp + allele) + '\n')
                        rsid_bp[0] = (rsid_bp[0].replace("_exon4", ""))
                    if HLA_EXON4_POSITION[i][0] in rsid_bp[0]:
                        continue

            if "HLA" in l:
                continue # 뭐던동 marker의 label이 'HLA_A_0101' 이런식이면 여기서 마무리가 됨.

            of.write(' '.join(rsid_bp + allele) + '\n') # 여기도 얘가 rs_id가지는 marker들을 담당하는 부분.


    return [outbgl, outmarker]



# def MakeExon234(__exonN__, markerchoicefile, inbgl, outbgl, inmarker, outmarker):
#
#
#     """
#
#     - Created by wschoi.
#     - Enhacned version of 'HLA2EXON234.py'
#
#     """
#
#
#
#     HLA_POSITION = {
#         'exon2': {'A': '30018647', 'C': '31347489', 'B': '31432578', 'DRB1': '32659998', 'DQA1': '32717189', 'DQB1': '32740687', 'DPA1': '33145518', 'DPB1': '33156558'},
#         'exon3': {'A': '30019161', 'C': '31346966', 'B': '31432060', 'DRB1': '32657452', 'DQA1': '32717867', 'DQB1': '32737862', 'DPA1': '33144914', 'DPB1': '33160845'},
#         'exon4': {'A': '30020015', 'C': '31346103', 'B': '31431210'}
#     }
#
#
#
#     ### Obtaining marker labels from sorted *.markers file(`sort_markers`)
#
#     selectMarkers = []
#     with open(markerchoicefile) as fin:
#         selectMarkers = [p.match(string=l).group(1) for l in fin]
#
#     # print(selectMarkers)
#
#
#     ## Changing rows for exon N in beagle file.
#
#     with open(inbgl) as pf, open(outbgl, 'w') as of:
#         for l in pf:
#
#             m = p.match(l)
#             header = [m.group(1), m.group(2)]
#
#             if header[0] == "M" and header[1] not in selectMarkers:
#                 # Excluding out 2-digit HLA alleles or Bizarre markers.
#                 continue
#
#
#             m2 = p_HLA_2field.match(header[1])
#
#             if m2:
#
#                 hla = m2.group(1)
#
#                 if __exonN__ == 'exon4' and not isClassI[hla]:
#                     # new_line = l
#                     continue
#                 else:
#                     # HLA alleles (with 4-digit(2-field))
#                     header2 = ' '.join([header[0], header[1]+'_{}'.format(__exonN__)])
#                     new_line = p.sub(repl=header2, string=l)
#                     # body = p.sub(repl='', string=l)
#                     # new_line = header2 + body
#
#                 of.write(new_line)
#
#             else:
#                 # Normal SNPs (ex. rs969931)
#                 of.write(l)
#
#
#     ### Changin rows for exon N in markers file(*.markers)
#     with open(inmarker) as mf, open(outmarker, 'w') as of:
#         for l in mf:
#
#             m = p.match(l)
#             rsid_bp = [m.group(1), m.group(2)]
#             # print(rsid_bp)
#
#             if rsid_bp[0] not in selectMarkers:
#                 continue
#
#
#             m2 = p_HLA_2field.match(rsid_bp[0])
#
#             if m2:
#                 hla = m2.group(1)
#
#                 if __exonN__ == 'exon4' and not isClassI[hla]:
#                     continue
#                 else:
#                     new_marker_name = rsid_bp[0]+'_{}'.format(__exonN__)
#                     new_bp = HLA_POSITION[__exonN__][hla]
#
#                     new_line = p.sub(repl=' '.join([new_marker_name, new_bp]), string=l)
#                     # print("Before : {}".format(l))
#                     # print("After : {}".format(new_line))
#
#                     of.write(new_line)
#             else:
#                 of.write(l)
#
#
#     return [outbgl, outmarker]


if __name__ == '__main__':

    """
    < Make_EXON234_Panel.py >
    
    INPUT : (0) 'exon2', 'exon3' or 'exon4', (1) Prefix of Reference Panel, (2) Prefix of Output file, 
            (3) Path to 'beagle2linkage.jar' file, (4) Path to Plink(v1.07) file.
            
    OUTPUT : Three copies of the reference panel of which the HLA markers have genomic position of middle point of Exon 2,3,4.
    
    """

    [__exonN__, infile, outfile, BEAGLE2LINKAGE, PLINK] = sys.argv[1:6]

    # ### Temporary Hardcoding
    # infile = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/HLA_PANEL/T1DGC/T1DGC_REF'
    # outfile = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190521_EXON234/T1DGC_REF_exon4'
    #
    # p_beagle2linkage = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/beagle2linkage.jar'
    # BEAGLE2LINKAGE = 'java -jar '+p_beagle2linkage
    #
    # p_plink = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/plink107/osx/plink'
    # PLINK = ' '.join([p_plink, '--noweb --allow-no-sex'])

    Make_EXON234_Panel(__exonN__, infile, outfile, BEAGLE2LINKAGE, PLINK)