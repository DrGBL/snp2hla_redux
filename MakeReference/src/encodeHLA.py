# -*- coding: utf-8 -*-

# (2017/11/27) recoded by Wanson Choi
import os, re
import argparse, textwrap
from math import floor


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "E", "F", "G", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5"]

# (2022.) WILL NOT USE HG19 AT ALL!!!!!!!!!!
# here I use the average position of exons 2,3 as the position for each gene (the mane/ensembl canonical transcripts if multiple choices of transcripts available)
genepos_hg = {"38": {"A": floor((29942757+29943543)/2), "B": floor((31356957+31356167)/2), "C": floor((31271868+31271073)/2), 
                     "E": floor((30489726+30490515)/2), "F": floor((29723658+	29724448)/2), "G": floor((29828047+29828818)/2),
                     "DMA": floor((32950803+32949611)/2), "DMB": floor((32938965+32937172)/2), "DOA": floor((33008261+33007311)/2), "DOB": floor((32815313+32814320)/2), 
                     "DPA1": floor((33073669+	33069641)/2), "DPB1": floor((33080672+33085231)/2), 
                     "DQA1": floor((32641310+32642253)/2), "DQB1": floor((32665067+32661967)/2),
                     "DRA": floor((32442448+32443466)/2),
                     "DRB1": floor((32584378+32581557)/2), "DRB3": floor((32455037+32452184)/2), "DRB4": floor((32547835+32544547)/2), "DRB5": floor((32522174+32519370)/2),
                     "DPA2": floor((33093137+33092119))/2), "H": floor((29887955+29888742))/2), "J": floor((30006800+30007593))/2), "K": floor((29926660+29927452))/2), "L": floor((30260397+30261479))/2),
                     "MICA": floor((31410543+31411359))/2), "MICB": floor((31505617+31506430))/2), "V": floor((29792437+29793136))/2)}}


def encodeHLA(_CHPED, _OUTPUT, __asSmallLetter=True, __addDummyMarker=False, __previous_version=False):



    ### Intermediate path.
    _OUTPUT = _OUTPUT if not _OUTPUT.endswith('/') else _OUTPUT.rstrip('/')
    if bool(os.path.dirname(_OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(_OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        INTERMEDIATE_PATH = "./"



    ### Acquiring `HLA_allele_sets`.

    HLA_allele_sets = {HLA_names[i]: [] for i in range(0, len(HLA_names))}

    p_1field = re.compile(r'\w+\*\d{2,3}')

    with open(_CHPED, 'r') as f_chped:

        count = 0

        for l in f_chped:

            """
            l[:6] := ("FID", "IID", "PID", "MID", "Sex", "Phe")
            l[6:8] := HLA-A
            l[8:10] := HLA-B
            ...
            l[20:22] := HLA-DRB1
            """

            t_line = re.split(r'\s+', l.rstrip('\n'))
            # print(t_line)

            for i in range(0, len(HLA_names)):

                idx1 = 2*i + 6
                idx2 = idx1 + 1

                al1 = t_line[idx1]
                al2 = t_line[idx2]


                # Allele 1
                if al1 != "0":

                    if al1 not in HLA_allele_sets[HLA_names[i]]:
                        HLA_allele_sets[HLA_names[i]].append(al1)

                    m = p_1field.match(al1)

                    if m:

                        al1_1field = m.group()

                        if al1_1field not in HLA_allele_sets[HLA_names[i]]:
                            HLA_allele_sets[HLA_names[i]].append(al1_1field)


                # Allele 2
                if al2 != "0":

                    if al2 not in HLA_allele_sets[HLA_names[i]]:
                        HLA_allele_sets[HLA_names[i]].append(al2)

                    m = p_1field.match(al2)

                    if m:

                        al2_1field = m.group()

                        if al2_1field not in HLA_allele_sets[HLA_names[i]]:
                            HLA_allele_sets[HLA_names[i]].append(al2_1field)


            count += 1
            # if count > 5 : break


    for i in range(0, len(HLA_names)):
        HLA_allele_sets[HLA_names[i]].sort()


    # # Result checking
    # print("\nHLA alleles.")
    # for k, v in HLA_allele_sets.items():
    #     print("{}: {}".format(k, v))



    ### Making a new *.HLA.map file.

    map_LABELS = ['_'.join(["HLA", HLA_allele_sets[HLA_names[i]][j]]) for i in range(0, len(HLA_names)) for j in range(0, len(HLA_allele_sets[HLA_names[i]]))]
    # print(map_LABELS)

    map_POS = [str(genepos_hg["38"][HLA_names[i]]) for i in range(0, len(HLA_names)) for z in range(0, len(HLA_allele_sets[HLA_names[i]]))]
    # print(map_POS)

    with open(_OUTPUT + ".map", 'w') as f_HLA_map:
        f_HLA_map.writelines(('\t'.join(["6", map_LABELS[i], "0", map_POS[i]]) + "\n" for i in range(0, len(map_LABELS))))

        if __addDummyMarker:
            f_HLA_map.write('\t'.join(["6", "dummy_marker", "0", "33999999"]) + "\n")



    ### Making a new *.HLA.ped file.

    with open(_OUTPUT + ".ped", 'w') as f_HLA_ped:
        f_HLA_ped.writelines(
            MakeHLAPed(_CHPED, HLA_allele_sets, __asSmallLetter=__asSmallLetter, __addDummyMarker=__addDummyMarker,
                       __previous_version=__previous_version))



    return [_OUTPUT+".ped", _OUTPUT+".map"]





# (2019. 1. 3.) Introduced for memory issues.
def PrintGenotypes4(_allele1, _allele2, _HLA_allele_sets_byHLA, __asSmallLetter=False, __previous_version=False):

    l_output = []

    _present_ = "p" if __asSmallLetter else "P"
    _absent_ = "a" if __asSmallLetter else "A"



    for i in range(0, len(_HLA_allele_sets_byHLA)):

        _ALLELE = _HLA_allele_sets_byHLA[i]

        G1 = "-1"
        G2 = "-1"

        if _allele1 != "0" and _allele2 != "0":

            if __previous_version:

                p_4digit = re.compile(r'\w+:(\d{2}):(\d{2})[A-Z]?$')
                p_2digit = re.compile(r'\w+:(\d{2})[A-Z]?$')


                # Allele 1
                if p_4digit.match(_allele1):
                    _allele1 = ''.join(p_4digit.findall(_allele1).pop())
                elif p_2digit.match(_allele1):
                    _allele1 = p_2digit.findall(_allele1).pop()

                if _allele1 == _ALLELE or _allele1[:2] == _ALLELE:
                    G1 = _present_
                elif len(_allele1) == 2 and len(_ALLELE) == 4 and _ALLELE[:2] == _allele1:
                    G1 = "0"
                else:
                    G1 = _absent_


                # Allele 2
                if p_4digit.match(_allele2):
                    _allele2 = ''.join(p_4digit.findall(_allele2).pop())
                elif p_2digit.match(_allele2):
                    _allele2 = p_2digit.findall(_allele2).pop()

                if _allele2 == _ALLELE or _allele2[:2] == _ALLELE:
                    G2 = _present_
                elif len(_allele2) == 2 and len(_ALLELE) == 4 and _ALLELE[:2] == _allele2:
                    G2 = "0"
                else:
                    G2 = _absent_


            else:
                # Dealing with generalized 4-field HLA alleles.

                p_1field = re.compile(r'\w+\*\d{2,3}')


                # Allele 1
                m = p_1field.match(_allele1)
                _al1_1field = m.group() # Just assume that given alleles is definitely in the form of r'\w+\*(\d{2,3}:?)+'

                # Determining 'Present' or 'Absent'.
                if _allele1 == _ALLELE or _al1_1field == _ALLELE:
                    G1 = _present_
                else:
                    G1 = _absent_


                # Allele 2
                m = p_1field.match(_allele2)
                _al2_1field = m.group()
                if _allele2 == _ALLELE or _al2_1field == _ALLELE:
                    G2 = _present_
                else:
                    G2 = _absent_


        else:
            # If at least one HLA allele which is given in *.chped is "0", then consider both of them are "0"
            G1 = "0"
            G2 = "0"




        if G1 == "0" or G2 == "0":
            l_output.append("0")
            l_output.append("0")
        else:
            l_output.append(G1)
            l_output.append(G2)



    return '\t'.join(l_output)




def MakeHLAPed(_CHPED, _HLA_allele_sets, __asSmallLetter=False, __addDummyMarker=False, __previous_version=False):

    with open(_CHPED, 'r') as f_chped:

        count = 0

        for l in f_chped:

            t_line = re.split(r'\s+', l.rstrip('\n'))

            """
            t_line[:6] := ("FID", "IID", "PID", "MID", "Sex", "Phe")
            t_line[6:8] := HLA-A
            t_line[8:10] := HLA-B
            ...
            t_line[20:22] := HLA-DRB1
            """

            t_iterator = range(0, len(HLA_names)) if not __previous_version else HLA_names

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                PrintGenotypes4(t_line[2 * i + 6], t_line[2 * i + 7], _HLA_allele_sets[HLA_names[i]],
                                __asSmallLetter=__asSmallLetter, __previous_version=__previous_version)
                for i in t_iterator if len(_HLA_allele_sets[HLA_names[i]]) > 0
            ])

            __return__ = '\t'.join([__ped_info__, __genomic_info__])


            if __addDummyMarker:
                dummy_markers = '\t'.join(['d', 'D'] if bool(count % 2) else ['D', 'd'])
                __return__ = '\t'.join([__return__, dummy_markers])


            yield __return__ + "\n"

            count += 1





if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        encodeHLA.py

        This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
        The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
                                            2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

    ###########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-chped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work like original version.\n\n",
                        action='store_true')
    parser.add_argument("--asSmallLetter", help="\n'P'resent and 'A'bsent to 'p'resent and 'a'bsent.\n\n",
                        action='store_true')
    parser.add_argument("--addDummyMarker", help="\nAdd dummy marker to prevent the glitch in work with plink(1.07).\n\n",
                        action='store_true')



    ##### <for Test> #####

    # (2019. 01. 06.)
    # # --previous-version
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-hg", "18",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190107/_3_encodeHLA/_Case_HAPMAP_CEU.HLA",
    #                           "--previous-version"])

    # # Generalized 4-field HLA alleles
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "18",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190107/_3_encodeHLA_4field/_Case_HAPMAP_CEU.HLA_4fieldTest",
    #                           "--asSmallLetter",
    #                           "--addDummyMarker"])





    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    encodeHLA(args.chped, args.o, __asSmallLetter=(not args.asSmallLetter),
              __addDummyMarker=args.addDummyMarker, __previous_version=args.previous_version)
