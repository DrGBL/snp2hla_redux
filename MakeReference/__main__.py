#-*- coding: utf-8 -*-
# python -m MakeReference

# import os, sys, re
import argparse, textwrap

from MakeReference.MakeReference_v2 import MakeReference_v2



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        MakeReference_v2.py

        

    #################################################################################################
                                     '''),
                                     add_help=False,
                                     prog='MakeReference')

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--variants", help="\nInput variants data file(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--chped", help="\nHLA Type Data(.chped)\n\n", required=True)
    parser.add_argument("--genes", help="\nComma separated list of HLA genes, in the same order as in the chped columns.\nIf not in the same order, there will be problems.\nDefault: A,B,C,E,F,G,H,J,K,L,V,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB4,DRB5,MICA,MICB.\n\n", default="A,B,C,E,F,G,H,J,K,L,V,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB4,DRB5,MICA,MICB")
    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg",
                        default="38")
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')
    parser.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')
    parser.add_argument("--mind", help="\nMaximum missing genotype rate per sample in reference. Default=0.3.\n\n", default=0.3)
    parser.add_argument("--hardy", help="\nMaximum acceptable Hardy-Weinberg equilibrium p-value for reference variants. Default=0.000001.\n\n", default=0.000001)
    parser.add_argument("--maf", help="\nMinimum allele frequency for reference variants. Default=0.01.\n\n", default=0.01)
    parser.add_argument("--miss", help="\nMaximum genotype missing rate for reference variants. Default=0.05.\n\n", default=0.05)
    parser.add_argument("--hla_maf", help="\nMinimum allele frequencies for HLA alleles to try to impute them. Default=0.0001.\n\n", default=0.0001)
    

    # Beagle 5.4
    parser.add_argument("--phasing", help="\nPerform phasing with Beagle 5.4.\n\n", action='store_true')
    parser.add_argument("--mem", help="\nJava Heap Memory size for Beagle 5.4. ex. 2g, 500m. Default: 4g.\n\n", default="4g")
    parser.add_argument("--tmp_folder", help="\nFolder where java temporary files are created. Default: /tmp.\n\n", default="/tmp")
    parser.add_argument("--nthreads", help="\nThe number of threads to use in Bealge 5.4. (ex. 2)\n\n", default=1, type=int)
    parser.add_argument("--burnin", help="\nNumber of burn-in MCMC iterations for Beagle 5.4. Default: 3.\n\n", default=3, type=int)
    parser.add_argument("--iter", help="\nNumber of MCMC iterations after burn-ins for Beagle 5.4. Default: 12.\n\n", default=12, type=int)
    parser.add_argument("--map", help="\nGenetic Map to use (optional).\n\n", default="null")
    parser.add_argument("--window", help="\nGenetic window to use. Needs to be at least 1.1 times the overlap. Default=2.8cM.\n\n", default=2.8)
    parser.add_argument("--overlap", help="Overlap between genetic windows. Default=1.25cM.\n\n", default=1.25)



    args = parser.parse_args()
    print(args)

    MakeReference_v2(_CHPED=args.chped, _OUT=args.out, _hg=args.hg, _genes=args.genes,
                     _variants=args.variants, _java_heap_mem=args.mem, _java_tmp_folder=args.tmp_folder,
                     _p_dependency=args.dependency, f_save_intermediates=args.save_intermediates,
                     f_phasing=args.phasing, _nthreads=args.nthreads,
                     _burnin=args.burnin, _iter=args.iter, _map=args.map,
                     _mind=args.mind, _hardy=args.hardy, _maf=args.maf, _miss=args.miss,
                     _hla_maf=args.hla_maf, _window=args.window, _overlap=args.overlap)