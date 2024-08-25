#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip


parser = argparse.ArgumentParser(description='Create annotations from BED files for all chromosomes')
parser.add_argument('--bedname', action="store",
                    dest="bedname", type=str,
                    help='The name of the BED file to conver to .annot file.')
parser.add_argument('--bedfile_path', action="store",
                    dest="bedfile_path", type=str,
                    help='The path containing the BED files.')
parser.add_argument('--bimfile_path', action="store",
                    dest="bimfile_path", type=str,
                    help='The path containing the BIM files.')
parser.add_argument('--annot_path', action="store",
                    dest="annot_path", type=str,
                    help='The path where to save the annotation files.')

args = parser.parse_args()
bedname = args.bedname
bedfile_path = args.bedfile_path
bimfile_path = args.bimfile_path
annot_path = args.annot_path



def make_annot_files(bed_for_annot, bimfile, annot_file):
    print('making annot file')
    df_bim = pd.read_csv(bimfile,
            delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2 - 1, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim) # all snps
    annotbed = bimbed.intersect(bed_for_annot) # snps annotated with 1
    bp = [x.start + 1 for x in annotbed]
    df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
    df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
    df_annot.fillna(0, inplace=True)
    
    # # using thinannot
    # df_annot = df_annot[['ANNOT']].astype(int)

    # not using thinannot
    temp = df_annot[['ANNOT']].astype(int)
    df_annot = pd.concat([df_bim.iloc[:,[0,3,1,2]], temp], axis = 1)
    
    if  annot_file.endswith('.gz'):
        with gzip.open(annot_file, 'wb') as f:
            df_annot.to_csv(f, sep = "\t", index = False)
    else:
        df_annot.to_csv(annot_file, sep="\t", index=False)



for numchr in range(1, 23, 1):
    bimfile = bimfile_path + "/" + "1000G.EUR.QC." + str(numchr) + ".bim"
    annot_file = annot_path + "/" + bedname + "." + str(numchr) + ".annot.gz"
    bedfile = bedfile_path + "/" + bedname + ".bed"
    bed_for_annot = BedTool(bedfile).sort()
    make_annot_files(bed_for_annot, bimfile, annot_file)
    print("We are at chrom : " + str(numchr))


