#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
from rna_tools_lib import RNAStructure
import sys

def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, 
                        prog='probing_bfactor_pdb.py', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--input_pdb", dest="pdb_in",   
                        help="Input PDB file.")
    parser.add_argument("-i2", "--input_pdb2", dest="pdb_in2",   
                        help="Input PDB file.")

    args = parser.parse_args()

    in_pdb = args.pdb_in
    in_pdb2 = args.pdb_in2
    
    return in_pdb, in_pdb2

def merge():

    with open(in_pdb) as f:
        pdb1 = f.readlines()

    #print(content)
    with open(in_pdb2) as f:
        pdb2 = f.readlines()

    pdb3 = ""
    for i in range(0,len(pdb1)):
        #print(pdb1[i])
        if pdb1[i][60:66] != "  1.25":
            pdb3 += pdb1[i]
         #   print("kupka")
        else:
            pdb3 += pdb2[i]

    print(pdb3)
if __name__ == '__main__':
    
    in_pdb, in_pdb2 = argument_parser()

    merge()