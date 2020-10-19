#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse


def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, 
                        prog='pdb_color_by_probing.py', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input_pdb", dest="pdb_in",   
                        help="Input PDB file.")
    parser.add_argument("-b", "--base_reactivtity_file", required=False, dest="base_in",
                        help="File with reactivities for bases. [DMS, CMCT]")
    parser.add_argument("-bb", "--backbone_reactivtity_file", required=False, dest="back_in",
                        help="File with reactivities for backbone. [shape]")
    parser.add_argument("-o", "--output_pdb", required=False, dest="pdb_out", default='',
                        help="Name of output PDB file.")


    args = parser.parse_args()

    in_pdb = args.pdb_in
    out_pdb = args.pdb_out
    base_react = args.base_in
    back_react = args.back_in

    return in_pdb, out_pdb, base_in, back_in
