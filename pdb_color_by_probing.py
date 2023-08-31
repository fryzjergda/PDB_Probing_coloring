#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import os


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

    return in_pdb, out_pdb, base_react, back_react


def check_color():


    backbone  = ["P  ", "OP1", "OP2", "O3'", "O5'", "O3*", "O5*", "C3'", 
                "C4'", "C5'", "C3*", "C4*", "C5*"]
    base = ["N9 ", "C8 ", "N7 ", "C5 ", "C4 ", "N3 ", "O2 ", "N2 ", "N1 ",
                 "C6 ", "O6 ", "C2 ", "N6 ", "O4 ", "N4 "]
    sugar = ["C3'", "C4'", "C5'", "C3*", "C4*", "C5*", "C1'", "C1*", "C2'",
                 "C2*", "O2'", "O2*", "O4'","O4*"]

    
    
def read_pdb():
    pdb_list = []
    with open(in_pdb, "rb") as f:
        for line in f:
            line = line.decode("UTF-8")
            if (line[:6] == "ATOM  " and line[16] != "B" and line[16] != "G") or \
                (line[:6] == "HETATM" and line[17:20] == "GTP" and \
                line[15] != "B" and line[15] != "G" and line[14] != "B" and \
                line[14] != "G" and line[16] != "B" and line[16] != "G"):
                pdb_list.append(line)
            elif line[:6] == "TER   ":
                break
    
    return pdb_list



def change_bfactor():
    
    back  = ["P  ", "OP1", "OP2", "O3'", "O5'", "O3*", "O5*", "C3'", 
                "C4'", "C5'", "C3*", "C4*", "C5*"]
    base = ["N9 ", "C8 ", "N7 ", "C5 ", "C4 ", "N3 ", "O2 ", "N2 ", "N1 ",
                 "C6 ", "O6 ", "C2 ", "N6 ", "O4 ", "N4 "]
    sugar = ["C3'", "C4'", "C5'", "C3*", "C4*", "C5*", "C1'", "C1*", "C2'",
                 "C2*", "O2'", "O2*", "O4'","O4*"]

    backbone = back + sugar
    
    pdb_out_list = pdb_in_list.copy()
    
    
    for i in range(0, len(reactivities_back_df)):
        #print(reactivities_base_df.iloc[i], "df") 
        for k in range(0, len(pdb_in_list)):
            if (int(pdb_in_list[k][22:26]) == i+1):
                line = list(pdb_out_list[k])
                if pdb_in_list[k][12:16].replace(" ","") not in base:
            
                    re = list(str("%.2f" % reactivities_back_df.iloc[i]["react"]))
                    put = [" ", " "]+re
                    if re == ['-', '9', '9', '9', '.', '0', '0']:
                        re = [" ", "-", "1", ".", "0", "0"]
                        put = re
                    line[60:66] = put
                    line = ''.join(line)
                elif pdb_in_list[k][12:16].replace(" ","") in backbone:
                    re = list(str("%.2f" % reactivities_back_df.iloc[i]["react"]))
                    put = [" ", " "]+re
                    if re == ['-', '9', '9', '9', '.', '0', '0']:
                        re = [" ", "-", "1", ".", "0", "0"]
                        put = re
                        
                    line[60:66] = put
                    line = ''.join(line)
                pdb_out_list[k] = line


    out_list = ""

    for i in range(0, len(pdb_out_list)):
        out_list += ''.join(pdb_out_list[i])
        
        
    write_output(out_list)
        
    
def write_output(pdb_out_list):

    if out_pdb == '':
        #outfile_name = reactivity.replace('.react','') + '_' + color_by + '.pdb'
        outfile_name = in_pdb.replace('.pdb','') + '_' + 'bfactor_ligand_binding_score' + '.pdb'
    else:
        outfile_name = out_pdb

    out= open(outfile_name, "w")

    out.write(pdb_out_list)
    out.close()
                    

def read_react(reactivity):

    with open(reactivity, "rb") as f:  # check number of lines in reactivity input file
        count = sum(1 for _ in f)

    if count >2:  # assume its reactivtity format
        react_df = pd.read_csv(reactivity, sep='\s+',index_col=False,header=None, encoding = "ISO-8859-1", 
                    names=["num","react"])
        react_df.set_index("num", inplace=True)
    else:  # assum its oneline format - copied from sheet
        react_df = pd.read_csv(reactivity, sep='\t',index_col=False,header=None, encoding = "ISO-8859-1")
        react_df = react_df.transpose()[0]
        react_df = pd.DataFrame(react_df)
        react_df.index.name="num"
        react_df.index += 1
        react_df.columns = ["react"]

    return react_df



if __name__ == '__main__':

        
    #try:
    #    os.system("dos2unix ./*")
    #except:
    #    pass
                     
    in_pdb, out_pdb, base_in, back_in = argument_parser()
    
    pdb_in_list = read_pdb()
    
    if base_in:
        reactivities_base_df = read_react(base_in)
    else:
        reactivities_base_df = ""
    
    if back_in:    
        reactivities_back_df = read_react(back_in)
        print(reactivities_back_df)
    else:
        print("else_back")
    
    
    change_bfactor()
    
    