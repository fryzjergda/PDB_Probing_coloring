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
    parser.add_argument("-r", "--reactivtity_file", required=True, dest="react_in",
                        help="File with reactivity. [copied from sheet, or ractivity format]")
    parser.add_argument("-o", "--output_pdb", required=False, dest="pdb_out", default='',
                        help="Name of output PDB file.")
    parser.add_argument("-c", "--color_by", required=False, dest="color", default='all', 
                        choices=['all', 'base', 'back', 'sugar','basesug', 'backsug'],
                        help="Which part of RNA will have raectivities set as b factor.")
    parser.add_argument("-s", "--sequence", required=False, dest="seq_in", default='',
                        help="Name of output PDB file.")


    args = parser.parse_args()

    in_pdb = args.pdb_in
    out_pdb = args.pdb_out
    reactivity = args.react_in
    color_by = args.color
    sequence_in = args.seq_in

    return in_pdb, out_pdb, reactivity, color_by, sequence_in


def change_bfactor():
    
#    print(reactivities_df)
    pdb_out_list = pdb_in_list.copy()
    
    first_res_pdb = int(pdb_in_list[0][22:26])
    print(first_res_pdb+pdb_start_index)
    last_res_pdb = int(pdb_in_list[-1][22:26])
    print(last_res_pdb, 'last')

    pdb_index_to_color=first_res_pdb+pdb_start_index
    pdb_stop_color = last_res_pdb-pdb_end_index
    print(pdb_stop_color, 'stop')
    #quit()
    pdb_range = list(range(pdb_stop_color, pdb_len))
    print(pdb_range)
    df_index_to_color=seq_start_index 
    print(pdb_range)
    #quit()
    print(df_index_to_color, reactivities_df.shape[0], seq_end_index)
#    quit()
    for i in range(df_index_to_color, reactivities_df.shape[0]-seq_end_index):  # go through the reactivity df
        for k in range(0, len(pdb_in_list)):  # go through the whole pdb
            if (int(pdb_in_list[k][22:26]) == pdb_index_to_color) and (int(pdb_in_list[k][22:26]) not in pdb_range):  # check residue number and get pdb rows of sthis residue 
                line = list(pdb_out_list[k])
                if pdb_in_list[k][13:16] in color:  # check if the atom is in color list, and change bfactor if yes 
                    re = list(str("%.2f" % reactivities_df.iloc[i]["react"]))
                    print("kupka", re)
                    
                    put = [" ", " "]+re
                    line[60:66] = put
                    line = ''.join(line)
                else:  # color other atoms to grey (b factor 1.25)
                    put = [" ", " ", "1", ".", "2", "5"]
                    line[60:66] = put
                    line = ''.join(line)
                pdb_out_list[k] = line
            elif (int(pdb_in_list[k][22:26]) in pdb_range) or (int(pdb_in_list[k][22:26]) <(first_res_pdb+pdb_start_index)):
                line = list(pdb_out_list[k])
                put = [" ", "-", "1", ".", "0", "0"]
                line[60:66] = put
                line = ''.join(line)
                pdb_out_list[k] = line

        pdb_index_to_color+=1
    
 #   for i in range(0, len(pdb_in_list)):
        
    
    
    
#    print(pdb_out_list)
    write_output(pdb_out_list)


def write_output(pdb_out_list):
    
    if out_pdb == '':
        outfile_name = reactivity.replace('.react','') + '_' + color_by + '.pdb'
    #    outfile_name = in_pdb.replace('.pdb','') + '_' + color_by + '.pdb'
    else:
        outfile_name = out_pdb    

    out= open(outfile_name, "w")
    for i in pdb_out_list:
        out.write(i)
    out.close()


def check_color():
    
    
    backbone  = ["P  ", "OP1", "OP2", "O3'", "O5'", "O3*", "O5*", "C3'", 
                "C4'", "C5'", "C3*", "C4*", "C5*"]
    base = ["N9 ", "C8 ", "N7 ", "C5 ", "C4 ", "N3 ", "O2 ", "N2 ", "N1 ",
                 "C6 ", "O6 ", "C2 ", "N6 ", "O4 ", "N4 "]
    sugar = ["C3'", "C4'", "C5'", "C3*", "C4*", "C5*", "C1'", "C1*", "C2'",
                 "C2*", "O2'", "O2*", "O4'","O4*"]

    if color_by == 'all':
        color = backbone + base + sugar
    elif color_by == 'base':
        color = base
    elif color_by == 'back':
        color = backbone
    elif color_by == 'sugar':
        color = sugar
    elif color_by == 'basesug':
        color = base + sugar
    elif color_by == 'backsug':
        color = backbone + sugar

    return color


def read_pdb():
    pdb_list = []
    with open(in_pdb) as f:
        for line in f:
            if (line[:6] == "ATOM  " and line[16] != "B" and line[16] != "G") or \
                (line[:6] == "HETATM" and line[17:20] == "GTP" and \
                line[15] != "B" and line[15] != "G" and line[14] != "B" and \
                line[14] != "G" and line[16] != "B" and line[16] != "G"):
                pdb_list.append(line)
            elif line[:6] == "TER   ":
                break    
    return pdb_list    

def read_sequence():
    
    with open(sequence_in, 'r') as f:  # read sequence from file
        seq=''
        for line in f:
          if line[0] !=">":
            #print f.readline()
            seq += line.strip().replace("T","U").replace("t","u").upper()  # sequence with SH
    
    
    userPDB = in_pdb
    s = RNAStructure(userPDB)
    pdbSeq = s.get_seq(compact=True).strip()
    #s.get_rnapuzzle_ready()
    #print(s)
    
    print(pdbSeq, 'tu')    
    #print(seq)
    lcs = ','.join(longest_common_substring(seq, pdbSeq))
    print(lcs)
    pdb_ind = pdbSeq.find(lcs)
    #print(pdb_ind, 'pdb_ind')
    seq_ind = seq.find(lcs)
    print(seq_ind,"seq_ind")
    #print(seq_ind, 'seq_ind')
    #print(lcs, 'lcs')
    seq_end = seq[::-1].find(lcs[::-1])
    print(seq[::-1])
    print(lcs[::-1])
#    if len(seq[::-1]) == len(lcs[::-1]):
#            seq_end = len(seq)
        
    print(seq_end,"seq_end")
    #print(seq[::-1])
    #print(lcs[::-1])
    #print(seq_end, 'seq_end')
    pdb_rna = ''.join([c for c in pdbSeq if c.isupper()])
    print(pdb_rna)
    pdb_end = pdb_rna[::-1].find(lcs[::-1])-1
#    if len(pdb_rna[::-1]) == len(lcs[::-1]):
#        pdb_end = len(seq)
    print(pdb_end)
    #print(pdb_rna[::-1])
    #print(lcs[::-1])
    #print(pdb_end, 'pdb_end')
    pdb_len = len(pdbSeq)
    #print(','.join(substr))
    print(pdb_ind, seq_ind, seq_end, pdb_end, pdb_len)
    return pdb_ind, seq_ind, seq_end, pdb_end, pdb_len

def longest_common_substring(S, T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

    return lcs_set
    
    
def read_react():
    
    with open(reactivity) as f:  # check number of lines in reactivity input file
        count = sum(1 for _ in f)
    
    if count >2:  # assume its reactivtity format
        react_df = pd.read_csv(reactivity, sep=' ',index_col=False,header=None, 
                    names=["num","react"])
        react_df.set_index("num", inplace=True)
    else:  # assum its oneline format - copied from sheet
        react_df = pd.read_csv(reactivity, sep='\t',index_col=False,header=None)
        react_df = react_df.transpose()[0]
        react_df = pd.DataFrame(react_df)
        react_df.index.name="num"
        react_df.index += 1
        react_df.columns = ["react"]

    return react_df


if __name__ == '__main__':

    in_pdb, out_pdb, reactivity, color_by, sequence_in = argument_parser()

    pdb_start_index, seq_start_index, seq_end_index, pdb_end_index, pdb_len = read_sequence()
#    quit()
    pdb_in_list = read_pdb()
    reactivities_df = read_react()
    color = check_color()
        
    change_bfactor()    
    
    print('\nB factor changed succesfully!\nGood job!\n')


