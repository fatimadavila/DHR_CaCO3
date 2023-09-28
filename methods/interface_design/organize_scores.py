#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import glob
import re


def grep_pattern_file(pattern, file):
    match_list = []
    for line in open((file), 'r'):
        line_match = re.match(pattern, line)
        if line_match:
            match_list.append(line)
    return (match_list)


def get_pose_info(path_pdb):
    cwd = os.getcwd()
    pdb_name = os.path.basename(path_pdb)
    bb = pdb_name.split('_side')[0]
    side = pdb_name.split('side')[1].split('_')[1]
    des_aas = pdb_name.split('alloaas')[1].split('_')[1]
    aasrft = pdb_name.split('aasurftail')[1].split('_')[1]
    if pdb_name.split('_')[-1] == 'uncapped.pdb':
        n_id = pdb_name.split('_')[-2]
        cap = 'uncapped'
    else:
        n_id = pdb_name.split('_')[-1].strip('.pdb')
        cap = 'capped'
    path = os.path.join(cwd, os.path.dirname(path_pdb))
    return(pdb_name, bb, side, des_aas, aasrft, n_id, cap, path)


def write_to_file(filename, input_list):
    with open(filename, 'a+') as fout:
        line = '\t'.join(map(str, input_list))
        print(line)
        fout.write(line + '\n')
    return fout


def main(argv=None):

    # ARGUMENT-HANDLING BLOCK
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--glob_pattern', '-gptn', type=str,
                        help='Recursively look for pdb files into this input directory pattern.',
                        default='lattice8_7_params_*side*.pdb')
    parser.add_argument('--output_scorefile_name', '-o_scf_name',
                        type=str, help='Name of output scorefile. Specify/Include extension.',
                        default='scorefile.tsv')
    args = parser.parse_args()

    # GET INPUT PATHS AND PDBS
    path_pdb_list = glob.glob(os.path.join('**/', args.glob_pattern),
                              recursive=True)

    # PROCESS HEADER
    header = str(grep_pattern_file('label ', path_pdb_list[0]))
    header_list = header.strip("\\n']").split(' ')[1:]
    extra_terms_header = ['vanilla', 'side', 'des_aas', 'aasurftail',
                          'backbone', 'nstruct_id', 'capping', 'pdb_name', 'path']
    header_list = header_list + extra_terms_header
    write_to_file(args.output_scorefile_name, header_list)

    # GET AND CALCULATE SCORES FROM EACH ITERABLE
    for path_pdb_element in path_pdb_list:
        score_str = str(grep_pattern_file('pose ',
                                          path_pdb_element)[0])
        score_list = score_str.strip("\n").split(' ')[1:]
        vanilla_score = sum([float(i) for i in score_list][
                            :17]) + sum([float(i) for i in score_list][18:-1])
        pdb_name, bb, side, des_aas, aasrft, n_id, cap, path = get_pose_info(
            path_pdb_element)
        score_list = score_list + [vanilla_score, side,
                                   des_aas, aasrft, bb, n_id, cap, pdb_name, path]
        write_to_file(args.output_scorefile_name, score_list)


if __name__ == '__main__':
    sys.exit(main())
