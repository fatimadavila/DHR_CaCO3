#! /usr/bin/env python3

import os
import sys
import glob
import shutil
import argparse
import subprocess


def make_directory(directory):
    '''Make a directory if it doesn't exist'''
    if not os.path.isdir(directory):
        os.mkdir(directory)


def main():
    parser = argparse.ArgumentParser(
        description='prepare_for_protocol.py -lsp 8.7')
    parser.add_argument('-gptn', '--glob_pattern', type=str,
                        help='Input pattern to glob in this current \
                                                working directory',
                        default='*pdb')
    parser.add_argument('--lattice_spacing', '-lsp', type=str,
                        help='Ideal lattice space distance to the 0.00 \
                                                Angstrom',
                        default=None, required=True)
    parser.add_argument('--constraint_tolerance', '-cst_tol', type=str,
                        help='Contraint tolerance to the nearest 0.00 \
                                                Angstrom',
                        default=0.05)
    args = parser.parse_args()

    cwd = os.getcwd()
    pdb_list = glob.glob(os.path.join(cwd, args.glob_pattern))

    for pdb in pdb_list:
        os.chdir(cwd)
        pdb_name = os.path.basename(pdb)
        dirname = pdb_name.split('.pdb')[0]
        make_directory(dirname)
        shutil.move(os.path.join(cwd, pdb_name),
                    os.path.join(cwd, dirname))
        os.chdir(os.path.join(cwd, dirname))
        for side in ('A', 'B'):
            template_string = 'python3 /home/fadh/scripts/surface_tools/lattice_csts/gen_latt_csts.py -pdb {0} -s {1} -lsp {2} -cst_tol {3}'
            cmd_string = template_string.format(pdb_name, side,
                                                args.lattice_spacing,
                                                args.constraint_tolerance)
            print(cmd_string)
            subprocess.check_output(cmd_string.split())

if __name__ == "__main__":
    sys.exit(main())
