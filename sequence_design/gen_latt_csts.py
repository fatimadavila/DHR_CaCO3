#! /usr/bin/env python3

from pyrosetta import *
import glob
init(extra_options='-mute=all')
import argparse
import sys
import re
import os


def get_sec_struct_dict(some_pose):
    some_pose.display_secstruct()
    # Split pose in dictionary containing helix number and span of residues
    # corresponding to it
    key_number = 1
    d = {}
    for i in range(some_pose.total_residue()):
        j = i + 1
        if j == 1:
            if some_pose.secstruct(j) is 'H':
                new_key_name = key_number
                l = []
                l.append(j)
                d[new_key_name] = l
                key_number = key_number + 1
        else:
            if some_pose.secstruct(
                    j) is 'H' and some_pose.secstruct(i) is 'H':
                l.append(j)
                d[new_key_name] = l
            elif some_pose.secstruct(j) is 'H' and some_pose.secstruct(i) != 'H':
                new_key_name = key_number
                l = []
                l.append(j)
                d[new_key_name] = l
                key_number = key_number + 1
    return d


def generate_lattice_csts(
        wd2, cst_file_name, lattice_space, cst_tolerance, all_hres_side, rep_len):
    cst_template = 'AtomPair CA %s CA %s HARMONIC %s %s'
    c = str('%.2f' % lattice_space).rjust(4)
    d = str('%.2f' % cst_tolerance).rjust(2)
    with open(wd2 + '/' + cst_file_name + '.cst', 'w') as fout:
        for res in all_hres_side:
            a = str(int(res)).rjust(2)
            b = str(int(rep_len + res)).rjust(3)
            line = cst_template % (a, b, c, d)
            print(line)
            fout.write(line + '\n')


def main(wd2, pdb, input_side, lattice_space, cst_tolerance):
    input_pose = pose_from_pdb(pdb)
    ss_dict = get_sec_struct_dict(input_pose)
    pos_len = input_pose.total_residue()
    rep_len = int(pos_len / (len(ss_dict) / 2))
    cst_file_name = (pdb.split('.pdb')[0]) + '_lattice_csts'
    if input_side == 'A':
        cst_file_name += '_A'
        all_hres_side = ss_dict[1]
        for k in range(len(ss_dict.keys()) - 2):
            if int(k) % 2 != 0:
                all_hres_side += ss_dict[k]
    elif input_side == 'B':
        cst_file_name += '_B'
        all_hres_side = []
        for k in range(len(ss_dict.keys()) - 2):
            if int(k) % 2 == 0:
                print(k)
                all_hres_side += ss_dict[k]
    generate_lattice_csts(
        wd2,
        cst_file_name,
        lattice_space,
        cst_tolerance,
        all_hres_side,
        rep_len)

if __name__ == "__main__":

    cwd = os.getcwd()

    ArgParser = argparse.ArgumentParser(
        description='Generate lattice constraints given pdb, cst params and side. DHR-specific.')
    ArgParser.add_argument(
        '-pdb',
        type=str,
        help='only one input pdb file',
        default=None)
    ArgParser.add_argument(
        '-lsp',
        type=str,
        help='Ideal lattice space distance to the 0.00 Anstrom',
        default=None)
    ArgParser.add_argument(
        '-cst_tol',
        type=str,
        help='Contraint tolerance to the nearest 0.00 Angstrom',
        default=None)
    ArgParser.add_argument("-s", "--side", type=str, choices=('A', 'B'),
                           help="One letter side identifier for the protein.E.g. 'A'",
                           default='A')
    Args = ArgParser.parse_args()

    main(
        wd2=cwd,
        pdb=Args.pdb,
        input_side=Args.side,
        lattice_space=float(
            Args.lsp),
        cst_tolerance=float(
            Args.cst_tol))
