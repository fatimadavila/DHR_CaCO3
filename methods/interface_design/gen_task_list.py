#!/usr/bin/env python3

import os
import glob

path_list = glob.glob('**/*.pdb')
cmd_f = 'cd {dir_name}; python3 ../mutate_EAD_surface.py -i {pdb_name} -s {side} -aas {aas} -aast {aast} -nstruct 30\n'.format
sides = ['A', 'B']
aas = ['EDA', 'ED', 'D']
aast = ['DNSTQKREDH', 'DNSTQEDH']
with open('task.list', 'w') as fout:
    for path_pdb_name in path_list:
        dire = os.path.dirname(path_pdb_name)
        pdb = os.path.basename(path_pdb_name)
        for s in sides:
            for aa in aas:
                for taa in aast:
                    print(cmd_f(dir_name=dire, pdb_name=pdb, side=s,
                                aas=aa, aast=taa))
                    line = cmd_f(dir_name=dire, pdb_name=pdb, side=s,
                                 aas=aa, aast=taa)
                    fout.write(line)
