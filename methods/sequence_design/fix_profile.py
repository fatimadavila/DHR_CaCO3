import sys
import os
import math
import pyrosetta
pyrosetta.init()

def get_loop_residues(p):
    pre_loops = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector('L')
    loops = pyrosetta.rosetta.core.select.residue_selector.PrimarySequenceNeighborhoodSelector(1,1,pre_loops)
    residues = pyrosetta.rosetta.core.select.get_residues_from_subset(loops.apply(p)) 
    return residues

def get_surface_residues(p):
    surface_sasa = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    surface_sasa.set_use_sc_neighbors(False)
    surface_sasa.set_ball_radius(2.0)
    surface_sasa.set_cutoffs(1, 40)
    surface_sasa.set_layers(False, False, True)
    residues = pyrosetta.rosetta.core.select.get_residues_from_subset(surface_sasa.apply(p))
    return residues

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False



pose = pyrosetta.pose_from_file(sys.argv[1])

loop_residues = get_loop_residues(pose)
surface_residues = get_surface_residues(pose)

with open('profile', 'r') as fin:
    profile_lines = fin.readlines()

with open('profile', 'w') as fout:
    residue_index = 0
    for line in profile_lines:
        line = line.strip()
        if line.startswith('aa'):
            fout.write(line)
            fout.write('\n')
        else:
            residue_index += 1
            split_line = line.split()
            for i in range(0, len(split_line)):
                if residue_index in loop_residues:
                    i = i # do nothing
                elif residue_index in surface_residues:
                    if i != 0 and isfloat(split_line[i]):
                        split_line[i] = '%.2f' % ( float(split_line[i]) / 4 )
                    if i == 1:
                        split_line[i] = '%.2f' % 0.75
                else:
                    if i != 0 and isfloat(split_line[i]):
                        split_line[i] = '%.2f' % ( float(split_line[i]) / 8 )
                    if i == 1:
		                            split_line[i] = '%.2f' % 0.37
                #print(split_line[i])
                #if isfloat(split_line[i]):
                #    #print('is float')
                #    if float(split_line[i]) == 0 or float(split_line[i]) >= 4.0:
                #        split_line[i] = '%.2f' % 4.00
                #        #print('value was 0')
            line = "    ".join(split_line)
            fout.write(line)
            fout.write('\n')
