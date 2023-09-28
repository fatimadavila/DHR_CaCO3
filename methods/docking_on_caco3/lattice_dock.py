#! /usr/bin/env python

# Rosetta modules
import pyrosetta
from pyrosetta import *
from rosetta import *
# from loadXMLandRun_w_params import MakeMoverFromXML

# Non-rosetta imports
import numpy as np
import argparse
import sys
import re

def MakeMoverFromXML(xmlfile, pose, script_vars):
   '''
   Main XML-parsing function. Returns mover.
   '''
   # Now actually create mover
   parser = protocols.rosetta_scripts.RosettaScriptsParser()
   
   # Necessary magic step:
   options = basic.options.process()
   #####################
   
   # Manage variable replacements:
   replaces = utility.vector1_string(0)
   ###DRH silent for now 
   for var in script_vars:
      replaces.append(var)
   ####################
   
   modified_pose = False
   tag = parser.create_tag_from_xml( xmlfile, replaces)
   in_mover = parser.generate_mover_for_protocol( pose, modified_pose, tag, options )

   return in_mover

def split_chain_and_residue(ResidueIds):
  ReturnList = []
  for Res in ResidueIds:
    ReturnList.append( re.sub(r'(\d+)[A-Z]+',r'\1',Res) )
    ReturnList.append( re.sub(r'\d+([A-Z]+)',r'\1',Res) )
  return tuple(ReturnList)

def calculate_direction(c1, c2):
   ref_coord = c1 - np.array([1,0,0])
   spin_disp = c2 - c1
   # print spin_disp
   direction = np.degrees(np.arctan(spin_disp[0]/spin_disp[1]))
   if spin_disp[0]<0 and spin_disp[1]<0:
      direction += -180
   elif spin_disp[1]<0:
      direction += 180
   return direction

def main():
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('-pdb', type=str, help='input pdb', required=True)
   parser.add_argument('-xml', type=str, help='docking xml to apply as mover', default='dock_on_surfaces.xml')
   parser.add_argument('-state', type=str, help='ON or OFF surface', required=True)
   parser.add_argument('-tag', type=str, help='tag for this thread', default='1A')
   parser.add_argument('-cst', type=str, help='cst file to load as script var', default=None)
   parser.add_argument('-nstruct', '-n', type=int, help='number of iterations for this thread', default=1)
   parser.add_argument('-dir_res1', type=int, help='1st of 2 residues to determine angle of the protein relative to Y-direction in XY ', required=True)
   parser.add_argument('-dir_atom1', type=str, help='1st of 2 atoms to determine angle of the protein relative to Y-direction in XY ', default='CA')
   parser.add_argument('-dir_res2', type=int, help='2nd of 2 residues to determine angle of the protein relative to Y-direction in XY ', required=True)
   parser.add_argument('-dir_atom2', type=str, help='2nd of 2 atoms to determine angle of the protein relative to Y-direction in XY ', default='CA')
   parser.add_argument('-params', '-extra_res_fa', type=str, nargs="+", help="params files to load" , default = [] )
   parser.add_argument('-repeat', '-rep', '-r', type=str, help="residues per repeat (currently this does not do anything)" , default=-1 )
   parser.add_argument('-dump', type=int, help="If set to 1, dump all pdbs. FOR DEBUGGING ONLY!!! (default= 0)", default= 0 )
   parser.add_argument('-mute', type=int, help='If set to 1, print all rosetta comments to screen (default= 0)', default=1 ) 
   args = parser.parse_args()

   if args.state == 'off': args.state = 'OFF'
   if args.state == 'on': args.state = 'ON'
   assert args.state in ['OFF', 'ON'], '-state  arguement must be set to OFF or ON (surface)'
   
   if len(args.params): params = '-extra_res_fa '+' '.join(args.params)
   else: params = ''

   if args.mute: pyrosetta.init(extra_options = "-corrections::beta_nov16 "+params+" -mute all -chemical:exclude_patches Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm    # Excluded to improve speed and memory usage " )
   else: pyrosetta.init(extra_options = "-corrections::beta_nov16 "+params+"  -chemical:exclude_patches Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm    # Excluded to improve speed and memory usage " )

   score_function = create_score_function("beta_nov16") 
   coord_cst_function = create_score_function("empty")
   coord_cst_function.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint, 1.0) 

   thread_name = '{0}_thread_{1}'.format(args.pdb.split('.pdb')[0], args.tag)
   low_nrg = 0

   silent_file_object = core.io.silent.SilentFileData( thread_name+".silent", False, False, "binary", core.io.silent.SilentFileOptions())

   with open(thread_name+'_score.tsv', 'w') as thread_log_file:
      print( '\t'.join(['total_score', 'spin', 'com_z', 'min_z', 'id_number', 'state', 'min_dist', 'med_dist', 'max_dist']), end="\n", file=thread_log_file)

   pose = Pose()
   pose_from_file(pose, filename=args.pdb)

   # Extract chain A (the surface) to shift residue numbers from position in protein to position in full system 
   surface_chain = pose.clone().split_by_chain()[1]

   res_in_surface = surface_chain.size()
   script_vars = ["surf_to_prot={0}".format(res_in_surface)] 

   if args.state == 'OFF':
      min_dist = 0
      med_dist = 0
      max_dist = 0
   
   mover_from_XML = MakeMoverFromXML( args.xml, pose, script_vars )

   for i in range(1, args.nstruct+1):
      num_string = str(i).rjust(4, '0')
      print (num_string)
      pose_instance = pose.clone()

      mover_from_XML.apply(pose_instance)

      # calculate direction: angle relative to the Y-axis in the XY-plane
      direction_coord1 = np.array(list(pose_instance.residue(args.dir_res1 + res_in_surface ).xyz(args.dir_atom1)))
      direction_coord2 = np.array(list(pose_instance.residue(args.dir_res2 + res_in_surface ).xyz(args.dir_atom2)))
      direction = calculate_direction(direction_coord1, direction_coord2)
      round_direction = int(round(direction))

      rosetta_score = score_function(pose_instance)
      pdbtag = '{0}_{1}_{2}spin'.format(thread_name, num_string, round_direction)
      
      if args.dump:
         pose_instance.dump_file( '{0}.pdb'.format(pdbtag) )

      if rosetta_score < low_nrg:
         # pose_instance.dump_file( '{0}_LOW_ENERGY.pdb'.format(thread_name) )
         low_nrg = rosetta_score
  
      'EXTRACT CHAIN B (protein) FOR SILENT FILE'
      protein_chain = pose_instance.split_by_chain()[2]
      # C-alpha carbon center of mass
      center_of_mass_z = float(tuple(pyrosetta.rosetta.protocols.geometry.center_of_mass( protein_chain, 1, protein_chain.size() ))[2])

      min_z = 1000
      for r in range(1, protein_chain.size()+1):
          z = float(list(protein_chain.residue(r).xyz('CA'))[-1])
          if z < min_z:
              min_z = z

      sfd_out = silent_file_object.create_SilentStructOP()
      sfd_out.fill_struct(protein_chain, '{0}_{1}_{2}spin'.format(thread_name, num_string, round_direction) )
      silent_file_object.write_silent_struct(sfd_out, filename=thread_name+".silent")

      if args.state == 'ON':
         calpha_coords = [ list(protein_chain.residue(resi).xyz('CA')) for resi in range(1, protein_chain.size()+1)]
         calpha_z_coords = [coord[-1] for coord in calpha_coords]

         min_dist = round(min(calpha_z_coords), 3)
         med_dist = round(np.median(calpha_z_coords), 3)
         max_dist = round(max(calpha_z_coords), 3)        

      with open(thread_name+'_score.tsv', 'a') as thread_log_file:
         print('\t'.join([str(x) for x in [round(rosetta_score, 3), round(direction, 3), round(center_of_mass_z, 3), round(min_z, 3), num_string, args.state, min_dist, med_dist, max_dist]]), end="\n", file=thread_log_file)


if __name__ == "__main__":
   sys.exit(main())

