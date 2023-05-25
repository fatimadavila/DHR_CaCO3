#! /usr/bin/env python3

import glob
# import os

# ori_dir = os.getcwd()

pdb_dir_list = glob.glob('/net/scratch/fadh/digs_submission/run_directory/straight_DHR_generation/run_0_0/X_h*/results/*pdb')
dir_out = '/net/scratch/fadh/digs_submission/run_directory/straight_DHR_generation/'
out_filename = 'jobs.file.design'

with open (dir_out+'/'+out_filename, 'w') as fout:
	for pdb_dir in pdb_dir_list:
		path = pdb_dir.split('/start')[0]
		pdb_name = pdb_dir.split('results/')[1]
		rep_info = pdb_dir.split('/')[8]
		h1 = float(rep_info.split('_')[1].strip('h'))
		l1 = float(rep_info.split('_')[2].strip('l'))
		h2 = float(rep_info.split('_')[3].strip('h'))
		l2 = float(rep_info.split('_')[4].strip('l'))
		rep_len_plus10 = str('%.0f' % (h1+l1+h2+l2+10))
		line = 'cd '+path+' ; bash /home/fadh/biomineralization/newDHRs/DHR_hematite/design_2/cmd '+pdb_name+' '+rep_len_plus10
		print(line)
		fout.write(line+'\n')
