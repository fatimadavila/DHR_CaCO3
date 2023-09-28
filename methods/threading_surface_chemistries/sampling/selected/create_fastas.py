import glob
import pyrosetta

pyrosetta.init()

pdbs = glob.glob('PEP*.pdb')
pdbs_sorted = sorted(pdbs)

name_seq_dict = {}

for p in pdbs_sorted:
    pose = pyrosetta.pose_from_pdb(p)
    seq = pose.sequence()
    rep_len = int(len(seq) / 4)  # Most sequences here have four repeats
    seq_1 = seq[:rep_len]  # Get first repeat
    seq_2to4 = seq[rep_len: rep_len * 2]  # Get second repeat
    seq_5 = seq[rep_len * 3:]  # Get last repeat
    prop_sequence = seq_1 + seq_2to4 * 3 + seq_5  # Create fasta with propagated inner repeats (5 repeats)
    final_seq = prop_sequence + 'SGGSGGENLYFQGS'
    name = p.strip('.pdb')
    name_seq_dict[name] = final_seq

with open('selected_5rep_pep.fa', 'w') as fout:
    for k, v in name_seq_dict.items():
        fout.write('>' + k + '\n')
        fout.write(v + '\n')
