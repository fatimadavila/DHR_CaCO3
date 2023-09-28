for a in $(seq 0 6 96); do
for b in $(seq 0 6 96); do
for c in $(seq -1.8 0.6 1.8); do
for d in $(seq -1.8 0.6 1.8); do
echo "/home/thuddy/th_home/TJ_rosetta_dev_backup/main/source/bin/rosetta_scripts.hdf5.linuxgccrelease -indexed_structure_store:fragment_store /home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h5 -s lattice8_7.pdb -parser:protocol des_test.xml -parser:script_vars z1="$c" z2="$d" rot1="$a" rot2="$b" @/home/thuddy/th_home/DHR/output/flags.txt -out:suffix _params_"$b"_"$c"_"$a"_"$d" "
done
done
done
done
