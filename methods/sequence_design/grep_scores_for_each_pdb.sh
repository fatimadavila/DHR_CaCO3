#!/usr/bin/env bash

mv score_file_for_all_pdbs.sc score_file_for_all_pdbs.sc.OLD
for dir_path in $(find run*/ -type d -name "X*")
do 
	for  silent_file_path in $(find $dir_path -maxdepth 1 -type f -name "results.out" )
	do 
		pdb_list=$(grep SCORE $silent_file_path | grep start | awk '{print $NF}')
		for pdb_name in $pdb_list
		do
			reps=$(grep SCORE $silent_file_path | grep $pdb_name | wc -l)
			for i in $(seq 1 $reps)
			do
				scores=$(grep SCORE $silent_file_path | grep $pdb_name | head -$i | tail -1)
				score_names=$(grep score $silent_file_path)
				echo "$scores $(dirname $silent_file_path) $(($i-1))" >> score_file_for_all_pdbs.sc
			done
		done
	done
done
echo "$score_names pdb_path rep" > tmp.file
cat score_file_for_all_pdbs.sc >> tmp.file
mv tmp.file score_file_for_all_pdbs.sc
