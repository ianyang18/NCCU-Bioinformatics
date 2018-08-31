fname=(sec.japan.female sec.japan.male sec.taiwan.female sec.taiwan.male sim.japan.female sim.japan.male)
mkdir diff

for name1 in `eval echo {0..$(( ${#fname[@]}-2))}`
do
	for name2 in `eval echo {$(( $name1+1 ))..$(( ${#fname[@]}-1 ))}`
	do
		index1=${fname[$name1]}
		index2=${fname[$name2]}
		echo "############ compare ${index1} ${index2} ##############"	
		cuffdiff -p 23 -o diff/diff.${index1}_${index2} -u -N -L ${index1},${index2}  merge/merge_${index1}_${index2}/merged.gtf bam.data/${index1}.bam bam.data/${index2}.bam
		
	done
done
