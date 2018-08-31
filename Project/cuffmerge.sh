fname=(sec.japan.female sec.japan.male sec.taiwan.female sec.taiwan.male sim.japan.female sim.japan.male)
mkdir manifest

for name1 in `eval echo {0..$(( ${#fname[@]}-2))}`
do
	for name2 in `eval echo {$(( $name1+1 ))..$(( ${#fname[@]}-1 ))}`
	do
		index1=${fname[$name1]}
		index2=${fname[$name2]}

		`ls $( pwd )/gtf.data/gtf.${index1}/gtf.${index1}.1/transcripts.gtf >> manifest/manifest_${index1}_${index2}`
		`ls $( pwd )/gtf.data/gtf.${index1}/gtf.${index1}.2/transcripts.gtf >> manifest/manifest_${index1}_${index2}`
		`ls $( pwd )/gtf.data/gtf.${index1}/gtf.${index1}.3/transcripts.gtf >> manifest/manifest_${index1}_${index2}`
		`ls $( pwd )/gtf.data/gtf.${index2}/gtf.${index2}.1/transcripts.gtf >> manifest/manifest_${index1}_${index2}`
		`ls $( pwd )/gtf.data/gtf.${index2}/gtf.${index2}.2/transcripts.gtf >> manifest/manifest_${index1}_${index2}`
		`ls $( pwd )/gtf.data/gtf.${index2}/gtf.${index2}.3/transcripts.gtf >> manifest/manifest_${index1}_${index2}`

		echo "start cuffmerge...."
		if [ "${index1:0:3}" == "${index2:0:3}" ]; then
			if [ "${index1:0:3}" == "sec" ]; then
				ref_seq=ref_seq/dsec-all-chromosome-r1.3.fasta
			else
				ref_seq=ref_seq/dsim-all-chromosome-r1.3.fasta
			fi
		else
			ref_seq=ref_seq/dsim_sec-all-chromosome-r1.3.fasta
		fi
		
		output=merge/merge_${index1}_${index2}
		cuffmerge -p 20 -o ${output} -s ${ref_seq} manifest/manifest_${index1}_${index2}
	done
done

