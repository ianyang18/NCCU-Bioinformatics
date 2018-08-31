type=(sechellia.taiwan sechellia.japan simulans.japan)

for t in ${type[@]}
do
	species=`echo $t | cut -c1-3`
	if [ ${species} == "sim" ];then
		gat=./data/annotation.table/dsim-all-r1.3.gtf
		fa_ref=./data/reference.genome/dsim-all-chromosome-r1.3.fasta
	else
		gat=./data/annotation.table/dsse-all-r1.3.gff
		fa_ref=./data/reference.genome/dsim-all-chromosome-r1.3.fasta
	fi

	for num in {1..6}
	do
		bam=${t}_${num}/accepted_hits.bam
		output=./cufflinks.output/${t}/${num}
		cufflinks -p 8 -N -o ${output} -G ${gat} -b ${fa_ref} -u ${bam}
	done
done
