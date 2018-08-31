type=(sechellia.taiwan sechellia.japan simulans.japan)

for t in ${type[@]}
do
	# reference genome indexing
	species=`echo $t | cut -c1-3`
	if [ ${species} == "sec"  ]; then
		reference=./data/annotation.table/dsim-all-chromosome-r1.3.fasta
		transcript_annotation=./data/annotation.table/dsec-all-r1.3.gtf
	else
		reference=./data/annotation.table/simulans.genome/dsim-all-chromosome-r1.3.fasta
		transcript_annotation=./data/annotation.table/dsim-all-r1.3.gff	
	fi
	reference_index=./data/genome.index/${species}_genome_bowtie2
	bowtie2-build ${reference} ${reference_index}

	n=0
	for f in `ls ./data/${t}.fastq/*_1.fastq`
	do
		fnumber=`echo $f | awk 'BEGIN {FS="_"} END {print $1}'`
		reads_1=${t}.fastq/${fnumber}_1.fastq
		reads_2=${t}.fastq/${fnumber}_2.fastq
		n=$(( $n+1 ))
		output=./data/tophat_mapping/${t}/simulans_$n

		tophat -p 20 -r 100 -G ${transcript_annotation} -o ${output} ${reference_index} ${reads_1} ${reads_2}
	done
done
