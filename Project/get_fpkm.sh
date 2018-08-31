type=(dsec_japan_female dsec_japan_male dsec_taiwan_female dsec_taiwan_male dsim_japan_female dsim_japan_male dmel_female dmel_male)

for t in ${type[@]}
do
	while read line;do
		species=`echo $t | awk 'BEGIN {FS="_"} {print $1}'`
		country=`echo $t | awk 'BEGIN {FS="_"} {print $2}'`
		sex=`echo $t | awk 'BEGIN {FS="_"} {print $3}'`

		# gene, dmel, dsec, dsim
		gene=`echo $line | awk 'BEGIN {FS=","} {print $1}'`
		dmel=`echo $line | awk 'BEGIN {FS=","} {print $2}'`
		dsec=`echo $line | awk 'BEGIN {FS=","} {print $3}'`
		dsim=`echo $line | awk 'BEGIN {FS=","} {print $4}'`

		# calculate dsec
		if [ "${species}" == "dsim" ]; then
			grep "$dsim" ./cufflinks.output/${t}/gene.fpkm_extract > tmp
		elif [ "${species}" == "dsec" ]; then
			grep "$dsec" ./cufflinks.output/${t}/gene.fpkm_extract > tmp
		else
			grep "$dmel" ./cufflinks.output/${t}/gene.fpkm_extract > tmp
		fi

		score=""
		while read line;do
			n=`echo ${line} | awk '{print $10}'`
			score="${score} ${n}"
		done < tmp

		echo "${t}*${gene}*${score}"
		
	done < gene_id.csv
done
