#################################################
#      Step 1. Download Raw Data From NCBI      #
#################################################

mkdir data
# download Drosophia Sechellia data collected from Taiwan
target_path=./data/sechellia.taiwan
mkdir $target_path
for num in {1952772..1952778}
do
  echo "download $target_path/SRR$num"
  wget -P $target_path ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP056/SRP056883/SRR$num/SRR$num.sra
done

# download Drosophia Sechellia data collected from Japan
target_path=./data/sechellia.japan
mkdir $target_path
for num in {1973486..1973491}
do
  echo "download $target_path/SRR$num"
  wget -P $target_path ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP057/SRP057153/SRR$num/SRR$num.sra
done

# download Drosophia Sechellia data collected from Japan
target_path=./data/simulans.japan
mkdir $target_path
for num in {1973492..1973497}
do
  echo "download $target_path/SRR$num"
  wget -P $target_path ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP057/SRP057154/SRR$num/SRR$num.sra
done

######################################################## 
#   Step 2. Convert raw data from SRA format to        #
#           paired-end FASTQ format using sratoolkit   #
########################################################

for target in sechellia.taiwan sechellia.japan simulans.japan
do
  for file in `ls ./data/${target}`
  do
    echo "processing ./data/${target}/${file}"
      ./sratoolkit.2.8.1/bin/fastq-dump --split-3 ./data/$target.fastq/$file
  done
  mv *.fastq ./data/$target.fastq/.
done

################################################################
#      Step 3. Download Genome Annotation Table from Flybase   #
################################################################

mkdir ./data/annotaion.table
# D. simulans
wget -P ./data/annotation.table ftp://flybase.net/genomes/Drosophila_sechellia/dsec_r1.3_FB2011_08/gff/dsec-all-r1.3.gff.gz
rm ./data/annotation/dsec-all-r1.3.gff.gz
# D. sechellia
wget -P ./data/annotation.table ftp://flybase.net/genomes/Drosophila_simulans/dsim_r1.3_FB2011_08/gff/dsim-all-r1.3.gff.gz
rm ./data/annotation/dsim-all-r1.3.gff.gz

#########################################################
#      Step 4. Download Reference Genome from Flybase   #
#########################################################

# D. sechellia
wget -P ./data/reference.genome ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r1.3_FB2011_08/fasta/dsim-all-chromosome-r1.3.fasta.gz
gzip -d dsim-all-chromosome-r1.3.fasta.gz
# D. simulans
wget -P ./data/reference.genome ftp://ftp.flybase.net/genomes/Drosophila_sechellia/dsec_r1.3_FB2011_08/fasta/dsec-all-chromosome-r1.3.fasta.gz
gzip -d dsec-all-chromosome-r1.3.fasta.gz
