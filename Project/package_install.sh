#############################################################
#   Required Package: 									    #
#		Tophat(v2.1.1), Bowtie(v2.2.8), Cufflinks(v2.2.1)   #
#		Samtools(v1.3.1), Sratoolkit(v2.8.1)			    #
#############################################################

mkdir package
# download Tophat(v2.1.1)
wget -P ./package http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -C ./package -zxvf ./package/tophat-2.1.1.Linux_x86_64.tar.gz
rm ./package/tophat-2.1.1.Linux_x86_64.tar.gz

# download Bowtie(v2.2.8)
wget -P ./package https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.8/bowtie2-2.2.8-linux-x86_64.zip/download
unzip ./package/download -d ./package
rm ./package/download

# download Cufflinks(v2.2.1)
wget -P ./package http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -C ./package -zxvf ./package/cufflinks-2.2.1.Linux_x86_64.tar.gz
rm ./package/cufflinks-2.2.1.Linux_x86_64.tar.gz

# download and install Samtools(v1.3.1)
wget -P ./package https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2/download
tar -C ./package -jxvf ./package/download
rm ./package/download
cd ./package/samtools-1.3.1
make
cd ..

# download Sratoolkit(v2.8.1)
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1/sratoolkit.2.8.1-ubuntu64.tar.gz
tar -zxvf sratoolkit.2.8.1-ubuntu64.tar.gz
rm sratoolkit.2.8.1-ubuntu64.tar.gz
mv sratoolkit.2.8.1-ubuntu64 sratoolkit.2.8.1

echo "--------------------------------------------------------------------------------"
echo "NOTICE: add path of Tophat, Bowtie, Cufflinks, Samtools, Sratoolkit to your PATH"
echo "--------------------------------------------------------------------------------"
