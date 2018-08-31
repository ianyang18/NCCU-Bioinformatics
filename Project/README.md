Tools Installed
---------------
>Download tools and install them.

### Tools
- TopHat(v2.1.1)
- Bowtie2(v2.2.8)
- Cufflinks(v2.2.1)
- Sratoolkit(v2.8.1) 
- Samtools(1.3.1)

```
bash package_install.sh
```

Data Sources
------------
### Datasets

- Raw data (D.Sec.Taiwan/ D.Sec.Japan/ D.Sim.Japan)
- Reference Genome (D.Sec_r1.3/ D.Sim_r1.3)
- Reference Genome Annotation (D.Sec_r1.3/ D.Sim_r1.3)

### Description
>Download raw data from NCBI and convert its from SRA format to FASTQ format. Also, download reference genome and its annotation table from flybase.


```
bash data_download_convert2fastq.sh
```

WorkFlow
--------
### Step1.
>Using Bowtie2 to index the reference genome and using TopHat to align RNA-seq reads to genome.

```
bash tophat_mapping.sh
```

### Step2.

>Assembling RNA-seq into transcriptomes and quantifing the gene expression by FPKM.

```
bash cufflnks.sh
```

### Option(Step3.

>Merging assemblies into a master transcriptome.

```
bash cuffmerge.sh
```

### Option(Step4.

>Comparing expression levels of genes and transcripts in RNA-seq.


```
bash cuffdiff.sh
```

### Step5.

>Differentially Gene Expression (FPKM)

```
bash get_fpkm.sh
```





