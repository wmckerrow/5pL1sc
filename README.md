## 5' scL1seq

This repository contains python scripts for counting UMIs from 10x genomics 5' targeted 
single cell RNA-seq with 100+ base pair paired end reads. Instructions are given in this
readme.

### Step 1: Build the custom cellranger index (human)
Building the cellranger custom index, will require bedtools and cellranger. Both are
available for install via anaconda:
-https://anaconda.org/bioconda/bedtools
-https://anaconda.org/hcc/cellranger

Before beginning you will need the ucsc genome browser version of the hg38 human genome,
which can be downloaded as follows:
```
curl -L -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
zcat hg38.fa.gz > hg38.fa
```

You will also need hg38 gene annotations in gtf format, which can be downloaded as follows:

```
curl -L -O https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
zcat Homo_sapiens.GRCh38.108.gtf.gz > Homo_sapiens.GRCh38.108.gtf
sed -i 's/^/chr/g' Homo_sapiens.GRCh38.108.gtf
sed -i 's/^chr#/#/g' Homo_sapiens.GRCh38.108.gtf
```

You can then use the provided shell script to build the index:
```
bash make_custom_cellranger_reference.sh /path/to/hg38.fa /path/to/5pL1sc/L1_annotation/L1HS_and_PA.bed /path/to/5pL1sc/L1_annotation/L1HS_and_dfam_L1PA.fa /path/to/Homo_sapiens.GRCh38.108.gtf L1HS_L1PA_seperated_hg38
```

This will create the custom cellranger index in a new directory called L1HS\_L1PA\_seperated\_hg38

### Step 1 v2: Build the custom cellranger index (mouse)
### Step 1: Build the custom cellranger index (human)
Building the cellranger custom index, will require bedtools and cellranger. Both are
available for install via anaconda:
-https://anaconda.org/bioconda/bedtools
-https://anaconda.org/hcc/cellranger

Before beginning you will need the ucsc genome browser version of the mm39 mouse genome,
which can be downloaded as follows:
```
curl -L -O http://hgdownload.cse.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
zcat mm39.fa.gz > mm39.fa
```

You will also need mm39 gene annotations in gtf format, which can be downloaded as follows:

```
curl -L -O https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz
zcat Mus_musculus.GRCm39.108.gtf.gz > Mus_musculus.GRCm39.108.gtf
sed -i 's/^/chr/g' Mus_musculus.GRCm39.108.gtf
sed -i 's/^chr#/#/g' Mus_musculus.GRCm39.108.gtf
```

You can then use the provided shell script to build the index:
```
bash make_custom_cellranger_reference.sh /path/to/mm39.fa /path/to/5pL1sc/L1_annotation/L1Md.bed /path/to/5pL1sc/L1_annotation/L1MdI.Consensus.fa /path/to/Mus_musculus.GRCm39.108.gtf L1Md_seperated_mm39
```

This will create the custom cellranger index in a new directory called L1Md\_seperated\_mm39


### Step 2: Cellranger

Then you can run cellranger count to align reads and count gene UMIs. Instructions to
obtain and run cellranger count can be found here.
- [Cell Ranger from 10x genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)

On a machine with 16 cores and 64gb of memory, cellranger count can be executed as follows:
```
cellranger count --sample=<sample name> --id=<output id> --transcriptome=/path/to/custom/index/L1HS_L1PA_seperated_hg38[L1Md_seperated_mm39] --fastqs=/path/to/folder/with/raw/fastqs --localcores=16 --localmem=64
```

If combining multiple samples, we recommend running [cellranger aggr](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate) to downsample reads
and avoid artifacts that may result from samples sequenced to different depths.

### Step 3: Count UMIs
LINE-1 UMIs are counted using:
- count\_properpairUMIs\_in\_range\_by\_nM.py
To count human LINE-1 in a single sample, you need only specify the bam file output by
Cell Ranger:
```
python count_properpairUMIs_in_range_by_nM.py --bam 5p_sc_sample/outs/possorted_genome_bam.bam > L1HS_properUMI_counts.txt
```

To count mouse LINE-1, you will need specify relevant ranges for the LINE-1 reads to fall.
We used:
  * --read1\_range L1MdTf\_I\_5end:1-600 --read2\_range L1MdTf\_I\_5end for L1MdTf
  * --read1\_range L1MdGf\_I\_5end:1-550 --read2\_range L1MdGf\_I\_5end for L1MdGf
  * --read1\_range L1MdA\_I\_5end:1-750 --read2\_range L1MdA\_I\_5end for L1MdA

If you intend to do differential expression, we recommend that you use cellranger aggr
to merge samples and down sample reads to a constant number per cell. The amount of
down sampling done by aggr can be found in the summary.json output file. You can then
pass that amount of down sampling forward using the --downsample_to option.

The output will be a two column tab delimited table, where column 1 is the cell barcode
and column two is the number of L1Hs UMIs with the cell barcode. Cells with 0 L1Hs UMIs
are not reported.

### Optional: add L1Hs UMIs to seurat object

With a Seurat object (seurat\_10x5p) and the output from above (L1HS\_properUMI\_counts.txt)
The UMIs can be added to the Seurat object in R as follows:
```
seurat_10x5p_L1_UMIs = read.table("~/Documents/LINE1_projects/sc_LINE1/vdj_v1_hs_nsclc_5gex/L1HS_properUMI_counts.txt",sep='\t')
seurat_10x5p_L1_UMIs = seurat_10x5p_L1_UMIs[seurat_10x5p_L1_UMIs[,1]%in%colnames(seurat_10x5p),]
L1_count = rep(0,ncol(seurat_10x5p))
names(L1_count) = colnames(seurat_10x5p)
L1_count[seurat_10x5p_L1_UMIs[,1]] = seurat_10x5p_L1_UMIs[,2]
seurat_10x5p$L1_count = L1_count
seurat_10x5p$L1_lognorm = log(seurat_10x5p$L1_count/seurat_10x5p$nCount_RNA*10^4+1)
```

### Important notes
- The Cell Ranger output will quantify an L1Hs transcript. This quantification does not
include the additional filtering in step 2. Only consider with extreme caution.
- The Cell Ranger index also includes Alu consensus that could be used to quantify Alu.
We chose not to include any Alu quantifications as we did not have a sample with a known
pattern of Alu expression that could serve as test case. Again extreme caution should be
excercised when using an Alu related output.
