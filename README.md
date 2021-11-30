## 5' scL1seq

This repository contains python scripts for counting UMIs from 10x genomics 5' targeted 
single cell RNA-seq with 100+ base pair paired end reads. Instructions are given in this
readme.

### Step 1: Cell Ranger
5' scL1seq uses a custom Cell Ranger database to find LINE-1 aligning reads. You can
download these at the following links.
- Human: [https://slice.med.nyu.edu/data/5scL1seq\_cellranger\_index\_human.tar.gz](https://slice.med.nyu.edu/data/5scL1seq_cellranger_index_human.tar.gz)
- Mouse: [https://slice.med.nyu.edu/data/5scL1seq\_cellranger\_index\_mouse.tar.gz](https://slice.med.nyu.edu/data/5scL1seq_cellranger_index_mouse.tar.gz)

Then you can run cellranger count to align reads and count gene UMIs. Instructions to
obtain and run cellranger count can be found here.
- [Cell Ranger from 10x genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)

On a machine with 16 cores and 64gb of memory, cellranger count can be executed as follows:
```
cellranger count --sample=<sample name> --id=<output id> --transcriptome=/path/to/custom/index/hg38_L1_Alu --fastqs=/path/to/folder/with/raw/fastqs --localcores=16 --localmem=64
```

If combining multiple samples, we recommend running [cellranger aggr](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate) to downsample reads
and avoid artifacts that may result from samples sequenced to different depths.

### Step 2: Count UMIs
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
