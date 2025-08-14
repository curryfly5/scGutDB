We can investigate the association between cell types and complex traits by calculating enrichment scores, such as scDRS.
all file you need is genotype, gwas summary, single cell(anndata.h5ad)


```python
import pandas as pd
import os
import scanpy as sc
import scdrs
import re
```

### 1.Get gene-based analysis from GWAS summary：magma analysis

##### (optional) first make gene annotation file from gtf, you can also download () , it's advisable to remake gene annotation for specific genotype. 


```python
### make gene loc file
gtf_file = "./genes.gtf"
output_file = "./pig.gene.loc"

with open(gtf_file, "r") as f:
    with open(output_file, "w") as out:
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    attributes = dict(re.findall(r"(\w+)\s\"([\w\.]+)\";", fields[8]))
                    chrom = fields[0]
                    start = fields[3]
                    end = fields[4]
                    gene_id = attributes["gene_id"]
                    gene_name = attributes.get("gene_name", "NA")
                    out.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{gene_name}\n")

```

##### 1.1 The previous step generated the pig.gene.loc file. Enter the command `magma --annotate` to generate the "pig.genes.annot" file, and specify the snpid （plink format bim file）.


```python
magma --annotate window=10,10\
    --snp-loc {snpid} \
    --gene-loc ./pig.gene.loc \
    --out ./pig
```

##### 1.2 Prepare GWAS summary statistics, rsid and pval are needed. (if meta-GWAS provided and a column of sample sizes per SNP is available in `[PVAL_FILE]`,  `ncol=[N_COLUMN_NAME]` for the --pval flag (instead of `N=[N]`).  )


```python
magma 
--bfile {genotype} \ ### plink format
--pval {gwas_summary} use=rs,p_wald N=2797 \ ### rs,pval were column names of GWAS summary statistics
--gene-annot {gene_anno} \ 
--out {outdir}
```

This step generates the out.genes.out file, which contains the information required for gene-based analysis.

##### 1.3 Extract magma zscore top1000 genes


```python
### 此处你可以根据实际情况做一些修改，目的是提取出一列基因，一列pval，而且顺序是按zscore排序的前1000个
df=pd.read_csv('./out.genes.out',delim_whitespace=True)
sorted_df = df.sort_values(by='ZSTAT', ascending=False)
new_df=sorted_df.iloc[:1000,:]
new_df.to_csv(f'./{Trait}_magma_gene_ztop1000.tsv',sep='\t',index=False)
```

### 2.Enrichment analysis:scDRS

##### 2.1 Imputation of gene expression levels for single-cell data.


```python
### you can download imputation.py from github
python imputation.py \
-i {filename} \ 
--h5ad ./singlecell.h5ad \ 
-o ./singlecell_impute
# i: Specify the custom filename
# --h5ad: Input the h5ad file for single-cell analysis data
# -o: Output directory
# After running, it will generate a file named {filenamei}_impute.h5ad, which is the h5ad file with imputed expression levels.
```

##### 2.2 To generate the gene set (gs) file using `scdrs munge-gs`, which creates a gene weight file based on the output from the previous step `magma_gene_ztop1000.tsv`, you can follow these steps:


```python
scdrs munge-gs \
    --out_file ./{Trait}.gs \
    --pval-file ./{Trait}_magma_gene_ztop1000.tsv \
    --weight zscore \
    --n-max 1000 \
    --fdr 0.05

```

##### 2.3 To generate enrichment scores for each cell in a single-cell atlas using `scdrs compute-score`, you can follow these steps:

In the command, the --h5ad option should include the annotated single-cell data in h5ad format, while the --gs-file option should specify the gene set file obtained from the previous step. 


```python
scdrs compute-score \
    --h5ad-file ./{filename}_impute.h5ad \
    --h5ad-species pig \
    --gs-file ./{Trait}.gs \
    --gs-species pig \
    --out-folder ./compute_score
```

##### 2.4 To perform statistical tests using `scdrs perform-downstream`, you can use the following command structure:

`--score-file` scores_file.h5ad: Input the score file obtained from the previous step, which contains the enrichment scores for each cell.
`--group-analysis` cell_type_annotation: Specify the name of the cell type annotation used during single-cell analysis, corresponding to an entry in adata.obs.


```python
scdrs perform-downstream \
        --h5ad-file ./{filename}_impute.h5ad \
        --score-file ./{Trait}.full_score.gz \
        --out-folder ./subc_score \
        --group-analysis Lineage \ # one factor from adata.obs
        --flag-filter-data True \
        --flag-raw-count True
```

### 3.Visualization

##### 3.1 Umap visualization


```python
adata=sc.read_h5ad(f'./{filename}_impute.h5ad')
dict_score = {
    trait: pd.read_csv(f"./{Trait}.full_score.gz", sep="\t", index_col=0)
    for trait in df.index
}

for trait in dict_score:
    adata.obs[trait] = dict_score[trait]["norm_score"]

sc.set_figure_params(figsize=[2.5, 2.5], dpi=150)
sc.pl.umap(
    adata,
    color="Lineage",
    ncols=1,
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
)

sc.pl.umap(
    adata,
    color=dict_score.keys(),
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    s=10,
    hspace=0.5,
    wspace=0.5)
```

##### 3.2 Statistical tests visualization


```python
dict_df_stats = {
    trait: pd.read_csv(f"./{Trait}.scdrs_group.Lineage", sep="\t", index_col=0)
    for trait in ['Trait']
}

fig, ax = scdrs.util.plot_group_stats(
    dict_df_stats={
        trait: df_stats
        for trait, df_stats in dict_df_stats.items()
    },
    plot_kws={
        "vmax": 0.2,
        "cb_fraction":0.12
    }
)
```
