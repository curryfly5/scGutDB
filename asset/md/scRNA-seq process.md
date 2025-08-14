## scRNA-seq process

### 一、上游处理

公司的下机数据为fastq格式，如果是10X平台，我们只需用cellranger count功能进行数据统计：完成细胞和基因的定量，产生我们用来做各种分析的基因表达矩阵。

```shell
#!/bin/bash
cellranger count --id=CE_180 --sample=A180_CE_1-1,A180_CE_1-2 --fastqs=/disk191_3/LabData/Pig/2_SingleCellRNA/Pro4_JH/X101SC22093502-Z01-J001-B1-1_10X_release_20221006/Rawdata/A180_CE_1 --transcriptome=/disk201/zz/JH/scRNA/duroc_reference_genome/SsDuroc --lanes=2 --localcores=64 --localmem=100
```

cellranger count也可以通过--sample命令和--lane命令对特定样品进行处理

```tips
# 1.所有样本
--fastqs=/PATH/TO/PROJECT_FOLDER
# 2. 某个样本的所有lanes数据
--fastqs=/PATH/TO/PROJECT_FOLDER \
--sample=SAMPLENAME
# 3. 某个样本的某个lane
--sample=SAMPLENAME \
--fastqs=/PATH/TO/PROJECT_FOLDER \
--lanes=1
# 4. 多个sample（补测数据可以使用这个命令）
--sample=A180_CE_1-1,A180_CE_1-2
```

### 二、下游分析

`cellranger count`计算的结果只能作为初步观测的结果，如果需要进一步分析聚类细胞，还需要进行下游分析，一般使用seurat（R）和scanpy（python）

具体分析脚本详见scanpy标准流程

工具：jupyter notebook

连接remote后可以调用我的公共分析环境，里面会有

#### 1、质控-标准化-归一化-降维-聚类

标准流程地址：`/disk212/yupf/database/scRNA-seq/scanpy_analysis/scanpy Standardized Process Script.ipynb`

《scanpy Standardized Process Script》中的前三部分

+ 预处理
+ cellcycle correct
+ reduction

### 2、细胞谱系注释

得到降维聚类的对象后下一步便是根据差异基因来鉴定细胞类型，细胞类型鉴定的思路是先鉴定大类的细胞谱系然后再根据细胞谱系进一步做亚群鉴定

经文献阅读和数据分析总结，肠道细胞基本可以分为以下七个谱系：

+ B系淋巴细胞谱系
+ T系/先天淋巴细胞谱系
+ 髓系免疫细胞谱系
+ 上皮细胞谱系
+ 内皮细胞谱系
+ 神经细胞谱系
+ 间充质细胞谱系

```python
#细胞谱系预测
first_marker={'B/plasma':['CD19', 'CD79B', 'MS4A1', 'JCHAIN'], # B/plasma cell genes
'T/ILC':['CD3E', 'CD3G', 'CD52', 'ZAP70'], # T/ILC genes
'myeloid lineage':['ENSSSCG00000028461', 'CD68', 'CSF2RB', 'ICAM1'], # non-lymphoid immune cell genes
'epithelial cell':['EPCAM', 'KRT8', 'ACTA2', 'COL18A1'],# epithelial cell genes
'endothelial cell':['PECAM1','CDH5'],# endothelial cell
'neuronal':['SOX10','ELAVL4','GAL'],# neuronal cell
'mesenchymal':['COL3A1','CALD1','COL1A1', 'COL1A2','MFAP4','ACTA2']# mesenchymal cell
}
cellt=['B/plasma','T/ILC','myeloid lineage','epithelial cell','endothelial cell','neuronal','mesenchymal']
for i in cellt:
    Glist=first_marker[i]
    for a in Glist:
        if a not in adata.var_names:
            first_marker[i].remove(a)
            print(a+' is not in adata.var_names')
sc.pl.dotplot(adata, first_marker, 'leiden',dendrogram=True, cmap='Reds')#此处分辨率可以选择较大的，划分较为精细
```

然后可以根据谱系结果进行提取子集做后续注释

#### 3、细胞自动注释

使用auto.py脚本

脚本地址`/disk212/yupf/database/scRNA-seq/scanpy_analysis/auto.py`

```shell
cd /disk212/yupf/database/scRNA-seq/scanpy_analysis/
nohup python auto.py --h5ad xxxx -o xxxx -c xxxx -m xxxx >auto.log 2>&1  &
```

参数说明

+ --h5ad'/'-h'

  type=str, default=None 

  说明：降维聚类后保存的h5ad文件

+ --outdir'/'-o'

  type=str, default='./' 

  说明：输出文件地址，会输出名为adata_prediction.h5ad的文件，其adata对象的obs中会包含一个名为’majority_voting‘的条目，包含了预测的cluster类别信息

+ --overclustering'/'-c'

  type=str, default='cluster_0.9'

  说明：adata.obs中'predicted_labels'是每个细胞的label信息，想和leiden聚类cluster进行拟合，需要给定参考的分辨率，才能进行majority_voting，一般选取分辨率较大值即可

+ --model'/'-m'

  type=str, default='/disk212/yupf/database/scRNA-seq/scanpy_analysis/atlas/ref_il/JH/JH.pkl'

  说明：选择细胞类型预测所使用的模型——如果要进行初步鉴定，则使用JH.pkl（默认值）；如果是选取T/B/M等谱系鉴定亚群，则选用其他模型

  | **model**             | **description**        |                                                              |
  | --------------------- | ---------------------- | ------------------------------------------------------------ |
  | model_Myeloid_ref.pkl | 髓系免疫细胞谱系参考集 | /disk212/yupf/database/scRNA-seq/scanpy_analysis/atlas/ref_il/Myeloid_ref/model_Myeloid_ref_22.11.27.pkl |
  | model_T_ref.pkl       | T/先天淋巴谱系参考集   | /disk212/yupf/database/scRNA-seq/scanpy_analysis/atlas/ref_il/T_ref/model_T_ref_22.11.27.pkl |
  | model_B_ref.pkl       | B淋巴谱系参考集        | /disk212/yupf/database/scRNA-seq/scanpy_analysis/atlas/ref_il/B_ref/model_B_ref_22.11.27.pkl |
  | JH.pkl                | 金华细胞类型参考集     | /disk212/yupf/database/scRNA-seq/scanpy_analysis/atlas/ref_il/JH/JH.pkl |

其余模型待训练——endothelial、neuronal、mesenchymal

#### 3、使用marker gene验证调整

详见《scanpy Standardized Process Script》

+ 细胞类型鉴定
+ 进一步鉴定

细胞类型鉴定可先采用celltypist训练模型的自动鉴定，辅助以dotplot、feature plot、overlap plot验证调整

1.feature plot中可以绘制各细胞类型的marker gene分布，marker gene包括以往文献收集的，也包括数据库公布的TOP50的DEG，做为参考

2.dotplot为以往文献中收集的marker gene

3.panglao、atlas、zju、human_marker分别为panglao数据库、nc文章、汪海峰老师文章、cmgh文章的overlap分析，均可以提供一定参考，但不准确，需要综合看，供参考

### 三、高级分析

#### 1、细胞轨迹分析

工具：monocle2（R）

#### 2、拟时序分析

工具：scVelo（python）

####  3、细胞通讯分析

工具：CellChat（R）

#### 4、GSEA与GSVA分析

工具：GSEA（R）、GSVA（R）

#### 5、转录调控网络分析（SCENIC）

工具：AUcell（R）、RcisTarget（R）、GENIE3（R）

#### ....









