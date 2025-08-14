# scGutDB  
**A comprehensive, interactive and user-friendly platform for exploring intestinal single-cell RNA-seq datasets across pig, human and mouse.**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)  
[![Platform](https://img.shields.io/badge/platform-Web-blue)](http://alphaindex.zju.edu.cn:8501/)  

---

## 📌 Overview
**scGutDB** is an online exploration platform dedicated to fully utilizing intestinal single-cell sequencing resources. It hosts the **Integrated Pig Gut Cell Atlas (IPGCA)** and additional human and mouse intestinal scRNA-seq datasets. With interactive visualizations and an intuitive query interface, scGutDB enables researchers to query the DEGs or specific genes in each atlas leverl datasets. And we provide data download including **DEG tables**, **three type of eQTL**, **IPGCA resource**, and **auto annotation pkl files**.

---

## 🚀 Key Features
- **Cross-species atlas** – Integrated single-cell landscapes of **human**, **pig** and **mouse** intestines.  
- **Interactive exploration** – Select species, cell lineages and sub-clusters; dynamically generate UMAPs, violin/dot plots and heatmaps.  
- **Gene-centric queries** – Rapidly search gene expression patterns across cell types and conditions.  
- **Community sharing** – Publish and share custom analysis pipelines directly on the platform.

---

## 🧬 Datasets
| Species | Tissue | # Cells | # Genes | Source |
|---------|--------|---------|---------|--------|
| Human   | Small & Large Intestine | 1.2 M | 24 k | From gutcellatlas.org |
| Pig     | Small & Large Intestine  | 1.2 M | 22 k | IPGCA |

*Visit the [live portal](http://alphaindex.zju.edu.cn:8501/) for the most recent numbers.*

---

## 🛠️ Quick Start
1. **Open the portal**  
   Navigate to [http://alphaindex.zju.edu.cn:8501/](http://alphaindex.zju.edu.cn:8501/) in any modern browser. No installation required.

2. **Select page**  
   If you want to query DEGs or marker of specific cell type, select **Marker query** page
   If you want to query Gene in IPGCA, select **Explore IPGCA** page, and choose the dataset **Human** or **Pig**
   If you want to query whether specific variant is an eQTL or ct-ieQTL(cell type interaction eQTL) or decon-eQTL(deconvolutinal cell type eQTL), select **eQTLs query** page
   If you want to download the auto annotation pkl file based on celltypist(https://www.celltypist.org/) , select **Auto annotation** page
   If you want to query the code of GWAS enrichment with cell type, select **Traits enrichment** page
   If you want to download the resource, select **Data download** page, and download all the resource you want.
   <img width="334" height="445" alt="image" src="https://github.com/user-attachments/assets/2673328a-c41a-4c0f-9dd2-ba4b9edef58c" />

4. **Visualize & query**  
   - Use the **Gene Search** box to highlight expression of any gene.  
  <img width="747" height="488" alt="image" src="https://github.com/user-attachments/assets/ed3167bf-81aa-4f1c-9760-6449bc351a87" />
   - Export high-resolution figures or raw count matrices with one click.

5. **Datadownload locally**  
   - If you want to download the resource, select **Data download** page, and download all the resource you want. 

---

## 📚 Citation
If you use scGutDB or the IPGCA data in your research, please cite:

> 
BibTeX entry:
```bibtex
```

---

## 🤝 Contributing
We welcome bug reports, feature requests and pull requests!  
- **Issues**: [GitHub Issues](https://github.com/zjupatholab/scGutDB/issues)  
- **Discussions**: [GitHub Discussions](https://github.com/zjupatholab/scGutDB/discussions)  
- **Email**: 12417015@zju.edu.cn

---

## 📄 License
MIT © 2024 Zhejiang University - Pathology & Omics Lab.

---

## 🔗 Links
- **Live portal**: https://alphaindex.zju.edu.cn/scgut/
- **Panlab**: https://alphaindex.zju.edu.cn/Panlab/
