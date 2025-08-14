import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sb
import streamlit as st
import logging
import time
import subprocess

@st.cache_resource
def read_h5ad(url):
    adata=sc.read_h5ad(url, backed='r')
    return adata

def runUmap(adata,color):
    # adata=read_h5ad(f'./asset/data/JH_CellLineage_{tissue}_{period}.h5ad')
    plot, ax = plt.subplots()
    with plt.rc_context({'figure.figsize': [10, 10]}):
        kwargs={'frameon': False }
        sc.pl.umap(adata,color=color,ax=ax,size=10,cmap='Reds',legend_loc='right margin' ,**kwargs)
    st.pyplot(plot)
    # del adata

#labdata所用featuremap   
def runFeaturemap(adata,gene,ncols=2):
    gene=gene.split(',')
    # adata=read_h5ad(f'./asset/data/JH_CellLineage_{tissue}_{period}.h5ad')
    try:
        with st.spinner('Plotting the graph......'):
            time.sleep(3)
            # plot, ax = plt.subplots()
            with plt.rc_context({'figure.figsize': [5, 5]}):
                kwargs={'frameon': False }
                plot=sc.pl.umap(adata,color=gene,return_fig=True,ncols=ncols,size=10,cmap='Reds',**kwargs)
            # st.pyplot(plot)
    except:
        st.write('one of these genes is not in the vars of anndata')
    return plot

# #marker所用featuremap   
# def runFeaturemap2(tissue,period,gene):
#     gene=gene.split(',')
#     color=['Lineage']
#     color=color+gene
#     adata=read_h5ad(f'./asset/data/JH_CellLineage_{tissue}_{period}.h5ad')
#     try:
#         with st.spinner('Plotting the graph......'):
#             time.sleep(3)
#             plot=sc.pl.umap(adata,color=color,return_fig=True,wspace=0.3,ncols=2)
#             # plot2=sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=4,min_logfoldchange=4)
#             st.pyplot(plot)
#             # st.pyplot(plot2)
#     except:
#         st.write('one of these genes is not in the vars of anndata')
    # del adata
# def runFeaturemap(tissue,period,gene,cluster=False):
    # gene=gene.split(',')
    # adata=sc.read_h5ad(f'./asset/data/JH_CellLineage_CE_0d.h5ad')
    
    # with st.spinner('Plotting the graph......'):
    #     time.sleep(3)
    #     plot=sc.pl.umap(adata,color=gene,return_fig=True,ncols=2)
    #     st.pyplot(plot)
    #     if cluster:
    #         plot2=sc.pl.umap(adata,color='Lineage',return_fig=True)
    #         st.pyplot(plot2)
    # del adata
    
#violinplot (Gene跨组织表达)
def runviolinplot(alldata,obs,gene):
    gene=gene.split(',')
    # alldata=read_h5ad('/disk212/yupf/database/web/xqqweb/asset/data/adata_CellLineage_rank_genes_groups.h5ad')  
    try:
        with st.spinner('Plotting the graph......'):
            time.sleep(3)
            plot2, ax = plt.subplots()
            sc.pl.violin(alldata, keys=gene, groupby=obs, stripplot=False, rotation=90,ax=ax)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            # plot=sc.pl.umap(adata,color=color,return_fig=True,wspace=0.3,ncols=2)
            st.pyplot(plot2, bbox_inches='tight')
    except:
        st.write('one of these genes is not in the vars of anndata')
    # del adata
    return plot2
#Comparison of marker genes using split violin plots
    
#dotplot 
def rundotplot(alldata,obs,gene):
    gene=gene.split(',')
    # alldata=read_h5ad('/disk212/yupf/database/web/xqqweb/asset/data/adata_CellLineage_rank_genes_groups.h5ad')  
    try:
        with st.spinner('Plotting the graph......'):
            time.sleep(3)
            plot3, ax = plt.subplots()
            sc.pl.dotplot(alldata, var_names=gene, ax=ax,groupby=obs)
            # plot=sc.pl.umap(adata,color=color,return_fig=True,wspace=0.3,ncols=2)
            st.pyplot(plot3, bbox_inches='tight')
    except:
        st.write('one of these genes is not in the vars of anndata')
    # del adata
    return plot3