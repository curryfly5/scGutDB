import streamlit as st 
import streamlit_lottie
import requests
import scanpy as sc


def GWAS_page():
    st.title('GWAS enrichment')
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    # 读取Markdown文件内容
    st.subheader('Tutorial')
    with open('./asset/md/markdown_magma_scDRS.md', 'r', encoding='utf-8') as file:
        markdown_content = file.read()

    # 在Streamlit中显示Markdown内容
    st.markdown(markdown_content)
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    file = st.file_uploader("Upload file", type="h5ad")

    if file is not None:
        adata = sc.read_h5ad(file)
        plot=sc.pl.umap(adata,color='majority_voting',return_fig=True)
        st.pyplot(plot)
        
if __name__ =='__main__':
    GWAS_page()
    
