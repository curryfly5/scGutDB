import streamlit as st
import streamlit.components.v1 as components
from streamlit_lottie import st_lottie
import pandas as pd
from polish import Polish
from PIL import Image
from json import load
import time

def render_metrics_cards():
    # 使用HTML和CSS创建React风格的卡片
    html_code = """
    <div style="display: flex; justify-content: space-between; gap: 10px;">
        <div style="
            flex: 1; 
            border: 1px solid #e0e0e0; 
            border-radius: 8px; 
            padding: 15px; 
            text-align: center; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            background-color: white;
        ">
            <h3 style="margin-bottom: 10px; color: #666;">Cells</h3>
            <div style="font-size: 24px; font-weight: bold; color: #333;">
                589,101
            </div>
        </div>
        
        <div style="
            flex: 1; 
            border: 1px solid #e0e0e0; 
            border-radius: 8px; 
            padding: 15px; 
            text-align: center; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            background-color: white;
        ">
            <h3 style="margin-bottom: 10px; color: #666;">Samples</h3>
            <div style="font-size: 24px; font-weight: bold; color: #333;">
                78
            </div>
        </div>
        
        <div style="
            flex: 1; 
            border: 1px solid #e0e0e0; 
            border-radius: 8px; 
            padding: 15px; 
            text-align: center; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            background-color: white;
        ">
            <h3 style="margin-bottom: 10px; color: #666;">Segments</h3>
            <div style="font-size: 24px; font-weight: bold; color: #333;">
                5
            </div>
        </div>
        
        <div style="
            flex: 1; 
            border: 1px solid #e0e0e0; 
            border-radius: 8px; 
            padding: 15px; 
            text-align: center; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            background-color: white;
        ">
            <h3 style="margin-bottom: 10px; color: #666;">Stages</h3>
            <div style="font-size: 24px; font-weight: bold; color: #333;">
                14
            </div>
        </div>
    </div>
    """
    
    # 在Streamlit中渲染HTML
    components.html(html_code, height=150)

# 高德地图的 API 密钥
amap_api_key = "	d6be0e8175dd5db3728b296fd8eb9dd5"  # 在这里替换为你的高德地图 API 密钥

# 定义高德地图 HTML 和 JavaScript 代码
map_html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="initial-scale=1.0, user-scalable=no, width=device-width">
    <style type="text/css">
    html, body, #container {{
        width: 100%;
        height: 100%;
        margin: 0px;
    }}
    </style>
    <script src="https://webapi.amap.com/maps?v=2.0&key={amap_api_key}"></script>
</head>
<body>
    <div id="container"></div>
    <script type="text/javascript">
    var map = new AMap.Map('container', {{
        zoom: 16,  // 缩放级别
        center: [120.08988, 30.29832]  // 初始中心点（经纬度），可替换为你想要的位置
    }});

    // 可选：添加标记
    var marker = new AMap.Marker({{
        position: new AMap.LngLat(120.08988, 30.29832)  // 标记位置
    }});
    map.add(marker);
    </script>
</body>
</html>
"""

def home_page():
    # st.markdown('<h1 style="font-size:42px;">Welcome to scGutDB</h1>', unsafe_allow_html=True)
    st.markdown("""
<style>
    .typing-text {
        font-family: 'Roboto', sans-serif;
        font-size: 42px;
        font-weight: 700;
        color: #2E4057;
        overflow: hidden;
        white-space: nowrap;
        margin: 0 0 5px 0;  /* Small bottom margin */
        letter-spacing: 1px;
        width: 0;
        animation: typing 1.5s steps(20, end) forwards;
    }
    
    @keyframes typing {
        to { width: 100%; }
    }
</style>

<h1 class="typing-text">Welcome to scGutDB</h1>
""", unsafe_allow_html=True)
    Polish.load_css('./style/style.css')     
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Introduction')
    with st.container():
        # lottie=Polish.load_lottie("https://assets10.lottiefiles.com/private_files/lf30_tzlxqk5k.json")#比较慢
        lottie=Polish.load_lottie("./asset/json/DNA_lottie.json",local=True)
        l_colum,r_colum=st.columns(2)
        with l_colum:
            # st.write('ScGutDB是一个在线的关于肠道单细胞数据的探索平台，为了充分利用单细胞数据资源，平台提供了以下功能人、猪、鼠肠道相关研究的单细胞数据，并提供了全面的猪小肠、大肠的单细胞细胞类型图谱，还辅助以单细胞GWAS映射、批次效应校正、marker基因查询、多层级自动注释等功能。')
            content1="""
            <p style="font-size:18px; text-indent: 1em;"><b>ScGutDB</b> is an online exploration platform for fully utilizing the resources of intestinal single-cell data which enables to visualize and explore some founding in the paper. It provides our comprehensive atlas IPGCA, as well as some human, pig, and mouse intestines scRNA-seq datasets. The platform also facilitates the sharing analysis code and offers functionalities for query and visualization.</p>
            """
            st.markdown(content1, unsafe_allow_html=True)
        with r_colum:
            st_lottie(lottie,key='coding',height=200)
        st.write('<p style="font-size:18px;"><b>IPGCA consists:</b></p>', unsafe_allow_html=True)
        col1, col2, col3 ,col4= st.columns(4)
        with col1:
            with st.container(border=True,height=120):
                st.metric("Cells", "589,101", "")
        with col2:
            with st.container(border=True,height=120):
                st.metric("Samples", "78", "")
        with col3:
            with st.container(border=True,height=120):
                st.metric("Segments", "5", "")
        with col4:
            with st.container(border=True,height=120):
                st.metric("Stages", "14", "")
        # render_metrics_cards()
        # content2="""
        # 1.Home page: provides an overview of the website and showcases the single-cell atlases of IPGCA; \n
        # 2.Marker query: query of the differential expressed gene (DEG) analysis result of IPGCA for each cell type and visualization tools; \n
        # 3.Explore data: interactively exploration of IPGCA and Human gut cell atlas, offering visualization for any gene of interest along with all associated covariates and dynamic gene expression patterns. \n 
        # 4.eQTL query: the interrogation of a comprehensive repository of eQTL data, encompassing eQTLs, ct-ieQTLs, and decon-ct-eQTLs. Users can perform searches by employing identifiers such as SNP ID, Gene ID, Tissue and cell type. \n 
        # 5.Auto annotation:a large collection of intestinal single-cell data to train automatic annotation models (celltypist.pkl files) for various cell lineages; \n 
        # 6.Batch correct: The batch effect correction results from IPGCA are provided, along with the processing code used for this analysis. \n 
        # 7.Traits enrichment: sharing the process of enrichment between intestinal cell types and relevant phenotypic data using software such as MAGMA and scDRS and performs heritability enrichment analysis with LDSC-seg. \n
        # """
        # st.markdown(content2, unsafe_allow_html=True)
    st.write('---')
    
    with st.container():
        st.subheader('Atlas of gut scRNA-seq')
        # content2="""
        # <main style="font-size:20px;">  Here, we present an integrated single-cell atlas of the gut, encompassing both porcine and human datasets. Users can select species of interest and delve into specific cell lineages to explore the diversity of cellular types, facilitating the investigation of cellular heterogeneity and function within the intestinal tract.</main>
        # """
        st.write('<p style="font-size:18px; text-indent: 1em;">Here, we present an integrated single-cell atlas of the gut, encompassing both <b>porcine</b> and <b>human datasets</b>. Users can select species of interest and delve into specific cell lineages to explore the diversity of cellular types, facilitating the investigation of cellular heterogeneity and function within the intestinal tract.</p>', unsafe_allow_html=True)
        # st.markdown(content2, unsafe_allow_html=True)
    colum1,colum2=st.columns(2)
    with colum1:
        opt_specice = st.selectbox('species',('pig','human'))
    if opt_specice == 'pig':
        with colum2:
            opt_Lineage = st.selectbox('Lineage',('All','T_ILC', 'B','Myeloid','Epithelial','Endothelial','Mesenchymal'))
        images=f'./asset/image/{opt_specice}_{opt_Lineage}.png'
        st.image(images)
    elif opt_specice == 'human':
        with colum2:
            opt_Lineage = st.selectbox('Lineage',('All','T_ILC', 'B','Myeloid','Epithelial','Endothelial','Mesenchymal','Neuron'))
        images=f'./asset/image/{opt_specice}_{opt_Lineage}.png'
        st.image(images)
    # elif opt_specice == 'mouse':
    #     return None
    st.write('---')
    st.subheader('Connections')
    with st.container(border=True):
        colum3,colum4=st.columns(2)
        with colum3:
            # 将 HTML 和 JavaScript 嵌入 Streamlit
            components.html(map_html, height=150)
        with colum4:
            content1="""
            <p style="line-height:1.5;font-size:16px;">
            <b style="font-size:20px;">Pengfei Yu</b><br>
            PhD Candidate, Yuchun Pan Lab, College of Animal Sciences, Zhejiang University, Hangzhou, Zhejiang 310058, China<br>
            Email:  <a href="mailto:pengfei_yu@zju.edu.cn">pengfei_yu@zju.edu.cn</a>
            </p>
            """
            st.markdown(content1, unsafe_allow_html=True)
    # content2="""
    # <p style="line-height:1.5;font-size:16px;">
    # <b>Qinqin Xie</b><br>
    # PhD Candidate, Yuchun Pan Lab, College of Animal Sciences, Zhejiang University, Hangzhou, Zhejiang 310058, China<br>
    # Email:  <a href="mailto:qinqin.xie@zju.edu.cn">qinqin.xie@zju.edu.cn</a>
    # </p>
    # """
    # st.markdown(content2, unsafe_allow_html=True)



if __name__ =='__main__':
    home_page()