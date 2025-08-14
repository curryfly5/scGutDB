import streamlit as st 
import streamlit_lottie
import requests
from script.featuremap import read_h5ad,runUmap,runFeaturemap,runviolinplot,rundotplot
import os
import matplotlib.pyplot as plt
from PIL import Image

def explore_page():
    
    st.markdown('<h1 style="font-size:42px;">Explore IPGCA</h1>', unsafe_allow_html=True)
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    # st.markdown(
    #     '<div style="font-size:32px; margin-bottom: 5px;"><b>Atlas</b></div>', 
    #     unsafe_allow_html=True
    # )
    # st.markdown('<h1 style="font-size:36px;">Introduction</h1>', unsafe_allow_html=True)
    st.subheader('Introduction')
    st.write("""
    <p style="font-size:18px; text-indent: 1em;">
    You can explore both human and pig intestine cell atlas here! Human data is from <a href="https://www.gutcellatlas.org" target="_blank">gutcellatlas.org</a>, and the pig data is from IPGCA(Pig Gut Cell Atlas), a half-million cell (589,101) atlas of the porcine intestine encompassing 5 intestinal segments, 14 developmental stages, and 7 breeds. 
    </p>
    """, unsafe_allow_html=True)
    st.write('---')
    with st.container():
        st.subheader('Visualization')
        st.write("""
        <p style="font-size:18px; text-indent: 1em;">
        Choose the datasets and different meta information you interested here, like segment,gender,stage...
        </p>
        """, unsafe_allow_html=True)
        # st.write('Choose the datasets and different meta information you interested here, like segment,gender,stage...')
    st.write('''
    <p style="line-height:1.5;font-size:14px;margin-bottom: 3px">
    Examples: pig_01.h5ad + CellTyped
    </p>
    ''', unsafe_allow_html=True)
    
    # 加载数据（假设 read_h5ad 已经使用 @st.cache_data）
    
    
## 可视化umap和featuremap      
    columA, columB = st.columns(2)
    files = os.listdir('/disk212/yupf/database/web/streamlit/asset/data/')
    target_element1 = 'pig_01.h5ad'
    new_files = [target_element1] + [item for item in files if item != target_element1]
    with columA:
        opt_1 = st.selectbox('Files',  tuple(new_files))
        adata = read_h5ad(f'./asset/data/{opt_1}')
    with columB:
        obs=adata.obs.columns
        # 将默认移到第一位的元素
        target_element2 = 'CellType'
        new_obs = [target_element2] + [item for item in obs if item != target_element2]
        opt_2 = st.selectbox('Dimensions', tuple(new_obs))
    # run UMAP
    if opt_1 == target_element1 and opt_2 == target_element2:
            st.image('./asset/image/static_CellType.png')
    else:
        umapplot=runFeaturemap(adata, opt_2)
        st.pyplot(umapplot)
    st.write('---')
    
    # 特征映射功能
    with st.container():
        st.subheader('Gene expression')
        if 'gene_input' not in st.session_state:
            st.session_state['gene_input'] = 'CD79B,CD3E'
            st.session_state['featuremap_initialized'] = False
        st.write('''
        <p style="line-height:1.5;font-size:14px;margin-bottom: 3px">
        Examples: CD79B,CD3E
        </p>
        ''', unsafe_allow_html=True)
        txt_3 = st.text_input('', placeholder='CD79B,CD3E', value=st.session_state['gene_input'], key="gene_input_input")

        if not st.session_state['featuremap_initialized']:
            static_feature=Image.open('./asset/image/static_umap_CD79BCD3E.png')
            # 将PNG转换为Matplotlib 图形
            plot, ax = plt.subplots()
            ax.imshow(static_feature)
            ax.axis('off')
            st.session_state['featuremap_plot'] = plot
            st.session_state['featuremap_initialized'] = True
        
        elif st.session_state['gene_input'] != txt_3:
            st.session_state['gene_input'] = txt_3
            plot = runFeaturemap(adata, txt_3)
            st.session_state['featuremap_plot'] = plot

        if 'featuremap_plot' in st.session_state:
            st.pyplot(st.session_state['featuremap_plot'])
    st.write('---')
## 可视化violin plot 
    with st.container():
        st.subheader('Dynamic expression')
        st.write('''
        <p style="line-height:1.5;font-size:14px;margin-bottom: 3px">
        Examples: Violin + CD79B + CellType
        </p>
        ''', unsafe_allow_html=True)
        cat_col,colum5, colum6 = st.columns(3)
        with cat_col:
            plt_cat = st.selectbox('', ( 'Violin Plot','Dot Plot'))
        with colum5:
            txt_4 = st.text_input('', placeholder='CD79B', value='CD79B', key="violin_gene_input")
        with colum6:
            opt_4 = st.selectbox('', tuple(new_obs))
        if plt_cat=='Violin Plot':
            if txt_4 == 'CD79B' and opt_4 == target_element2:
                st.image('./asset/image/static_violion_CD79B.png')
            else:
                runviolinplot(adata, opt_4, txt_4)
        if plt_cat=='Dot Plot':
            rundotplot(adata, opt_4, txt_4)

if __name__ =='__main__':
    explore_page()