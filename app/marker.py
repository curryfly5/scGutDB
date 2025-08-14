import streamlit as st 
import streamlit_lottie
import requests
from script.featuremap import runUmap,runFeaturemap,read_h5ad
import pandas as pd
import re

def marker_page():
    st.title('Marker gene query')
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Introduction')
    content1='<p style="font-size:18px; text-indent: 1em;">Here is marker gene database extracted from <b>IPGCA</b>, you can search by specific cell lineage and specific gene for both human and pig intestine.</p>'
    st.markdown(content1, unsafe_allow_html=True)
    # with st.container():
    #     st.subheader('Atlas of gut scRNA-seq')
    st.write("""
    <ul>
        <li><b>By Lineage</b>: exhibit DEGs (Differential expressed genes) among cell types of this Lineage, default is ALL.</li>
        <li><b>By Gene</b>: Search for the genes that interest you here. The terms of differentially expressed genes (DEGs) containing this gene will be retained, and the corresponding feature map will be plotted below.</li>
    </ul>
    """, unsafe_allow_html=True)
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Query')
    # 选择物种、组织、谱系
    colum1,colum2=st.columns(2)
    with colum1:
        opt_specice=st.selectbox('species',('pig','human'))
    with colum2:
        opt_lineage=st.selectbox('Lineage',('All','B lineage', 'Endothelial lineage', 'Epithelial lineage', 'Mesenchymal lineage', 'Myeloid lineage', 'Neuronal lineage', 'T/ILC lineage'))
   
   # 猪数据（test 
    if opt_specice=='pig':
        with st.container():
            colum3,colum4=st.columns([7,1])
            with colum3:
                search_term=st.text_input('input the gene you interested')
            with colum4:
                button_container = st.container()
                with button_container:
                    # st.write('<div data-testid="stMarkdownContainer" class="css-184tjsw e16nr0p34"><p>Download</p></div>',unsafe_allow_html=True)
                    st.write('<div style="margin-bottom: 28px;"></div>', unsafe_allow_html=True)
                    button=st.button('Submit') 
        df=pd.read_csv('./asset/table/celltype_DEG.csv',index_col=0)
        #实现选择谱系和对应dataframe
        if opt_lineage=='All':
            st.checkbox("Use container full width", value=True, key="use_container_width")
            st.dataframe(df,use_container_width=st.session_state.use_container_width) 
        else:  
            st.checkbox("Use container full width", value=True, key="use_container_width")
            st.dataframe(df[df['lineage'].str.contains(opt_lineage)].copy(),use_container_width=st.session_state.use_container_width)    
        if button:
            st.empty()
            ser=re.sub(r'[,]{1,50}','|',search_term)
            result=df[df['names'].str.contains(ser)].copy()
            st.write('Query result:')
            st.dataframe(result,use_container_width=st.session_state.use_container_width)
            adata=read_h5ad(f'./asset/data/{opt_specice}_01.h5ad')
            st.write('Featuremap visualization:')
            search_term='CellType,'+search_term
            plot=runFeaturemap(adata,search_term,ncols=1)
            st.pyplot(plot)
    else:
        return None
        # opt_1 = st.selectbox(
        # 'tissue',
        # ('CE', 'CO', 'JE','IL','DU'))
        # opt_2 = st.selectbox(
        # 'time',
        # ('0d', '60d','90d','180d','240d'))
        # runUmap(opt_1,opt_2)
        # with st.container():
        #     st.subheader('Feature map of the sample')
        #     st.write('choose Gene you intersted here!')
        #     txt_3= st.text_input('Text your gene to analyze')#('Text your gene to analyze', '''CD79''')
        #     if txt_3:
        #         runFeaturemap(opt_1,opt_2,txt_3)

if __name__ =='__main__':    
    marker_page()