import streamlit as st 
import streamlit_lottie
import requests
import pandas as pd
import base64
from streamlit.components.v1 import iframe

def auto_page():
    st.title('Auto annotation')
    st.markdown('<hr style="margin-top: 10px; margin-bottom: 20px;">', unsafe_allow_html=True)
    st.subheader('Introduction')
    st.write("""
    <p style="font-size:18px; text-indent: 1em;">
    Here is auto annotation part, we uploaded our pre-trained cell atlas pkl modules here, you can download the pkl and predicti with software <a href="https://www.celltypist.org/" target="_blank">celltypist</a>.
    </p>
    """, unsafe_allow_html=True)
    # 添加一个下载按钮的方法
    def button(file,name):
        with open(file, 'rb') as f:
            data = f.read()
        st.download_button('Download', data,file_name=name)
    # pkl简介
    data = [
    ['IPGCA_lineage.pkl','IPGCA reference at cell lineage level'],
    ['IPGCA_celltype.pkl','IPGCA reference at cell type level'],
    ['human_linage.pkl', 'Human gut reference'],
    ['B_lineage.pkl','IPGCA reference for B cell lineage'],
    ['Endothelial_lineage.pkl','IPGCA reference for Endothelial cell lineage'],
    ['Epithelial_lineage','IPGCA reference for Epithelial cell lineage'],
    ['Mesenchymal_lineage.pkl','IPGCA reference for Mesenchymal cell lineage'],
    ['Myeloid_lineage.pkl','IPGCA reference for Myeloid cell lineage'],
    ['Neuron_Lineage.pkl','IPGCA reference for Neuron cell Lineage'],
    ['T_ILC_lineage.pkl','IPGCA reference for T/ILC cell lineage'],
    ]
    df = pd.DataFrame(data, columns=['File Name', 'Description'])
    st.table(df)
    #添加一些选项并绑定下载按钮
    colum1,colum3,colum4=st.columns(3)
    with colum1:
        opt_specice=st.selectbox('species',('pig','human'))
    # with colum2:
    #     opt_tissue=st.selectbox('tissue',('small intestine','colon'))
    #     if opt_tissue=='small intestine':
    #         opt_tissue='SI'
    with colum3:
        opt_lineage=st.selectbox('Lineage',("IPGCA_Lineages",'B lineage', 'Endothelial lineage', 'Epithelial lineage', 'Mesenchymal lineage', 'Myeloid lineage', 'Neuronal lineage', 'T/ILC lineage'))
    with colum4:
      button_container = st.container()
      with button_container:
        # st.write('<div data-testid="stMarkdownContainer" class="css-184tjsw e16nr0p34"><p>Download</p></div>',unsafe_allow_html=True)
        st.write('<div style="margin-bottom: 28px;"></div>', unsafe_allow_html=True)
        opt_lineage=opt_lineage.replace(' ','_')
        opt_lineage=opt_lineage.replace('/','_')
        button(f'./asset/pkl/{opt_lineage}.pkl', f'{opt_lineage}.pkl')
      # Adjust the position of the container using CSS
      
    # button('asset/pkl/human_linage.pkl', 'human_linage.pkl')
    
    st.write('''
             ''',unsafe_allow_html=True)
    st.markdown('<hr style="margin-top: 10px; margin-bottom: 10px;">', unsafe_allow_html=True)
    with open('./asset/md/query.md', 'r', encoding='utf-8') as file:
        markdown_content = file.read()

    # 在Streamlit中显示Markdown内容
    st.markdown(markdown_content)
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    # with open("asset/md/celltypist.md", "r",encoding="utf-8") as f:
    #     markdown_text = f.read()
    # st.subheader('celltypist简介')
    # with st.container():
    #     colum1,colum2=st.columns(2)
    #     with colum1:
    #         st.markdown(markdown_text)
    #     with colum2:
    #         st.image('asset/image/celltypist.png')
        
    # def get_binary_file_downloader_html(bin_file, file_label='File'):
    #     with open(bin_file, 'rb') as f:
    #         data = f.read()
        # bin_str = base64.b64encode(data).decode()
        # href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{file_label}.pkl">Download {file_label}</a>'
        # return href
        

if __name__ =='__main__':    
    auto_page()
