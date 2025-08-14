import streamlit as st 
import streamlit_lottie
import requests


def batch_page():
    st.title('Batch effect evaluation')
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    with open('/disk212/yupf/database/scRNA-seq/benchmark/CE69/plot.md', 'r', encoding='utf-8') as file:
        markdown_content = file.read()

    # 在Streamlit中显示Markdown内容
    st.markdown(markdown_content)
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
if __name__ =='__main__':
    batch_page()