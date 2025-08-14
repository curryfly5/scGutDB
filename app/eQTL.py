import streamlit as st 
import streamlit_lottie
import sqlite3
import requests
# from script.featuremap import runUmap,runFeaturemap,runFeaturemap2
import pandas as pd
import re

def eQTL_page():
    st.title('eQTL/ieQTL query')
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Introduction')
    st.write("""
    <p style="font-size:18px; text-indent: 1em;">
    Here is gut eQTL database, we uploaded our eQTL/celltype-interaction eQTL(ct-ieQTL)/deconvolution celltype eQTL(decon-ct-eQTL) here.
    </p>
    """, unsafe_allow_html=True)
    # st.write('here is gut eQTL database, we uploaded our eQTL/celltype-interaction eQTL(ct-ieQTL)/deconvolution celltype eQTL(decon-ct-eQTL) here.')
    # with st.container():
    #     st.subheader('Atlas of gut scRNA-seq')
    # st.write('search the genes or celltype ...you intersted here!We provide celltype interaction eQTL and deconvolution celltype eQTL here!')
    st.write("""
    <p style="font-size:18px;">
    <ul>
        <li><b>QTL type</b>: choose eQTL/ieQTL/decon-ct-eQTL</li>
        <li><b>Tissue</b>: choose specific tissue for eQTL explore</li>
        <li><b>Gene ID/Name</b>: query specific gene to check if this gene is an eGene for any tissue/celltype.</li>
        <li><b>SNP ID</b>: query specific SNP to check if this SNP is an eQTL for any eGene.</li>
    </ul>
    </p>
    """, unsafe_allow_html=True)
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Query')
    # 选择物种、组织、谱系
    colum1,colum2=st.columns(2)
    with colum1:
        opt_QTL=st.selectbox('QTL type',('eQTL', 'ct-ieQTL', 'decon-eQTL'))
    with colum2:
        opt_tissue=st.selectbox('tissue',('ALL','Duodenum','Jejunum','IPEC_J2','Ileum', 'Colon','smallint','largeint'))
    # with colum3:
    #     opt_lineage=st.selectbox('celltype',('Enterocytes', 'B cells', 'Crypt', 'BEST4 enterocytes', 'Plasma cells','Goblet', 'Endothelial', 'Mesenchymal', 'Microfold'))
   # 猪数据（test 
    if opt_QTL=='eQTL':
        conn = sqlite3.connect('./sqlite/gut_eqtl.db')
        # 创建索引（如果尚未创建）
        cursor = conn.cursor()
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_gene_name ON [{opt_tissue}] (gene_name)")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_rsid ON [{opt_tissue}] (rsid)")
        conn.commit()
        # Streamlit 页面布局
        with st.container():
            col1, col2, col3 = st.columns([3, 3, 1])
            gene = col1.text_input('Gene ID/Name (optional)')
            snp_id = col2.text_input('SNP ID (optional)')
            with col3:
                button_container = st.container()
                with button_container:
                    st.write('<div style="margin-bottom: 28px;"></div>', unsafe_allow_html=True)
                    button = st.button('Submit', help='here')

        # 构建查询条件
        query = f"SELECT * FROM [{opt_tissue}] WHERE gene_name = ?"
        params = [gene]
        if snp_id:
            query += " AND rsid = ?"
            params.append(snp_id)

        # 按钮执行查询
        if button:
            try:
                # 执行查询并读取结果
                result_df = pd.read_sql(query, conn, params=params)
                if not result_df.empty:
                    st.dataframe(result_df)  # 显示查询结果
                else:
                    st.write("No results found for the given criteria.")
            except Exception as e:
                st.write(f"Error in querying the database: {e}")

        # 关闭连接
        conn.close()

    elif opt_QTL=='ct-ieQTL':
        conn = sqlite3.connect('./sqlite/gut_ieqtl.db')
        cursor = conn.cursor()
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_gene_name ON [{opt_tissue}] (gene_name)")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_rsid ON [{opt_tissue}] (rsid)")
        conn.commit()
        with st.container():
            col1, col2, col3, col4 = st.columns([1.5,1.5,1.5,1])
            gene = col1.text_input('Gene ID/Name (optional)')
            snp_id = col2.text_input('SNP ID (optional)')
            cell_type = col3.selectbox('Cell Type', ['Any','Enterocytes', 'B cells', 'Crypt', 'BEST4 enterocytes', 'Plasma cells','Goblet', 'Endothelial', 'Mesenchymal', 'Microfold'], index=0)
            with col4:
                button_container = st.container()
                with button_container:
                    # st.write('<div data-testid="stMarkdownContainer" class="css-184tjsw e16nr0p34"><p>Download</p></div>',unsafe_allow_html=True)
                    st.write('<div style="margin-bottom: 28px;"></div>', unsafe_allow_html=True)
                    button=st.button('Submit',help='here') 
                
                # # 用户输入 P 值和 FDR
                # col4, col5 = st.columns(2)
                # p_value = col4.text_input('P-value (immunQTL) (optional)')
                # fdr_value = col5.text_input('FDR (immunQTL) (optional)')
                
                # 构建查询条件
                query = f"SELECT * FROM [{opt_tissue}] WHERE gene_name = '{gene}'"
                if snp_id:
                    query += f" AND variant_id='{snp_id}'"
                if cell_type != 'Any':
                    query += f" AND CellType='{cell_type}'"
                # if p_value:
                #     query += f" AND P_value < {p_value}"
                # if fdr_value:
                #     query += f" AND FDR < {fdr_value}"

                # 按钮执行查询
            if button:
                try:
                    # 执行查询并读取结果
                    result_df = pd.read_sql(query, conn)
                    if not result_df.empty:
                        st.dataframe(result_df)  # 显示查询结果
                    else:
                        st.write("No results found for the given criteria.")
                except Exception as e:
                    st.write(f"Error in querying the database: {e}")
        conn.close()
        
    elif opt_QTL=='decon-ct-eQTL':
        conn = sqlite3.connect('./sqlite/gut_decon_eqtl.db')
        cursor = conn.cursor()
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_gene_name ON [{opt_tissue}] (gene_name)")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_rsid ON [{opt_tissue}] (rsid)")
        conn.commit()
        with st.container():
            col1, col2, col3, col4 = st.columns([1.5,1.5,1.5,1])
            gene = col1.text_input('Gene ID/Name (optional)')
            snp_id = col2.text_input('SNP ID (optional)')
            cell_type = col3.selectbox('Cell Type', ['Any','Enterocytes', 'B cells', 'Crypt', 'BEST4 enterocytes', 'Plasma cells','Goblet', 'Endothelial', 'Mesenchymal', 'Microfold'], index=0)
            with col4:
                button_container = st.container()
                with button_container:
                    # st.write('<div data-testid="stMarkdownContainer" class="css-184tjsw e16nr0p34"><p>Download</p></div>',unsafe_allow_html=True)
                    st.write('<div style="margin-bottom: 28px;"></div>', unsafe_allow_html=True)
                    button=st.button('Submit',help='here') 
                
                # # 用户输入 P 值和 FDR
                # col4, col5 = st.columns(2)
                # p_value = col4.text_input('P-value (immunQTL) (optional)')
                # fdr_value = col5.text_input('FDR (immunQTL) (optional)')
                
                # 构建查询条件
                query = f"SELECT * FROM [{opt_tissue}] WHERE gene_name = '{gene}'"
                if snp_id:
                    query += f" AND variant_id='{snp_id}'"
                if cell_type != 'Any':
                    query += f" AND CellType='{cell_type}'"
                # if p_value:
                #     query += f" AND P_value < {p_value}"
                # if fdr_value:
                #     query += f" AND FDR < {fdr_value}"

                # 按钮执行查询
            if button:
                try:
                    # 执行查询并读取结果
                    result_df = pd.read_sql(query, conn)
                    if not result_df.empty:
                        st.dataframe(result_df)  # 显示查询结果
                    else:
                        st.write("No results found for the given criteria.")
                except Exception as e:
                    st.write(f"Error in querying the database: {e}")
        conn.close()
if __name__ =='__main__':
    eQTL_page()