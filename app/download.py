import streamlit as st
import pandas as pd
import os
import base64
import io
import scanpy as sc

def download_page():
    st.title('Data Download')
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Introduction')
    content1 = '''
    <p style="font-size:18px; text-indent: 1em;">Welcome to the data download center of IPGCA (Integrated Pig Gut Cell Atlas). This page provides comprehensive resources for researchers interested in intestinal single-cell genomics. Our database includes:</p>
    <ul style="font-size:18px;">
        <li><b>Differential Expression Gene (DEG) Tables</b>: Cell type-specific and lineage-specific DEGs from our integrated analysis</li>
        <li><b>Single-cell Data</b>: Complete h5ad format datasets of pig and human gut cell atlases</li>
        <li><b>eQTL Resources</b>: Three types of expression quantitative trait loci data (bulk eQTL, cell type-interaction eQTL, and deconvolution-based cell type eQTL)</li>
    </ul>
    <p style="font-size:18px; text-indent: 1em;">All datasets are freely available for academic research. For commercial use, please contact us through the information provided at the bottom of this page.</p>
    '''
    st.markdown(content1, unsafe_allow_html=True)
    
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('DEG Tables')
    
    # Create download function
    def get_table_download_link(df, filename, text):
        """Generate table download link"""
        csv = df.to_csv(index=True)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="{filename}" style="display: inline-block; padding: 0.5em 1em; text-decoration: none; border: 1px solid #ccc; border-radius: 4px; margin: 5px 0; color: black;">{text}</a>'
        return href
    

    
    # Table data download section
    st.write("The following table data files are available for download:")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Load celltype_DEG.csv
        celltype_deg = pd.read_csv('./asset/table/celltype_DEG.csv', index_col=0)
        st.markdown(get_table_download_link(celltype_deg, 'celltype_DEG.csv', 'Download Cell Type DEGs'), unsafe_allow_html=True)
        
        # Show data preview
        with st.expander("Preview Cell Type DEGs"):
            st.dataframe(celltype_deg.head(10))
    
    with col2:
        # Load lineage_DEG.csv
        lineage_deg = pd.read_csv('./asset/table/lineage_DEG.csv', index_col=0)
        st.markdown(get_table_download_link(lineage_deg, 'lineage_DEG.csv', 'Download Lineage DEGs'), unsafe_allow_html=True)
        
        # Show data preview
        with st.expander("Preview Lineage DEGs"):
            st.dataframe(lineage_deg.head(10))
    
    # End of DEG Tables section
    
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Single-cell Data (h5ad format)')

    st.write("The following h5ad format single-cell data files are available for download:")
    st.write("Note: h5ad files are large and may take some time to download.")

    col5, col6 = st.columns(2)

    with col5:
        st.write("Pig Gut Single-cell Data")
        st.markdown(f'<a href="https://example.com/pig_01.h5ad" download style="display: inline-block; padding: 0.5em 1em; text-decoration: none; border: 1px solid #ccc; border-radius: 4px; margin: 5px 0; color: black;">Download Pig Gut Data</a>', unsafe_allow_html=True)

    with col6:
        st.write("Human Gut Single-cell Data")
        st.markdown(f'<a href="https://example.com/human_01.h5ad" download style="display: inline-block; padding: 0.5em 1em; text-decoration: none; border: 1px solid #ccc; border-radius: 4px; margin: 5px 0; color: black;">Download Human Gut Data</a>', unsafe_allow_html=True)
    
    # eQTL Data Section
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('eQTL Data Resources')
    
    st.write("The following eQTL datasets are available for download:")
    st.write("Note: These datasets provide valuable insights into gene expression regulation in intestinal tissues.")
    
    # Use 3 columns in one row for eQTL data
    col7, col8, col9 = st.columns(3)
    
    with col7:
        # Bulk eQTL data
        bulk_eqtl_path = './asset/table/eQTL/bulk_eQTL.csv'  # Placeholder path
        try:
            # For demonstration - in reality would load the actual file
            bulk_eqtl_sample = pd.DataFrame({
                'SNP_ID': ['rs123456', 'rs234567', 'rs345678'],
                'Gene': ['GENE1', 'GENE2', 'GENE3'],
                'P_value': [1.2e-8, 5.6e-7, 3.4e-6],
                'Effect': [0.45, -0.32, 0.28]
            })
            st.markdown(get_table_download_link(bulk_eqtl_sample, 'bulk_eQTL.csv', 'Download Bulk eQTL Data'), unsafe_allow_html=True)
            
            # Show data preview
            with st.expander("Preview Bulk eQTL Data"):
                st.dataframe(bulk_eqtl_sample)
        except Exception as e:
            st.write(f"Bulk eQTL data is being prepared and will be available soon.")
    
    with col8:
        # Cell type-interaction eQTL (ct-ieQTL) data
        ct_ieqtl_path = './asset/table/eQTL/ct_ieQTL.csv'  # Placeholder path
        try:
            # For demonstration - in reality would load the actual file
            ct_ieqtl_sample = pd.DataFrame({
                'SNP_ID': ['rs123456', 'rs234567', 'rs345678'],
                'Gene': ['GENE1', 'GENE2', 'GENE3'],
                'Cell_Type': ['B cells', 'T cells', 'Enterocytes'],
                'P_value': [1.2e-8, 5.6e-7, 3.4e-6],
                'Interaction_Effect': [0.45, -0.32, 0.28]
            })
            st.markdown(get_table_download_link(ct_ieqtl_sample, 'ct_ieQTL.csv', 'Download ct-ieQTL Data'), unsafe_allow_html=True)
            
            # Show data preview
            with st.expander("Preview ct-ieQTL Data"):
                st.dataframe(ct_ieqtl_sample)
        except Exception as e:
            st.write(f"Cell type-interaction eQTL data is being prepared and will be available soon.")
    
    with col9:
        # Deconvolution-based cell type eQTL (decon-eQTL) data
        decon_eqtl_path = './asset/table/eQTL/decon_eQTL.csv'  # Placeholder path
        try:
            # For demonstration - in reality would load the actual file
            decon_eqtl_sample = pd.DataFrame({
                'SNP_ID': ['rs123456', 'rs234567', 'rs345678'],
                'Gene': ['GENE1', 'GENE2', 'GENE3'],
                'Cell_Type': ['B cells', 'T cells', 'Enterocytes'],
                'P_value': [1.2e-8, 5.6e-7, 3.4e-6],
                'Effect': [0.45, -0.32, 0.28]
            })
            st.markdown(get_table_download_link(decon_eqtl_sample, 'decon_eQTL.csv', 'Download decon-eQTL Data'), unsafe_allow_html=True)
            
            # Show data preview
            with st.expander("Preview decon-eQTL Data"):
                st.dataframe(decon_eqtl_sample)
        except Exception as e:
            st.write(f"Deconvolution-based cell type eQTL data is being prepared and will be available soon.")
    
    # Contact section
    st.markdown('<hr style="margin-top: 20px; margin-bottom: 5px;">', unsafe_allow_html=True)
    st.subheader('Contact')
    st.write("For questions, additional data requests, or to report issues:")
    st.write("📧 Email: pengfei_yu@zju.edu.cn")
    # st.write("🔗 GitHub: https://github.com/example/IPGCA")
    st.write("🏢 Institution: College of Animal Sciences, Zhejiang University, Hangzhou, China")

if __name__ == '__main__':
    download_page()