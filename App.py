import streamlit as st
# from app.batch import batch_page
from app.home import home_page
from app.marker import marker_page
from app.scGWAS import GWAS_page
from app.explore import explore_page
from app.auto import auto_page
from app.eQTL import eQTL_page
from app.download import download_page
from polish import Polish
# import streamlit_option_menu
from streamlit_option_menu import option_menu

# config setting
# st.set_page_config(page_title='scGUTdb',page_icon=':pig2:',layout="wide")
st.set_page_config(page_title='scGUTdb',page_icon=':pig2:')

# st.markdown("""
#     <style>
#     body {
#         zoom: 1.18; /* 修改此值调整缩放比例，1.0 为默认比例 */
#     }
#     </style>
#     """, unsafe_allow_html=True)

# navbar
with st.sidebar:
    st.image('./asset/image/LOGO.png')
    # st.title('scRNA-seq online web for gut')
    selected = option_menu(
    menu_title = "scGutDB",
    options = ["Home page","Marker query","Explore data","eQTLs query","Auto annotation","Traits enrichment","Data download"],
        # options = ["Home page","Marker query","Explore data","eQTLs query","Auto annotation","Batch correct","Traits enrichment"],
    icons = ["house-door","search","easel","cursor","bezier2","body-text","download"],
    # icons = ["house-door","search","easel","cursor","bezier2","bar-chart","body-text"],
    menu_icon = "list",
    default_index = 0,
        styles={
        "container": {"padding": "0!important", "background-color": "#F0F2F6"},
    }
)

# pages
if selected=="Home page":
    home_page()
if selected == "Marker query":
    marker_page()
if selected=="Explore data":
    explore_page()
if selected=="eQTLs query":
    eQTL_page()
if selected=="Auto annotation":
    auto_page()
# if selected=="Batch correct":
#     batch_page()
if selected=="Traits enrichment":
    GWAS_page()
if selected=="Data download":
    download_page()