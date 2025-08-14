import streamlit as st 
import streamlit_lottie as sl
from streamlit_lottie import st_lottie
import requests
import json


class Polish:
    def load_lottie(url,local=False):
        if local==False:
            r =requests.get(url)
            if r.status_code !=200:
                return None
            return r.json()
        else:
            with open(url) as J:
             return json.load(J)
        
    
    def load_css(css_file):
        with open(css_file) as f:
            st.markdown(f'<style>{f.read()}</style>',unsafe_allow_html=True)