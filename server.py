import os
os.chdir('/disk212/yupf/database/web/streamlit/')
cmd='streamlit run App.py --server.port 8501 --browser.serverAddress alphaindex.zju.edu.cn'
os.system(cmd)
