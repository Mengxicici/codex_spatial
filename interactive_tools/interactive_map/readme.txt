evn requirement is in requirements_streamlit.txt
In the current directory ：streamlit run streamlit_app.py --server.fileWatcherType=none
if default:
the parameter used will be saved in Map_param.yaml file
please save yaml before running map
the selected cell info will be saved in selected_1 (or name you assigned)
and can be import to dataframe exploration step again
so that can make subset analysis and can be compared with whole sample

bokeh server allows only one map,
 please close current map before making new map
