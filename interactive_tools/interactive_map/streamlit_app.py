import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import yaml
import bokeh
import streamlit as st
import os
import subprocess
import webbrowser
from bokeh import events
from bokeh.io import output_file,output_notebook, show
from bokeh.models import ColumnDataSource, LinearColorMapper, Whisker,Selection,Select, CustomJS,Button, CustomJS, Button,Div,Legend
from bokeh.plotting import figure, show
from bokeh.sampledata.penguins import data
from bokeh.transform import factor_cmap,jitter
from bokeh.models.tools import BoxSelectTool
from bokeh.palettes import Spectral6, RdBu,Purples256,Spectral3,Spectral4,Plasma256,Viridis256,Turbo256,Inferno256 #TolRainbow3,
from bokeh.models.glyphs import VBar,Scatter,Circle
from bokeh.layouts import row, column,gridplot
from bokeh.server.server import Server
from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler

from density_tools import *

st.write("Feb 3 2023")
st.header('test streamlit')
if st.button("user guide"):
    """
    1. Please uncheck unnecessary visualization to speed up the calculation
    2. Interested cell subsets found in Map Exploration can be reloaded to Dataframe Exploration to conduct spatial analysis and visualization of subsets.
    3. There are some interactive plots, you can try to click the legend
    4. May need 10-20 secs to load big interlocked map from bokeh server
    """
if st.button("cell sociology"):
    """
    # welcome to spactial cell sociology workflow

    ## Spatial cell sociology hypothesis​:

    1. Cell-cell communication can be partially deciphered​,with the position and density information of cell and sub-cellular proteomic/scRna-seq information 

    2. Cell-cell communication have a common pattern within one population while different among various cell populations​

    3. The conversation can be affected by TME(neighborhood/density/regions), and the communication have impact on TME
    """

# """
# For BioHPC:
# path1=r"/Volumes/archive/bioinformatics/Jamieson_lab/shared/CODEX/cytomap_share/leiden2p5_spatial_matrix/density/cell_dens_r40/CytoMAP_Sample_P7T1_leiden_25_testing.csv"
# path2=r"/Volumes/archive/bioinformatics/Jamieson_lab/shared/CODEX/cytomap_share/leiden2p5_spatial_matrix/density/cell_dens_r40/CytoMAP_Sample_P7T2_leiden_25_testing.csv"



# #for DELL:
# path1=r"Z:\Shared\cytomap_share\leiden2p5_spatial_matrix\density\cell_dens_r40\CytoMAP_Sample_P7T1_leiden_25_testing.csv"
# path2=r"Z:\Shared\cytomap_share\leiden2p5_spatial_matrix\density\cell_dens_r40\CytoMAP_Sample_P7T2_leiden_25_testing.csv"
# """

    

select_event = st.sidebar.selectbox(
    "which step to go?",
    ("Dataframe exploration", "Map exploration", "Calculation tools")
)


st.write("Data loading")
if os.path.exists('/endosome/'):
    print("loading from bioHPC")
    
    path1=r"/endosome/archive/bioinformatics/Jamieson_lab/shared/CODEX/cytomap_share/leiden2p5_spatial_matrix/density/cell_dens_r40/CytoMAP_Sample_P7T1_leiden_25_testing.csv"
    path2=r"/endosome/archive/bioinformatics/Jamieson_lab/shared/CODEX/cytomap_share/leiden2p5_spatial_matrix/density/cell_dens_r40/CytoMAP_Sample_P7T2_leiden_25_testing.csv"
    
# Verify that the network drive is accessible
elif os.path.exists('Z:\\'):
    st.write("loading from DELL")
    path1=r"Z:\Shared\cytomap_share\leiden2p5_spatial_matrix\density\cell_dens_r40\CytoMAP_Sample_P7T1_leiden_25_testing.csv"
    path2=r"Z:\Shared\cytomap_share\leiden2p5_spatial_matrix\density\cell_dens_r40\CytoMAP_Sample_P7T2_leiden_25_testing.csv"

else:
    st.write("loading from MAC")
    path1=r"/Users/mengxi/Documents/utsw/codex/codex_spatial/cytomap/r60_reg20/csv/cell_density/cell_density_r40/CytoMAP_Sample_P7T1_leiden_25_testing.csv"
    path2=r"/Users/mengxi/Documents/utsw/codex/codex_spatial/cytomap/r60_reg20/csv/cell_density/cell_density_r40/CytoMAP_Sample_P7T2_leiden_25_testing.csv"



# uploaded_file1 = st.file_uploader("Choose a cvs file pretreatment")
# if uploaded_file1 is not None:
#     path1=uploaded_file1
        

# uploaded_file2 = st.file_uploader("Choose a cvs file aftertreatment")
# if uploaded_file2 is not None:
#     path2=uploaded_file2

output_dir=r"/archive/bioinformatics/Jamieson_lab/shared/CODEX/my/codex_datasets/anotated_rois/"
path1=output_dir+'P7T1+dist.csv'
path2=output_dir+'P7T2+dist.csv'
densP7T1=pd.read_csv(path1)
densP7T2=pd.read_csv(path2)
# if st.checkbox("load ROI data"):
#     pre=
#     post=
#     densP7T1=pd.read_csv(pre)
#     densP7T2=pd.read_csv(post)


 

st.write("PreRT:")
st.write(densP7T1.shape)
st.write("PostRT:")
st.write(densP7T2.shape)     

    

# densP7T1=pd.read_csv(path1)
# densP7T2=pd.read_csv(path2)
# datasets_list=[densP7T1,densP7T2]

if select_event ==  "Dataframe exploration":

    st.write("pre_RT info shape: ",densP7T1.shape)#,"pre_RT info list: ",densP7T1.columns)
    st.write("post_RT info shape:",densP7T2.shape)#,"post_RT info list:",densP7T2.columns)


    regions_maker = st.multiselect(
        'chose regions_maker',
        densP7T1.columns,
        [densP7T1.columns[0]])[0]


    cell_name = st.multiselect(
        'chose cell type name',
        densP7T1.columns,
        [ densP7T1.columns[1]])[0]


    dist_target = st.multiselect(
        'the target cell to calculate the distance',
        densP7T1.columns,
        [densP7T1.columns[0]])[0]

    T_col=[regions_maker, cell_name,dist_target] 
    st.write('You selected:', T_col)
    """
    distance unit changed from pixel to um
    """
    dist_to_tar1=get_data(densP7T1,T_col,regions_maker,cell_name)
    dist_to_tar2=get_data(densP7T2,T_col,regions_maker,cell_name)

    st.write("explore distance - cell composition: ")

    """
    define a near end/middle end/far end distance for analysis
    default value is 0.25-0.5 quantile    
    """

    def slider_setting(dist_to_tar,dist_target):
        threshmax = int(dist_to_tar[dist_target].max())
        threshmin = int(dist_to_tar[dist_target].min())
        default=[int(dist_to_tar[dist_target].quantile(0.25)),int(dist_to_tar[dist_target].quantile(0.5))]
        print(threshmin,threshmax,default)
        return [threshmin,threshmax,default]
   
    
    slider_s1=slider_setting(dist_to_tar1,dist_target)
    thresh1 = st.slider(
        'Select a range of spanning values for PreRT',
        slider_s1[0], slider_s1[1], (slider_s1[2][0], slider_s1[2][1]))
    slider_s2=slider_setting(dist_to_tar2,dist_target)
    thresh2 = st.slider(
        'Select a range of spanning values for PostRT',
        slider_s2[0], slider_s2[1], (slider_s2[2][0], slider_s2[2][1]))

    threshed1=make_dis_thres(dist_to_tar1,thresh1,dist_target)
    threshed2=make_dis_thres(dist_to_tar2,thresh2,dist_target)

    if st.checkbox("show cell precentage pre RT"):
        st.write(show_grouped_cell_with_counts(threshed1,dist_target))
    if st.checkbox("show cell percentage post RT"):
        st.write(show_grouped_cell_with_counts(threshed2,dist_target))

    if st.checkbox("dist-composition plots of PreRT"):        
        image=st.plotly_chart(show_pie_plot_of_dist_subplot(threshed1,dist_target,"pre"))
    if st.checkbox("dist-composition plots of PostRT"):
        st.plotly_chart(show_pie_plot_of_dist_subplot(threshed2,dist_target,"post"))

    if st.checkbox("distance boxenplot PreRT"):
        title="Pre RT cells distance to " + dist_target
        st.pyplot(boxen_dist_to_X(threshed1,dist_target,title))
    if st.checkbox("distance boxenplot PostRT"):
        tittle="Post RT cells distance to " + dist_target
        st.pyplot(boxen_dist_to_X(threshed2,dist_target,tittle))

    if st.checkbox("parallel categories of cell distance : pre RT"):
        parallel_categories2(get_tidy(threshed1,dist_target),style="pre",figtitle="PRE treatment: "+dist_target)
        #st.plotly_chart(parallel_categories2(get_tidy(threshed1,dist_target),style="pre",figtitle="PRE treatment: distance to tumor"))
    if st.checkbox("parallel categories of cell distance : POST RT"):
        parallel_categories2(get_tidy(threshed2,dist_target),style="post",figtitle="Post treatment: "+dist_target)
        #st.plotly_chart(parallel_categories2(get_tidy(threshed2,dist_target),style="pre",figtitle="POST treatment: distance to tumor"))

    if st.checkbox("compare distance to "+ dist_target+ "pre and post treatment"):
        #violoin_plot_distto_X(threshed1,threshed2,dist_target)
        st.pyplot(violoin_plot_distto_X(threshed1,threshed2,dist_target))




























if select_event ==  "Map exploration":

    datasetoption=st.selectbox("which dataset",["pre","post"])
    if datasetoption=="pre":
        densP7=densP7T1
        map_path=path1
    elif datasetoption=="post":
        densP7=densP7T2
        map_path=path2

    xmax=int(densP7.X_centroid.max())
    ymax=int(densP7.Y_centroid.max())

    valuesx = st.slider(
        'Select a range of values',
        0, xmax, (15000, 21000))
    valuesy = st.slider(
        'Select a range of values',
        0, ymax, (1000, 7000))
    st.write('Valuesx:', valuesx,'Valuesy:', valuesy)

  
    regions_maker = st.multiselect(
        'chose regions_maker',
        dens.columns,
        [densP7T1.columns[1]])

    # xy_ratio=(int((x_end-x_start)/7),int((y_end-y_start)/7))
    # subset_P7T1=densP7T1[(densP7T1.X_centroid>x_start) & (densP7T1.X_centroid<x_end)&(densP7T1.Y_centroid>y_start)&(densP7T1.Y_centroid<y_end)].copy()
    # print(subset_P7T1.shape)

    cell_name = st.multiselect(
        'chose cell name',
        dens.columns,
        [densP7T1.columns[1]])


    cell_ID = st.multiselect(
        'chose cell_ID',
        dens.columns,
        [densP7T1.columns[1]])



    cell_dist_to1 = st.multiselect(
        'chose cell_dist_to1',
        dens.columns,
        [densP7T1.columns[1]])


    cell_dist_to2 = st.multiselect(
        'chose cell_dist_to2',
        dens.columns,
        [densP7T1.columns[1]])


    celltype_dist_to = st.multiselect(
        'chose celltype_dist_to',
        dens.columns,
        [densP7T1.columns[1]])


    region_dist_to = st.multiselect(
        'chose region_dist_to',
        dens.columns,
        [densP7T1.columns[1]]) 
    

    selected_filename = st.text_input('selected subset file name (Duplicate names will be overwritten)', 'selected_1')+".csv"
    st.write('Your selected sebset cell info will be saved in', selected_filename)
    # save variables in yaml file for map plotting to use multiselect to get the default 
    # have to add number otherwise the list will be saved in yaml cause error in bokeh map!!     
    map_dict={"map_type":datasetoption,
              "map_path":map_path,
              "outpath": r"./interlock_map_streamlit_test.html",
              'xstart': valuesx[0],
              'xend': valuesx[1],
              'ystart': valuesy[0],
              'yend': valuesy[1],
              "region_dist_to": region_dist_to[0],
              "regions_maker": regions_maker[0],
              "cell_name" : cell_name[0],
              "cell_ID" : cell_ID[0],
              "cell_dist_to1": cell_dist_to1[0],
              "cell_dist_to2": cell_dist_to2[0],
              "celltype_dist_to": celltype_dist_to[0],
              "region_dist_to": region_dist_to[0],
              "selected_filename":selected_filename,
              
              }
    yaml_filename = "Map_param.yaml"
    # yaml_filename = st.text_input('yaml filename', 'MAP_para_1')+".yaml"
    # st.write('Your analysis parameters will be saved in', yaml_filename+".yaml")



    if st.button(":star: please save yaml before making plots :star:"):
        with open(yaml_filename, "w") as f:
            yaml.dump(map_dict, f)

    with open(yaml_filename) as f:
        yamldata = yaml.safe_load(f)

# Display the YAML data in the Streamlit app
    st.write("YAML Data:", yamldata)
    #dens_dataset = subset_P7T1
    # regions_maker="r60_20reg_allpro_noM"
    # cell_name="fulllabel_leiden2p5"
    # cell_ID="fulllabel_leiden2p5_ID"
    # cell_dist_to1="DistTo_All_T_reg"
    # cell_dist_to2="DistTo_All_tumor"
    # region_dist_to="DistTo_All_Blood_vessels"
    # celltype_dist_to="DistTo_All_Blood_vessels"
    # outpath=r"./interlock_map_streamlit_test.html"

    # def run_bokeh_server():
    #     #subprocess.Popen(["bokeh", "serve","--args","yaml_filename", "bokeh_serve.py"])
    #     process=subprocess.Popen(["bokeh", "serve", "bokeh_serve.py"])
    #     webbrowser.open_new("http://localhost:5006/bokeh_serve")

    # def stop_bokeh_server():
    #     #subprocess.run(["bokeh", "serve", "-k"])
    #     #subprocess.run([ "pkill","bokeh"])
    #     process.kill()
    #     df = pd.read_csv(selected_filename)
    #     return df
    # Read the HTML content from the file



    # chart_data = pd.DataFrame(
    #      np.random.randn(20, 3),
    #      columns=['a', 'b', 'c'])

    #st.line_chart(dens_dataset.loc[:,options])
    port="5006"
    if st.button("Show Map (may take 10-20 secs!)"):
        process=subprocess.Popen(["bokeh", "serve", "bokeh_serve.py","--port", port])
        webbrowser.open_new("http://localhost:5006/bokeh_serve")
        if st.button("Stop current bokeh to open a new one"):
            process.kill()
    if st.button("get selected data"):
        df = pd.read_csv(selected_filename)
        st.write(df)
    else:
        st.write('Byebye~')

if select_event ==  "Calculation tools":
    uploaded_calcfile = st.file_uploader("Choose a cvs subset file to be calculate")
    if uploaded_calcfile is not None:
        calc = pd.read_csv(uploaded_calcfile)
        st.write(calc) 