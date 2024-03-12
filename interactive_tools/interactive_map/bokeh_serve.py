import bokeh
import yaml
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, Button, BoxSelectTool

from bokeh.plotting import figure, curdoc
from bokeh.layouts import row
from bokeh import events
from bokeh.io import output_file,output_notebook, show
from bokeh.models import ColumnDataSource, LinearColorMapper, Whisker,Selection,Select, CustomJS,Button, CustomJS, Button,Div,Legend
from bokeh.plotting import figure, show
from bokeh.sampledata.penguins import data
from bokeh.transform import factor_cmap,jitter
from bokeh.models.tools import BoxSelectTool
from bokeh.palettes import Spectral6, RdBu,Purples256,Spectral3,Spectral4,Plasma256,Viridis256,Turbo256,Inferno256 #,TolRainbow3,
from bokeh.models.glyphs import VBar,Scatter,Circle
from bokeh.layouts import row, column,gridplot
import pandas as pd
import numpy as np

#from streamlit_app import *
import sys
#yaml_filename = sys.argv[1]
yaml_filename = "Map_param.yaml"

with open(yaml_filename) as f:
    data = yaml.safe_load(f)

print(data)

x_start=data["xstart"]
x_end=data["xend"]
y_start=data["ystart"]
y_end=data["yend"]
path=data["map_path"]
densP7T1=pd.read_csv(path)
xy_ratio=(int((x_end-x_start)/20),int((y_end-y_start)/20))
subset_P7T1=densP7T1[(densP7T1.X_centroid>x_start) & (densP7T1.X_centroid<x_end)&(densP7T1.Y_centroid>y_start)&(densP7T1.Y_centroid<y_end)].copy()

# from streamlit_app import regions_maker
# from streamlit_app import cell_name
# from streamlit_app import cell_ID
# from streamlit_app import cell_dist_to1
# from streamlit_app import cell_dist_to2
# from streamlit_app import region_dist_to
# from streamlit_app import celltype_dist_to
# from streamlit_app import outpath


regions_maker = data["regions_maker"]
cell_name = data["cell_name"]
cell_ID=data["cell_ID"]
cell_dist_to1=data["cell_dist_to1"]
cell_dist_to2=data["cell_dist_to2"]
region_dist_to=data["region_dist_to"]
celltype_dist_to=data["celltype_dist_to"]
outpath=data["outpath"]
selected_filename=data["selected_filename"]


#print(subset_P7T1)

def make_interlock_lasso_map(dens_dataset,regions_maker,cell_name,cell_ID,cell_dist_to1,cell_dist_to2,region_dist_to,celltype_dist_to,outpath,xy_ratio=(1200,800)):


    """data preparation """
    dens_datasetc=dens_dataset.copy()
    dens_datasetc["y_centroid"]=dens_datasetc['Y_centroid'].apply(lambda x: dens_dataset.Y_centroid.max()-x)
    dens_datasetc["regions_maker_ID"]=dens_datasetc[regions_maker].astype(int)
    dens_datasetc[regions_maker]=dens_datasetc[regions_maker].astype(str)
    dens_datasetc[cell_name]=dens_datasetc[cell_name].astype(str)

    df=dens_datasetc
    #Regions = sorted(dens_datasetc[regions_maker].unique())


    df[regions_maker]=df[regions_maker].astype(str)
    region_makers = list(sorted(df[regions_maker].unique()))
    #df.sort_values(by=region_maker, axis=0, ascending=True, inplace=True)
    Markers =  df[regions_maker].unique()


    df.sort_values(by=cell_ID, axis=0, ascending=True, inplace=True)
    Labels =  df[cell_name].unique()
    cell_names = list(sorted(df[cell_name].unique()))
    source = ColumnDataSource(df)
    options = dict(tools="pan,wheel_zoom,box_zoom,box_select,poly_select,reset,save")


    #color_map_region=factor_cmap(regions_maker, Spectral4, Regions)

    color_map_label=factor_cmap(cell_name, Spectral4, Labels)


    mapper_region = LinearColorMapper(palette=Viridis256, low=min(dens_dataset[regions_maker]), high=max(dens_dataset[regions_maker]))
    color_mapper_region={'field': regions_maker, 'transform': mapper_region}

    mapper_label = LinearColorMapper(palette=Plasma256, low=min(dens_datasetc[cell_ID]), high=max(dens_datasetc[cell_ID]))
    color_mapper_label={'field': cell_ID, 'transform': mapper_label}



    p0 = figure(width=xy_ratio[0]+280,height=xy_ratio[1],title="cell map", **options)
    p0.add_layout(Legend(),'right')
    #p1.scatter("X_centroid", "y_centroid",color={'field': 'r60_20reg_allpro_noM', 'transform': color_mapper},marker="square",source=source)
    scatter=p0.scatter(x="X_centroid",y= "y_centroid",fill_color=color_mapper_label,size=2,line_color=color_mapper_label,
            alpha=0.5, 
            selection_color='darkgreen', 
            hover_color="red",
            selection_alpha=1.0, hover_alpha=1.0, 
            legend_group=cell_name,
            source=source
            )



    p1 = figure(width=xy_ratio[0]+100,height=xy_ratio[1],title="region map", **options)
    p1.add_layout(Legend(),'right')
    #p1.scatter("X_centroid", "y_centroid",color={'field': 'r60_20reg_allpro_noM', 'transform': color_mapper},marker="square",source=source)
    scatter=p1.circle(x="X_centroid",y= "y_centroid",fill_color=color_mapper_region,size=2,line_color=color_mapper_region,
            alpha=0.5, 
            selection_color='darkred', 
            hover_color="red",
            selection_alpha=1.0, hover_alpha=1.0, 
            legend_group=regions_maker,
            source=source      
            )


    # source = ColumnDataSource(data=df)
    # renderer = p1.add_glyph(source, scatter)


    p2 = figure(width=600,height=600,
    x_axis_label=cell_dist_to1, y_axis_label=cell_dist_to2, 
    title=cell_dist_to1+" vs "+cell_dist_to2, 
    **options
    #toolbar_location="right"
    )
    p2.circle(x=cell_dist_to1,y=cell_dist_to2,  color="lightblue", line_color="white",
            alpha=0.4, 
    selection_color='darkblue', 
    hover_color="red",
    selection_alpha=1.0, hover_alpha=1.0, 
    source=source)
    # p2.width = 600
    # p2.height = 600
    div2 = Div(width=300)
    button2 = Button(label="mean dist " + cell_dist_to1 +" vs " +cell_dist_to2, width=300, height=200)

    button2.js_on_event(events.ButtonClick,  CustomJS(args=dict(div=div2), code="""
    div.text = "if poly selection: <br> one click each time to select and double click on seleted area to Calculate mean dist from selected cells to tumor / T-reg on fig3!";
    """)) 

    p2s2 = ColumnDataSource(data=dict(x=[df.loc[:,cell_dist_to1].min(), df.loc[:,cell_dist_to1].max()], ym=[df.loc[:,cell_dist_to2].mean(), df.loc[:,cell_dist_to2].mean()]))
    p2.line(x='x', y='ym', color="orange", line_width=5, alpha=0.6, source=p2s2)
    p2s3 = ColumnDataSource(data=dict(xm=[df.loc[:,cell_dist_to1].mean(), df.loc[:,cell_dist_to1].mean()], y=[df.loc[:,cell_dist_to2].min(), df.loc[:,cell_dist_to2].max()]))
    p2.line(x='xm', y='y', color="orange", line_width=5, alpha=0.6, source=p2s3)

    source.selected.js_on_change('indices', CustomJS(args=dict(s=source, s2=p2s2,s3=p2s3,cell_dist_to2=cell_dist_to2,cell_dist_to1=cell_dist_to1,div=div2), code="""
    const inds = s.selected.indices
    if (inds.length > 0) {
            const ym = inds.reduce((a, b) => a + s.data[cell_dist_to2][b], 0) /inds.length 
            const xm = inds.reduce((c, d) => c + s.data[cell_dist_to1][d], 0) /inds.length 
            s2.data = { x: s2.data.x, ym: [ym, ym]};
            s3.data = { xm: [xm,xm], y: s3.data.y};
            div.text = "Mean distance to cell y label cell:" + ym* 377.43 / 1000 +"um"+ "<br>Mean Distance to x label cell: " + xm* 377.43 / 1000 +"um"+ " <br>Number of selected points: " + inds.length;
    }      
    """))

    # Events with no attributes

    # p0.js_on_event(events.SelectionGeometry, CustomJS(args=dict(s=source, div=div2,cell_dist_to2=cell_dist_to2,cell_dist_to1=cell_dist_to1), code="""
    # var indices = s.selected.indices;
    # var sum = 0;
    # var sum1 = 0;
    # for (var i = 0; i < indices.length; i++) {
    #         sum += s.data[cell_dist_to1][indices[i]];
    #         sum1 += s.data[cell_dist_to2][indices[i]];
    # }
    # var xm = sum/indices.length;
    # var ym = sum1/indices.length;
    # div.text = "Mean of selected y-coordinates: " + ym + "Mean of selected x-coordinates: " + xm + " Number of selected points: " + indices.length;
    # """))










    p3 = figure(width=600,height=600,x_axis_label="region types",x_range=region_makers, background_fill_color="#efefef", y_axis_label=region_dist_to,title=region_dist_to+" vs type of region", **options)
    #p4.xgrid.grid_line_color = None


    #list(set(df[region_maker]))


    g3 = df.groupby(regions_maker)
    upper = g3[region_dist_to].quantile(0.80)
    lower = g3[region_dist_to].quantile(0.20)

    source3 = ColumnDataSource(data=dict(base=region_makers, upper=upper, lower=lower))

    error = Whisker(base="base", upper="upper", lower="lower", source=source3,
                    level="annotation", line_width=2)
    error.upper_head.size=30
    error.lower_head.size=30
    p3.add_layout(error)


    p3.circle(jitter(regions_maker, 0.8,range=p3.x_range),y=region_dist_to,  size=3, #color=color_mapper_label, 
            alpha=0.35, 
            selection_color='darkblue', 
            #selection_color='tomato', 
            hover_color="red",
            selection_alpha=1.0, hover_alpha=1.0, 
            source=source,
            nonselection_fill_alpha=0.2,
            color=factor_cmap(regions_maker, "Iridescent23", region_makers)
            )



    # p3 = figure(width=600,height=600,
    #         x_axis_label=regions_maker, y_axis_label=region_dist_to,
    #         title=region_dist_to+" vs regions",
    #         **options,
    #         toolbar_location="right"
    #         )
    # p3.circle(x="regions_maker_ID", y=region_dist_to, size=3, color=color_mapper_region, 
    #         fill_alpha=0.35, 
    #         selection_color='tomato', 
    #         hover_color="red",
    #         selection_alpha=1.0, hover_alpha=1.0, 
    #         source=source,
    #         nonselection_fill_alpha=0.2)
    # p3.xaxis.ticker = sorted(df["regions_maker_ID"].unique())
    #p3.scatter(celltype_dist_to, regions_maker,source,**options)


    """plot the celltype distance to some other cells"""
    p4 = figure(width=600,height=600,x_axis_label="cell types",x_range=cell_names, background_fill_color="#efefef", y_axis_label=celltype_dist_to,title=celltype_dist_to+" vs type of cell", **options)
    #p4.xgrid.grid_line_color = None
    """add wisker plot"""
    g4 = df.groupby(cell_name)
    upper = g4[celltype_dist_to].quantile(0.80)
    lower = g4[celltype_dist_to].quantile(0.20)
    source4 = ColumnDataSource(data=dict(base=cell_names, upper=upper, lower=lower))
    error = Whisker(base="base", upper="upper", lower="lower", source=source4,
                    level="annotation", line_width=2)
    error.upper_head.size=30
    error.lower_head.size=30
    p4.add_layout(error)

    p4.circle(jitter(cell_name, 0.8,range=p4.x_range),y=celltype_dist_to,  size=3, #color=color_mapper_label, 
            alpha=0.35, 
            selection_color='darkblue', 
            #selection_color='olivedrab', 
            hover_color="red",
            selection_alpha=1.0, hover_alpha=1.0, 
            source=source,
            nonselection_fill_alpha=0.2,
            color=factor_cmap(cell_name, "Category20_20", cell_names)
            )

    div4 = Div(width=300)
    button4 = Button(label="mean distance to : "+ celltype_dist_to, width=300,height=200)

    button4.js_on_event(events.ButtonClick,  CustomJS(args=dict(div=div4), code="""
    div.text = "if poly selection: <br> one click each time to select and double click on seleted area to Calculate mean distance from selected cells to your target cell type on fig4!";
    """)) 

    p4s2 = ColumnDataSource(data=dict(x=[cell_names[0],cell_names[-1]], ym=[dens_dataset.loc[:,celltype_dist_to].mean(), dens_dataset.loc[:,celltype_dist_to].mean()]))
    p4.line(x='x', y='ym', color="orange", line_width=5, alpha=0.6, source=p4s2)


    source.selected.js_on_change('indices', CustomJS(args=dict(s=source, s2=p4s2,celltype_dist_to=celltype_dist_to,div=div4), code="""
    const inds = s.selected.indices
    if (inds.length > 0) {
            const ym = inds.reduce((a, b) => a + s.data[celltype_dist_to][b], 0) / inds.length 
            
            s2.data = { x: s2.data.x, ym: [ym, ym]};
            
            div.text = "Mean distance to target cell type :<br> " + ym* 377.43/1000 + "um"+ "<br> Number of selected points: " + inds.length;
    }      
    """))





    #p4.xaxis.ticker = cell_IDs

    # p4.xaxis.major_label_overrides = dict(zip(cell_IDs,Labels))
    p4.xaxis.major_label_orientation = "vertical"

    def save_selected_data():
            # Get the selected data from the ColumnDataSource
            selected = df.iloc[source.selected.indices]
            # Convert the selected data to a pandas dataframe
            file_name = selected_filename
    
            # Save the selected data to a CSV file
            selected.to_csv(f"{file_name}", index=False)


    # Create a button widget and add it to the plot
    button = Button(label="Save Selected Data")
    button.on_click(save_selected_data)

    #p= gridplot([[p0],[p1],[p2,p3,p4]], toolbar_location="right") #, p2, p3

    row_1 = row(p0,p1)
    row_2 = row(p2,button2,div2)
    # Arrange the narrow plots in a single column
    row_3 = row(p4,button4,div4)
    row_4 = row(p3,button)

    # Arrange the plots in a grid
    p = column(row_1,row_2,row_3,row_4)

    #output_file(outpath)
    curdoc().add_root(p)

    show(p)




make_interlock_lasso_map(subset_P7T1,regions_maker,cell_name,cell_ID,cell_dist_to1,cell_dist_to2,region_dist_to,celltype_dist_to,outpath,xy_ratio)
