import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import os
import streamlit as st
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from itertools import cycle


palette_B = cycle(px.colors.qualitative.Bold)
#palette = cycle(['black', 'grey', 'red', 'blue'])
palette_2 = cycle(px.colors.sequential.Viridis)
palette_1 = cycle(px.colors.sequential.Sunset)

"""
make subset and change pixel into um 
"""
def get_data(dens_P7T1,T_col,regions_maker,cell_name):
    targeted_sebset=dens_P7T1.loc[:,T_col].copy()
    targeted_sebset.rename(columns={cell_name:'cellType',regions_maker:'regionType'},inplace=True)
    # findout dist
    dists=targeted_sebset.filter(regex='^DistTo').columns
    # change pixel into um
    if len(dists)>0:
        targeted_sebset.loc[:,dists]=targeted_sebset.loc[:,dists].multiply(0.377)
    st.write(targeted_sebset.groupby("cellType").describe())
    return targeted_sebset

"""
put threshood = distol border spanning
"""

def make_dis_thres(targeted_sebset,P7T1_thres,dist_target):
    dfcell_T1=targeted_sebset.copy()
    distol=dfcell_T1[dist_target]>P7T1_thres[1]
    border=dfcell_T1[dist_target]<P7T1_thres[0]
    spanning=(dfcell_T1[dist_target]>P7T1_thres[0])&(dfcell_T1[dist_target]<P7T1_thres[1])

    dfcell_T1.loc[distol,dist_target]='distol'
    dfcell_T1.loc[border,dist_target]='border'
    dfcell_T1.loc[spanning,dist_target]='spanning'
    dfcell_T1[dist_target+'_num']=targeted_sebset[dist_target]
    return dfcell_T1

def show_grouped_cell_with_counts(dfcell_T1,tar_column_T):
    dfcellwithcounts_1=dfcell_T1.copy()
    dfcellwithcounts_1["counts"]=np.ones((dfcell_T1.shape[0],1))

    groupedcellwithcounts1=dfcellwithcounts_1.groupby([tar_column_T,'cellType']).agg(len)["counts"]/dfcell_T1.shape[0]*100
    return(groupedcellwithcounts1)

"""show pie plot"""
def show_single_pie_plot_of_dist(dfcell_T1,tar_column_T):
    groupedcell1=dfcell_T1.groupby([tar_column_T,'cellType']).agg(len)
    border1=groupedcell1.loc["border",:]
    border1["celltype"]=border1.index
    #border1.columns=["cellType","regionType","num"]
    spanning1=groupedcell1.loc["spanning",:]
    spanning1["celltype"]=spanning1.index
    distol1=groupedcell1.loc["distol",:]
    distol1["celltype"]=distol1.index

    

    # fig = px.pie(pd.DataFrame(border1), values='num', names='cellType', color_discrete_sequence=px.colors.sequential.RdBu)

    # fig.show()

    fig = px.pie(border1, values=tar_column_T+'_num', names="celltype", color_discrete_sequence=px.colors.sequential.Sunset,width=900, height=600)

    fig.show()
    fig = px.pie(spanning1, values=tar_column_T+'_num', names="celltype", color_discrete_sequence=px.colors.sequential.Sunset,width=900, height=600)
    fig.show()
    fig = px.pie(distol1, values=tar_column_T+'_num', names="celltype", color_discrete_sequence=px.colors.sequential.Sunset,width=900, height=600)
    fig.show()
    values = groupedcell1.index.get_level_values(1)
    fig = px.pie(groupedcell1, values=tar_column_T+'_num', names=values, color_discrete_sequence=px.colors.sequential.Sunset,width=800, height=600)
    fig.show()



    # marker_colors


def show_pie_plot_of_dist_subplot(dfcell_T1,tar_column_T,style):
    
    groupedcell1=dfcell_T1.groupby([tar_column_T,'cellType']).agg(len)
    border1=groupedcell1.loc["border",:]
    border1["celltype"]=border1.index
    #border1.columns=["cellType","regionType","num"]
    spanning1=groupedcell1.loc["spanning",:]
    spanning1["celltype"]=spanning1.index
    distol1=groupedcell1.loc["distol",:]
    distol1["celltype"]=distol1.index

    

    # fig = px.pie(pd.DataFrame(border1), values='num', labels='cellType', color_discrete_sequence=px.marker_colors.sequential.RdBu)

    # fig.show()
    if style=="pre":
        markerstyle=px.colors.sequential.Sunset
    else:
        markerstyle=px.colors.sequential.Viridis
    fig = go.Figure(make_subplots(rows=1, cols=4,
            specs=[[{"type": "pie"}, {"type": "pie"},{"type": "pie"}, {"type": "pie"}]]))
    fig1 = go.Pie(labels=border1.loc[:,"celltype"],values=border1.loc[:,tar_column_T+'_num'].round(2), marker_colors = markerstyle,textfont=dict(size=10),title="border composition")
    fig2 = go.Pie(labels=spanning1.loc[:,"celltype"],values=spanning1.loc[:,tar_column_T+'_num',].round(2),marker_colors = markerstyle,textfont=dict(size=10),title="spanning composition")    
    fig3 = go.Pie(labels=distol1.loc[:,"celltype"],values=distol1.loc[:,tar_column_T+'_num'].round(2),marker_colors = markerstyle,textfont=dict(size=10),title="distol composition")   
    values1 = groupedcell1.index.get_level_values(1)
    fig0= go.Pie(labels=list(values1),values=groupedcell1.loc[:,tar_column_T+'_num'].round(2),marker_colors = markerstyle,textfont=dict(size=10),title="overall cell composition")
    fig.add_trace(fig0,row=1, col=1)
    fig.add_trace(fig1,row=1, col=2)
    fig.add_trace(fig2,row=1, col=3)
    fig.add_trace(fig3,row=1, col=4)
    fig.update_layout(legend=dict(title="Categories", x=0, y=-0.7, traceorder="normal"))
    fig.update_layout(showlegend=True,legend=dict(font=dict(size=10),orientation='h'))
    
    plt.figure(figsize=(100,50))
    # Show the figure
    #fig.show()
    return(fig)



    """boxen plot of all cell distance to tumor """

def boxen_dist_to_X(dfcell_T1,tar_column_T,title):
    plt.figure(figsize=(50,10))

    # g = sns.JointGrid(data=dfcell_T1, x="cellTyep",y="'DistTo_tumor'", space=0, ratio=17)
    # g.plot_joint(sns.scatterplot,  sizes=(30, 120),
    #         color="g", alpha=.6, legend=False)
    # g.plot_marginals(sns.rugplot, height=1, color="g", alpha=.6)

    means = dfcell_T1.loc[:,tar_column_T+"_num"].mean()
    stds = dfcell_T1.loc[:,tar_column_T+"_num"].std()

    # exclude rows that fall outside of the range mean +/- 3 * standard deviation for each column

    dfcell_T1_noOutlier = dfcell_T1.loc[(dfcell_T1[tar_column_T+"_num"] > means - 2 * stds) & (dfcell_T1[tar_column_T+"_num"] < means + 2 * stds)]
    #dfcell_T1_noOutlier[tar_column_T+"_num"]=dfcell_T1_noOutlier[tar_column_T+"_num"]*377.43/1000
    boxen1=sns.catplot(data=dfcell_T1_noOutlier.round(3), x="cellType",y=tar_column_T+"_num",kind='boxen' ,height=4, aspect=4,)# hue="variable",s=5，hue="name",
    boxen1.set_ylabels(label=tar_column_T+" (um)",size=10)
    boxen1.set_xticklabels(rotation=45,ha="right",size=10)
    boxen1.set_yticklabels(size=10)
    boxen1.set_titles(title,size=15)
    plt.title(title,fontsize=15)
    plt.show()
    return boxen1



    """violin iolin of all cell distance to xxx"""

def violoin_plot_distto_X(dfcell_T1,dfcell_T2,tar_column_T):
    plt.figure(figsize=(50,10))

    # g = sns.JointGrid(data=dfcell_T1, x="cellTyep",y="'DistTo_tumor'", space=0, ratio=17)
    # g.iolin_joint(sns.scatteriolin,  sizes=(30, 120),
    #         color="g", alpha=.6, legend=False)
    # g.iolin_marginals(sns.rugiolin, height=1, color="g", alpha=.6)

    means1 = dfcell_T1.loc[:,tar_column_T+"_num"].mean()
    stds1 = dfcell_T1.loc[:,tar_column_T+"_num"].std()

    # exclude rows that fall outside of the range mean +/- 3 * standard deviation for each column

    range_control=1 # chose 1* std or 2*std

    dfcell_T1_noOutlier = dfcell_T1.loc[(dfcell_T1[tar_column_T+"_num"] > means1 - range_control * stds1) & (dfcell_T1[tar_column_T+"_num"] < means1 + range_control * stds1)]
    # dfcell_T1_noOutlier[tar_column_T+"_num"]=dfcell_T1_noOutlier[tar_column_T+"_num"]*377.43/1000
    dfcell_T1_noOutlier = dfcell_T1_noOutlier.assign(sample_name="P7T1")


    means2 = dfcell_T2.loc[:,tar_column_T+"_num"].mean()
    stds2 = dfcell_T2.loc[:,tar_column_T+"_num"].std()

    # exclude rows that fall outside of the range mean +/- 3 * standard deviation for each column

    dfcell_T2_noOutlier = dfcell_T2.loc[(dfcell_T2[tar_column_T+"_num"] > means2 - range_control * stds2) & (dfcell_T2[tar_column_T+"_num"] < means2 + range_control * stds2)]
    # dfcell_T2_noOutlier[tar_column_T+"_num"]=dfcell_T2_noOutlier[tar_column_T+"_num"]*377.43/1000
    dfcell_T2_noOutlier = dfcell_T2_noOutlier.assign(sample_name="P7T2")

    df_all=pd.concat([dfcell_T1_noOutlier,dfcell_T2_noOutlier],axis=0)

    violin=sns.violinplot(data=df_all.round(3), x="cellType",y=tar_column_T+"_num",hue='sample_name' ,split=True,height=4, aspect=4,)# hue="variable",s=5，hue="name",
    violin.set_ylabel(tar_column_T+" um",fontsize=30)
    violin.set_xlabel("Type of cells",fontsize=30)
    violin.tick_params(axis='y', labelsize=30)
    violin.legend(fontsize=36)
    for item in violin.get_xticklabels():  # 获取X轴上的标签
        item.set_rotation(45)
        item.set_ha("right")
        item.set_size(30)
    for item in violin.get_yticklabels():  # 获取X轴上的标签
        
        
        item.set_size(30)
    # violin.set_xticklabels(rotation=45,ha="right",size=10)
    #violin.set_yticklabels(size=20)
    #plt.show()
    return(violin.get_figure())


def get_tidy(dfcell1,tar_column_T):
    dfcell1['value']=np.ones((dfcell1.shape[0],1))
    dfcell1['index1']=dfcell1.index
    dfcell_tidy1=pd.pivot(dfcell1,columns=tar_column_T,values='cellType').reset_index()

    dfcell_tidy1['cellType']=dfcell1['cellType']
    dfcell_tidy1['regionType']=dfcell1['regionType']
    dfcell_tidy1=dfcell_tidy1.fillna(0)
    return dfcell_tidy1

def parallel_categories2(dfcell_tidy1,figtitle,style):
    if style=="pre":
        fig = px.parallel_categories(dfcell_tidy1, dimensions=['cellType','border',  'spanning','distol'],
                        color='regionType', color_continuous_scale=px.colors.sequential.Inferno,
        
                        width=1000, height=800)
    elif style=="post":
        fig = px.parallel_categories(dfcell_tidy1, dimensions=['cellType','border',  'spanning','distol'],
                        color='regionType', color_continuous_scale=px.colors.sequential.Viridis, #Inferno,
        
                        width=1000, height=800)


    fig.update_layout(
        autosize = False , 
        title=figtitle,title_font=dict(size=20,color='#333333'),
        title_x=0.5, title_y=0.95,
        margin=dict(l=150, r=150, ),
        paper_bgcolor="white"#"LightSteelBlue",
    )


    fig.update_layout(coloraxis_colorbar=dict(
        # title="Number of Bills per Cell",
        # thicknessmode="pixels", thickness=50,
        # lenmode="pixels", len=200,
        # yanchor="top", y=1,
        # ticks="outside", ticksuffix=" bills",
        # dtick=5
        orientation='h'
    ))
    
    
    fig.show()
    return(fig)