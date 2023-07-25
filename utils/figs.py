import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.pyplot import gcf
import matplotlib.pyplot as plt

def draw_heatmap(exprs,sample_set_dict,
                 annot=pd.DataFrame(),
                 color_dict={},
                 figsize = (20,10),dendrogram_ratio=(0.05,0.1),
                 no_legend=False,cluster_rows=True,
                 col_labels = True,row_labels = False,col_range=(-3,3)):
    '''* exprs - expressions of genes to plot
       * sample_set_dict = {"bic1":set(s1),"bic2":set(s2)}
       * annot - annotation for samples, columns - subtypes, rows - samples 
       * color_dict - how to color each label'''
    # show only variables in annotation
    #cols = [x for x in color_dict.keys() if x in list(annot.columns.values)]
    cols = list(annot.columns.values)
    s_names = []

    for s_name in sample_set_dict.keys():
        s = sample_set_dict[s_name]
        annot[s_name] = "white"
        annot.loc[list(s),s_name] = "black"
        s_names.append(s_name)
        
    sample_order = annot.index.values
    col_colors = annot.loc[:,cols+s_names]
            
    for col in reversed(cols):
        col_color_map = color_dict[col]
        col_colors[col]= col_colors[col].apply( lambda x: col_color_map[x])
        for subt in list(col_color_map.keys()):
            subt_samples = annot.loc[annot[col]==subt,:].index
            new_sample_order = [x for x in sample_order if x not in subt_samples] + [x for x in sample_order if x in subt_samples] 
            sample_order = new_sample_order
    
    for col in reversed(s_names):
        for subt in ["black","white"]:
            subt_samples = annot.loc[annot[col]==subt,:].index
            new_sample_order = [x for x in sample_order if x not in subt_samples] + [x for x in sample_order if x in subt_samples] 
            sample_order = new_sample_order
    
    vmin, vmax = col_range
    g = sns.clustermap(exprs.loc[:,sample_order],figsize=figsize,
                       col_cluster=False,row_cluster=cluster_rows,
                       dendrogram_ratio=dendrogram_ratio,colors_ratio=0.03,
                       cmap=sns.color_palette("coolwarm", as_cmap=True),
                       vmin=vmin,vmax=vmax,
                       xticklabels=col_labels, yticklabels=row_labels ,
                       col_colors=col_colors)
    ax = g.ax_heatmap
    ax.set_ylabel("")
    ax.set_xlabel("samples")
    
    g.ax_row_dendrogram.set_visible(False)
    #g.cax.set_position([.10, .2, .03, .45])
    # from https://stackoverflow.com/questions/47350879/seaborn-clustermap-subplots-adjust-cancels-colour-bar-relocation
    dendro_box = g.ax_row_dendrogram.get_position()
    dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
    g.cax.set_position(dendro_box)
    # Move the ticks to the left (https://stackoverflow.com/a/36939552/1878788)
    g.cax.yaxis.set_ticks_position("left")
    
    legends = []
    i = 0
    n_patches =0
    for col in cols:
        patches = []
        col_color_map = color_dict[col]
        # add patches only for groups found in annotation
        plot_groups  = [x for x in col_color_map.keys() if x in set(annot[col].values)]
        for group in plot_groups:
            p = g.ax_row_dendrogram.bar(0, 0, color=col_color_map[group], 
                                            label=group, linewidth=0)
            patches.append(p)
        if not no_legend:
            if len( set(annot[col].values))<=10:
                # add the legend
                legends.append(plt.legend(patches, plot_groups, loc="upper left", title=col,
                                      ncol = 10, bbox_to_anchor=(0.05*i+0.06*(n_patches), 0.95), 
                                      bbox_transform=gcf().transFigure))
                n_patches+=len(patches)
                if i>0:
                    plt.gca().add_artist(legends[i-1])
                i+=1
            else:
                # add the legend to the side
                legends.append(plt.legend(patches, plot_groups, loc='upper right', title=col,
                                      ncol = 2, #bbox_to_anchor=(0.05*i+0.06*(n_patches), 0.95), 
                                      bbox_transform=gcf().transFigure))
                if i>0:
                    plt.gca().add_artist(legends[i-1])
                i+=1
    return g, sample_order

def draw_heatmap2(exprs,biclusters,
                 annot=None,
                 color_dict=None,
                 figsize = (20,10),dendrogram_ratio=(0.05,0.1),
                 colors_ratio=(0.005,0.02),
                 no_legend=False,no_cbar=False,
                 cluster_rows=True,
                 xlabel = "samples",
                 col_labels = True,row_labels = False,
                 col_range=(-3,3),
                 bic_prefix = "bic_"):
    '''* exprs - expressions of genes to plot
       * biclusters in UnPaSt format
       * annot - annotation for samples, columns - subtypes, rows - samples 
       * color_dict - how to color each label'''
    s_names = []
    ordered_genes = []
    if type(annot)!=type(None):
        cols = list(annot.columns.values)
        row_colors = pd.DataFrame(data="white",
                                 index=exprs.index.values,
                                 columns=[bic_prefix+str(x) for x in biclusters.index.values])
    else:
        annot = pd.DataFrame(index = exprs.columns)
        row_colors = pd.DataFrame(index=exprs.index.values)
        cols = []
    if type(color_dict)!=type(None):
        pass
    
    
    for row in biclusters.iterrows():
        bic_id = bic_prefix+str(row[0])
        s = row[1]["samples"]
        g = sorted(row[1]["genes_up"])+sorted(row[1]["genes_down"])
        annot[bic_id] = "white"
        annot.loc[list(s),bic_id] = "black"
        row_colors.loc[g,bic_id] = "black"
        ordered_genes+= g
        s_names.append(bic_id)
    
    ordered_genes =  ordered_genes + sorted(set(exprs.index.values).difference(set(ordered_genes)))
        
    sample_order = annot.index.values
    col_colors = annot.loc[:,s_names+cols]
            
    for col in reversed(cols):
        col_color_map = color_dict[col]
        col_colors[col]= col_colors[col].apply( lambda x: col_color_map[x])
        for subt in list(col_color_map.keys()):
            subt_samples = annot.loc[annot[col]==subt,:].index
            new_sample_order = [x for x in sample_order if x not in subt_samples] + [x for x in sample_order if x in subt_samples] 
            sample_order = new_sample_order
    
    for col in reversed(s_names):
        for subt in ["black","white"]:
            subt_samples = annot.loc[annot[col]==subt,:].index
            new_sample_order = [x for x in sample_order if x not in subt_samples] + [x for x in sample_order if x in subt_samples] 
            sample_order = new_sample_order
    
    vmin, vmax = col_range
    g = sns.clustermap(exprs.loc[ordered_genes,sample_order],figsize=figsize,
                       col_cluster=False,row_cluster=cluster_rows,
                       dendrogram_ratio=dendrogram_ratio,colors_ratio=colors_ratio,
                       cmap=sns.color_palette("coolwarm", as_cmap=True),
                       vmin=vmin,vmax=vmax,
                       xticklabels=col_labels, yticklabels=row_labels,
                       col_colors=col_colors,
                      row_colors = row_colors)
    ax = g.ax_heatmap
    ax.set_ylabel("")
    ax.set_xlabel(xlabel)
    
    g.ax_row_dendrogram.set_visible(False)
    #g.cax.set_position([.10, .2, .03, .45])
    # from https://stackoverflow.com/questions/47350879/seaborn-clustermap-subplots-adjust-cancels-colour-bar-relocation
    dendro_box = g.ax_row_dendrogram.get_position()
    dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
    g.cax.set_position(dendro_box)
    # Move the ticks to the left (https://stackoverflow.com/a/36939552/1878788)
    g.cax.yaxis.set_ticks_position("left")
    
    if no_cbar:
        g.ax_cbar.set_visible(False)
    
    legends = []
    i = 0
    n_patches =0
    for col in cols:
        patches = []
        col_color_map = color_dict[col]
        # add patches only for groups found in annotation
        plot_groups  = [x for x in col_color_map.keys() if x in set(annot[col].values)]
        for group in plot_groups:
            p = g.ax_row_dendrogram.bar(0, 0, color=col_color_map[group], 
                                            label=group, linewidth=0)
            patches.append(p)
        if not no_legend:
            if len( set(annot[col].values))<=10:
                # add the legend
                legends.append(plt.legend(patches, plot_groups, loc="upper left", title=col,
                                      ncol = 10, bbox_to_anchor=(0.05*i+0.06*(n_patches), 0.95), 
                                      bbox_transform=gcf().transFigure))
                n_patches+=len(patches)
                if i>0:
                    plt.gca().add_artist(legends[i-1])
                i+=1
            else:
                # add the legend to the side
                legends.append(plt.legend(patches, plot_groups, loc='upper right', title=col,
                                      ncol = 2, #bbox_to_anchor=(0.05*i+0.06*(n_patches), 0.95), 
                                      bbox_transform=gcf().transFigure))
                if i>0:
                    plt.gca().add_artist(legends[i-1])
                i+=1
    return g, sample_order, (row_colors, col_colors)

def order_one(exprs,s0,subt_dict,
             subt_order = ["Her2","Basal","LumA","LumB","Normal"]):
    all_samples = set()
    for subt in subt_dict.keys():
        all_samples = all_samples.union(subt_dict[subt]) 
    
    s0 = set(s0)
    bg = all_samples.difference(s0)

    ordered_samples =[]
    ordered_sets = [s0,bg]
    for s_i in ordered_sets:
        for subt in subt_order:
            ordered_subset = list(s_i.intersection(subt_dict[subt]))
            #ordered_subset = list(exprs.loc[:,ordered_subset].sum().sort_values().index)
            ordered_samples+= ordered_subset
    return ordered_samples
    
def order_two(s0,s1, subt_dict,
             subt_order = ["Her2","Basal","LumA","LumB","Normal"]):
    all_samples = set()
    for subt in subt_dict.keys():
        all_samples = all_samples.union(subt_dict[subt]) 
    s0 = set(s0)
    s1 = set(s1)

    overlap = s0.intersection(s1)
    s0_ = s0.difference(s1)
    s1_ = s1.difference(s0)
    bg = all_samples.difference(s0).difference(s1)

    ordered_samples =[]
    ordered_sets = [overlap,s0_,s1_,bg]
    for s_i in ordered_sets:
        for subt in subt_order:
            ordered_subset = list(s_i.intersection(subt_dict[subt]))
            ordered_samples+= ordered_subset
    return ordered_samples