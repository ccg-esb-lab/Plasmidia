import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib


def plot_matrix(df2di,xs,ys,title,cmapi,cbarlabel,plot_variables,vmini,vmaxi):
    cb_kwargs = {"shrink" : 0.75,
    "orientation" : "vertical",
    "fraction" : 0.1,
    "pad" : 0.05,
    "aspect" : 30,
    'label': cbarlabel
    }

    
    #vmin=-1
    fig=plt.figure(figsize=(len(xs)*1.0,len(ys)*0.5))
    
    ax=sns.heatmap(df2di,cmap=cmapi,vmin=vmini,vmax=vmaxi,cbar_kws=cb_kwargs)  

    ax.figure.axes[-1].yaxis.label.set_size(30)
    #xs2=[round(x,3) for x in xs]
    ax.set_xticklabels(xs,rotation='270',fontsize=20)
    #ys2=[int(x) for x in ys]
    ax.set_yticklabels(ys,fontsize=20,rotation='0' )
    plt.xlabel(plot_variables[0],fontweight='bold',fontsize=30)
    plt.ylabel(plot_variables[1],fontweight='bold',fontsize=30)
    plt.title(title,fontweight='bold',fontsize=30)
    fig

def generate_matrix_data(df,xs,ys,fixed_variablesi,fixed_valuesi,plot_variables,zval="Extinction"):
    mtx=np.zeros((len(ys),len(xs)))
    subdf=df[(df[fixed_variablesi[0]]==fixed_valuesi[0])&(df[fixed_variablesi[1]]==fixed_valuesi[1])]
    #print(subdf.head())
    for yi,y in enumerate(ys):
        for xi,x in enumerate(xs):
            #x=round(x,2)
            #y=round(y,2)
            #print(y,x)
            c=float(list(subdf[(subdf[plot_variables[0]]==y)&(subdf[plot_variables[1]]==x)][zval])[0])
            mtx[yi,xi]=c
    
    return mtx



def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
