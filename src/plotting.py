import itertools
import json
import os

import colorcet as cc
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import ternary
from cycler import cycler
from scanpy import api as sc

from downstream.src.visualization import visualize
from icat.src import utils

loc = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(loc, 'configs/figures.mplstyle'))
plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))

method_dictionary = {
    'icat': 'icat',
    'seurat233': 'Seurat 2.3',
    'seurat.aligned': 'Seurat 2.3 - Aligned',
    'seurat311': 'Seurat 3.1',
    'scanorama': 'scanorama',
    'icat_scan': 'scanorama + icat',
    'seurat_icat': 'Seurat 3.1 + icat'
}

method_colors = ["#0a5e62", "#4bbfb8", "#fbcf5b",
                 "#ff5959", "#4f8522", "#a2d665"]

method_colors1 = ["#3c4347", "#5e6769", "#719192",
                  "#e1cec5", "#518413", "#a3c541"]

method_colors2 = ["#f77855", "#584b43", "#537d90",
                  "#a5d2c9", "#518413", "#a3c541"]

metric_dictionary = {
    'adjusted.mutual.info': 'AMI',
    'adjusted.rand': 'ARI',
    'completeness': 'Completeness',
    'fowlkes.mallows': 'Fowlkes-Mallows',
    'homogeneity': 'Homogeneity'
}

def stacked_barplot(df, label, cluster, xlabel=''):
    """
    Plot a stacked barplot.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing two columns of sample labels.
    label : str
        Name of column in `df` containing known labels.
    cluster : str
        Name of column in `df` containing predicted labels.
    xlabel : str
        Description of group plotted along the x-axis.
    """
    counts = df.groupby([cluster, label]).size().unstack().fillna(0)
    cluster_totals = counts.sum(axis=0)
    # percentages of cells belonging to each cluster for each known label 
    # (e.g. 87% of cells in known label X are in cluster Y)
    percentages = (counts / cluster_totals).T * 100
    labels = percentages.index.values
    clusters = percentages.columns.values
    xticks = range(percentages.shape[0])
    totals = np.zeros(percentages.shape[0])
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    # plot barplots up to calculated percentage for each cluster in known labels
    for each, color in zip(clusters, colors()):
        new_percentages = percentages.loc[:, each].values
        plt.bar(xticks, new_percentages, color=color['color'],
                width=0.85, bottom=totals, label=each)
        # update cumulate percentage for starting points
        totals += new_percentages
    legend_cols = int(np.ceil(len(clusters) / 20))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
               ncol=legend_cols)
    plt.xticks(xticks, labels, rotation=90)
    plt.xlabel(xlabel, fontsize=24, labelpad=5)
    plt.ylabel("Percentage of Cells", fontsize=24)
    yticks, __ = plt.yticks()
    ylabels = ["{:d}%".format(int(y)) for y in yticks]
    plt.yticks(yticks, ylabels)
    plt.ylim(0, 100)
    plt.tight_layout()

def ranked_heatmap(ranks):
    """
    Plot a heatmap of method ranks across various method.

    Low values are considered good (e.g. 1 is better than 2). Negative values
    are used during plotting for visual purposes.
    
    Parameters
    ----------
    ranks : pd.DataFrame
        A dataframe containing methods as index values and metrics as columns.
        Each cell represents the rank of the specified method in the respective
        metric compared to other methods.
    
    Returns
    -------
    mpl.Figure
        Heatmap summarizing method performance.
    """
    figsize = np.array(plt.rcParams['figure.figsize']) * 1.25
    labelsize = plt.rcParams['axes.titlesize']
    fig, ax = plt.subplots(figsize=figsize)
    plot_ranks = -1 * ranks.loc[ranks.mean(axis=1).sort_values().index, :]
    sns.heatmap(plot_ranks, cmap='viridis',
                cbar_kws={'shrink': 0.95})
    ax.set_yticklabels([method_dictionary[x.get_text()]\
                        for x in ax.get_yticklabels()],
                       fontsize=labelsize*0.75, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=labelsize*0.75,
                       rotation=90)
    fig.axes[1].annotate('Worse', (-1, -0.05), xycoords='axes fraction',
                         fontsize=labelsize*0.7)
    fig.axes[1].annotate('Better', (-1, 1.01), xycoords='axes fraction',
                         fontsize=labelsize*0.7)
    fig.axes[1].set_yticklabels([''] * len(fig.axes[1].get_yticks()) )
    plt.ylabel('')
    plt.tight_layout()
    return fig


def trendplot(results, x, y, hue=None, xlabel=None):
    """
    Plot performance trends across some variable.
    
    Parameters
    ----------
    results : pd.DataFrame
        Dataframe containing performance measures and other features across
        datasets and methods.
    x : string
        Column in `results` representing the independent variable in cluster
        performance.
    y : string
        Column in `results` measuring cluster performance -- presumably affected
        by `x`.
    hue : string, optional
        Column in `results` to separate on. Default is None, and all points rows
        will be used in a single plot.
    xlabel : string, optional
        Label for x-axis, by default None, and `x` will be used.
    
    Returns
    -------
    mpl.Axes
        Plot of linear trends between `x` and `y` separated on `hue`.
    """
    lm = sns.lmplot(x=x, y=y, hue=hue, col=hue, col_wrap=3, data=results)
    xmin = min(results[x])
    xmin -= 0.1*xmin
    xmax = max(results[x])
    xmax += 0.1*xmax
    labelsize= plt.rcParams['axes.titlesize']
    plt.ylim(-0.05, 1)
    for ax in lm.fig.axes:
        ax.set_title(ax.get_title().replace("{} = ".format(hue), ""),
                     fontsize=int(labelsize * 0.75))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xlim(xmin, xmax)
    if xlabel is None:
        xlabel = x
    lm.fig.text(0.5, 0.03, xlabel, ha='center', va='center',
                fontsize=labelsize)
    lm.fig.text(0.02, 0.5, y, ha='center', va='center', rotation='vertical',
                fontsize=labelsize)
    return lm


def flip(items, ncol):
    """
    Flips matplotlib legends to increment by rows before columns.

    Taken from here:
    https://stackoverflow.com/questions/10101141/matplotlib-legend-add-items-across-columns-instead-of-down
    
    Parameters
    ----------
    items : [type]
        [description]
    ncol : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def metric_plot(scores, errors=None):
    """
    Plot performance metrics across methods. 
    
    Parameters
    ----------
    scores : pd.DataFrame
        Performance data frame where each column is a different performance
        metric, and each row represents a different method.
    
    errors : pd.DataFrame, optional
        Errors associated with each method-performance pair in `scores`. Should
        be the same shape/size as `scores`. 
    """
    if errors is not None:
        if not isinstance(errors, pd.DataFrame):
            raise ValueError("Expected DataFrame for `errors` parameter")
        if np.all(errors.index != scores.index):
            raise ValueError("`scores` and `errors` dataframes should have the"
                             " same index.")
        if np.all(errors.columns != scores.columns):
            raise ValueError("`scores` and `errors` dataframes should have the"
                             " same columns.")
    indices = np.arange(scores.shape[0]).astype(float)
    width = np.min(np.diff(indices)) / scores.shape[1]
    indices = indices + indices * width
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    starts = indices - width * scores.shape[1] / 2
    for method, color in zip(sorted(scores.columns), colors()):
        yerr = None
        if errors is not None:
            yerr = errors[method]
        plt.bar(starts + width, scores[method], width, color=color['color'],
                label=method_dictionary[method],
                yerr=yerr)
        starts = starts + width
    indices = indices + width / 2
    plt.xticks(indices, labels=[metric_dictionary[x] for x in scores.index.values])
    plt.title("Method Performance", loc='left')
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    legend_cols = int(np.ceil(scores.shape[1] / 2))
    plt.legend(flip(handles, legend_cols), flip(labels, legend_cols),
               loc='upper center', bbox_to_anchor=(0.5, -0.05),
               ncol=legend_cols, frameon=False, fancybox=False)
    plt.ylim(0, 1)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()

def close_plot():
    plt.cla()
    plt.clf() 
    plt.close()


def ternary_color_point(point, scale=9):
    """
    Transform a ternary point into an associated RGBV value.
    
    Parameters
    ----------
    point : list-like, numpy.ndarray
        A ternary point. Expected to have length equal to three.
    scale : int, optional
        Total value points sum to, by default 9.
    
    Returns
    -------
    np.ndarray
        Color associated with the provided point in RGBV space.
    """
    point = np.array(point)
    assert point.size == 3
    scaled = (point + 3) / scale
    scaled[scaled < 0] = 0
    cmaps = [plt.get_cmap('Blues'), plt.get_cmap("Reds"), plt.get_cmap('Greens')]
    rgba = np.array([cmap(p) for cmap, p in zip(cmaps, scaled)])
    out = (rgba.mean(axis=0)) * np.hstack((scaled, [1])) 
    out[-1] = 1
    return out


# [9, 0, 0] 'H1975', -> bottom right 
# [0, 9, 0] 'H2228', -> top
# [0, 0, 9] 'HCC827' -> bottom left
def ternary_plot(data, column, label, scale=9):
    """
    Summarize `benchmark` clusters by plotting cluster medians
    
    Parameters
    ----------
    data : pd.DataFrame
        Dateframe containing assigned cluster labels along with known cell
        mixtures.
    column : string
        Column in `data` containing assigned cluster labels.
    label : string
        Column in `data` containing known cell mixtures.
    scale : int, optional
        Total number of cells used in each mixture, by default 9.
    
    Returns
    -------
    mpl.Axes
        Ternary plot of cluster medians of cell mixtures.
    """
    data['c1'] = data.apply(lambda x: int(x[label].split(',')[0]), axis=1)
    data['c2'] = data.apply(lambda x: int(x[label].split(',')[1]), axis=1)
    data['c3'] = data.apply(lambda x: int(x[label].split(',')[2]), axis=1)
    by_column = data.groupby(column)
    plot_data = data.groupby(column)['c1', 'c2', 'c3'].median()
    size = by_column.size()
    plot_data['size'] = ((950 - 175) * (size - 45) /\
                         (950 - 45) + 175).astype(int)
    __, tax = ternary.figure(scale=scale)
    colors = plot_data.apply(lambda x: ternary_color_point(
                                                   np.array([x.c1, x.c2, x.c3]),
                                                   scale),
                             axis=1).values
    labelsize=plt.rcParams['axes.titlesize']
    tax.gridlines(multiple=1)
    tax.boundary(scale=scale, linewidth=1.25)
    tax.scatter(plot_data[['c1', 'c2', 'c3']].values, c=colors,
                s=plot_data['size'], alpha=1)
    tax.left_corner_label('HCC827', position=[-0.1, 0.075, 0],
                          fontsize=labelsize)
    tax.top_corner_label('H2228', position=[-0.02, 1.17, 0],
                          fontsize=labelsize)
    tax.right_corner_label('H1975', position=[1.02, 0.02, 0],
                          fontsize=labelsize)
    tax.ticks(axis='lbr', multiple=1, offset=0.025,
              fontsize=int(labelsize * 0.75), linewidth=1.25)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    return tax
