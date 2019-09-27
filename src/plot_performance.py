from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
from scanpy import api as sc
import numpy as np
from cycler import cycler

method_dictionary = {
    'icat': 'icat',
    'seurat': 'Seurat',
    'scanorama': 'scanorama',
    'icat_scan': 'scanorama + icat'
}

metric_dictionary = {
    'adjusted.mutual.info': 'AMI',
    'adjusted.rand': 'ARI',
    'completeness': 'Completeness',
    'fowlkes.mallows': 'Fowlkes-Mallows',
    'homogeneity': 'Homogeneity'
}


def stacked_barplot(df, label, cluster):

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
    for each, color in zip(clusters, colors()):
        new_percentages = percentages.loc[:, each].values
        plt.bar(xticks, new_percentages, color=color,
                width=0.85, bottom=totals, label=each)
        totals += new_percentages
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(xticks, labels, rotation=90)
    plt.tight_layout()


def plot_performance(df):
    indices = np.arange(df.shape[0]).astype(float)
    width = np.min(np.diff(indices)) / df.shape[1]
    indices = indices + indices * width
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    starts = indices - width * df.shape[1] / 2
    for method, color in zip(df.columns, colors()):
        plt.bar(starts + width, df[method], width, color=color['color'],
                label=method_dictionary[method])
        starts = starts + width
    indices = indices + width / 2
    plt.xticks(indices, labels=[metric_dictionary[x] for x in df.index.values])
    plt.title("Method Performance", loc='left')
    plt.legend()
#     plt.legend(loc='upper center', bbox_to_anchor=(0, -0.05),
#                ncol=df.shape[1])
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
#     plt.tight_layout()