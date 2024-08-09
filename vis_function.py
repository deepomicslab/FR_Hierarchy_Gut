import seaborn as sns
import matplotlib.pyplot as plt
import copy
import numpy as np


def fr_distribution_plot(eigen_dict, xtitle, title):
    sns.set()
    labels = []
    for key, value in eigen_dict.items():
        labels.append('cluster' + str(key))
        fr_distribution(value[1].values, xtitle)
    plt.legend(labels=labels)
    plt.title(title)


def fr_distribution(fr_matrix, xtitle):
    plt.figure()
    n_sp = fr_matrix.shape[0]
    # draw scatter plot
    plot_data = []
    for i in range(n_sp):
        for j in range(i+1, n_sp):
            if fr_matrix[i][j] > 0:
                plot_data.append(fr_matrix[i][j])
    # draw kdeplot
    ax = sns.kdeplot(data=plot_data, shade=True, legend=True)
    ax.set(xlabel=xtitle)
    plt.show()

def stack_hist(df, title, xlabel, ylabel):
    # columns is on stack
    # each hist is a row
    plt.figure((15, 15))
    fig, ax = plt.subplots()
    for col in list(df.columns):
        ax.bar(list(df.index), df[col], label=col)

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend()
    plt.show()

def sort_avg(fr_df):
    df = copy.deepcopy(fr_df)
    df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)
    df = df.reindex(df.columns, axis=0)
    return df

def cutoff_matrix(sorted_df, cutoff):
    df_new = copy.deepcopy(sorted_df)
    df_new[df_new>=cutoff] = 1
    df_new[df_new<cutoff] = 0
    df_new = sort_avg(df_new)
    return df_new

def heatmap_fr_binary(sorted_df, title, cutoff):
    plt.figure()
    df_new = copy.deepcopy(sorted_df)
    df_new[df_new>=cutoff] = 1
    df_new[df_new<cutoff] = 0
    df_new = sort_avg(df_new)
    sns.heatmap(df_new, cmap='YlGnBu_r', vmin=0, vmax=1, xticklabels=False, yticklabels=False, square=True)
    plt.title(title)
    plt.show()
    return df_new

def heatmap_fr_continue(sorted_df, title):
    plt.figure()
    df_new = copy.deepcopy(sorted_df)
    df_new = sort_avg(df_new)
    sns.heatmap(df_new, cmap='YlGnBu_r', vmin=0, vmax=1, xticklabels=False, yticklabels=False, square=True)
    plt.title(title)
    plt.show()
    return df_new

def fr_frequency(fr_matrix, xtitle):
    plt.figure()
    n_sp = fr_matrix.shape[0]
    # draw scatter plot
    plot_data = []
    for i in range(n_sp):
        for j in range(i+1, n_sp):
            if fr_matrix[i][j] > 0:
                plot_data.append(fr_matrix[i][j])
    # draw frequency
    bins = np.linspace(0, 1, 51)
    plt.hist(plot_data, bins, density=True)
    plt.ylabel('density')
    plt.title(xtitle)
    plt.show()