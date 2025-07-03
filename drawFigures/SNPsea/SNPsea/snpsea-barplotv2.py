#!/usr/bin/env python3
"""
snpsea-barplot

Create a bar plot of p-values for each condition in --gene-matrix with an
adjacent heatmap all pairwise sample Pearson correlations of conditions.

Usage:
    snpsea-barplot <out> [--top INT --gene-matrix FILE --title STRING --width INT --fontsize INT --shorten INT --format STRING]

Options:
    -h, --help          Show this message and exit.
    <out>               Directory with snpsea output files.
    --top INT           Only plot the top N results [default: 25].
    --gene-matrix FILE  Override the --gene-matrix argument in args.txt.
    --title STRING      Title printed at the top of the figure.
    --width INT         Width of the figure [default: 12].
    --fontsize INT      Font size of y-axis labels [default: 12].
    --shorten INT       Shorten lengthy condition names [default: 48].
    --format STRING     Format for saving the figure (pdf, svg, png) [default: pdf].

Author:
    Kamil Slowikowski <slowikow@broadinstitute.org>
"""

from docopt import docopt
import itertools
import matplotlib as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import os
import re


def main():
    args = docopt(__doc__)

    # Set backend for headless servers
    mp.use('Agg')
    
    # Set all fonts to Times New Roman, including math text
    mp.rc('font', family='Times New Roman')
    mp.rc('mathtext', fontset='custom')
    mp.rc('mathtext', rm='Times New Roman')
    mp.rc('mathtext', it='Times New Roman:italic')
    mp.rc('mathtext', bf='Times New Roman:bold')

    out = lambda x: os.path.join(args['<out>'], x)
   
    width = float(args['--width'])
    height = width / 2.0

    matrix = args['--gene-matrix'] or find_gene_matrix(out('args.txt'))

    if not args['--title']:
        args['--title'] = os.path.basename(args['<out>'].rstrip('/'))
    
    # Get the output format
    output_format = args['--format'].lower()
    
    # Construct the output filename based on format
    output_file = f"condition_pvalues_barplot.{output_format}"

    barplot(out('condition_pvalues.txt'),
            matrix,
            out('snp_genes.txt'),
            out(output_file),
            title=args['--title'],
            figsize=(width, height),
            fontsize=float(args['--fontsize']),
            namelen=int(args['--shorten']),
            top=int(args['--top']))


def barplot(f_pvalues, f_matrix, f_genes, f_plot, figsize=(5, 5), fontsize=10,
            namelen=48, cluster_method='all', title=None, alpha=0.05, top=50):
    # Read results from snpsea
    pvalues = pd.read_table(f_pvalues, index_col=0)

    # The matrix must be in GCT format
    comp = 'gzip' if f_matrix.endswith('.gz') else None
    matrix = pd.read_table(f_matrix, 
                          compression=comp,
                          skiprows=2,
                          index_col=0).drop('Description', axis=1)

    # Bonferroni threshold calculation
    bonf = lambda a, n: 1 - (1 - a) ** (1 / float(n))

    n_samples = matrix.shape[1]
    bonferroni_threshold = -np.log10(bonf(alpha, n_samples))
    
    # Calculate uncorrected threshold (p=0.05)
    uncorrected_threshold = -np.log10(0.05)

    # Select top results
    sig = np.sum(pvalues['pvalue'] > bonferroni_threshold)
    pvalues = pvalues.sort_values('pvalue')[:top]
    matrix = matrix.loc[:, pvalues.index]

    # Handle clustering based on user genes if specified
    if cluster_method == 'user':
        snp_genelists = pd.read_table(f_genes)
        splitcomma = lambda x: str(x).split(',')
        genelists = snp_genelists['genes'].dropna().apply(splitcomma)
        snp_genes = list(set(flatten(genelists)))
        matrix = matrix.loc[snp_genes, :]

    # Create figure
    fig = plt.figure(figsize=figsize)

    # Axes for the heatmap triangle
    ax = fig.add_subplot(121, frame_on=False, aspect=2.0)

    # Get the heatmap triangle's axes and clustered sample order
    cax, order = heatmap_triangle(matrix, ax)

    # Adjust spacing
    fig.subplots_adjust(wspace=0, hspace=0, left=0, right=0.4)

    # Axes for the barplot
    ax = fig.add_subplot(122, frame_on=False)
    ax.set_axisbelow(True)

    # Order p-values by clustering
    pvalues = pvalues.loc[order]
    
    # Shorten lengthy names
    pvalues.index = [shorten(x, namelen) for x in pvalues.index]

    # Negative log10 transform
    pvalues['pvalue'] = -np.log10(pvalues['pvalue'])

    # Mark significant p-values at different thresholds
    idx_bonf_signif = pvalues['pvalue'] >= bonferroni_threshold
    idx_uncorr_signif = (pvalues['pvalue'] >= uncorrected_threshold) & (pvalues['pvalue'] < bonferroni_threshold)
    
    # Create color series with three levels
    colors = pd.Series('#aaaaaa', index=pvalues.index)  # default gray
    colors[idx_uncorr_signif] = '#99ccff'  # light blue for uncorrected significant
    colors[idx_bonf_signif] = '#5599ff'    # darker blue for Bonferroni significant

    # Get indices for labeling
    signif = pvalues.index[idx_bonf_signif | idx_uncorr_signif]

    # Create horizontal barplot
    pvalues['pvalue'].plot(ax=ax,
                          kind='barh',
                          #title=title,
                          linewidth=0,
                          grid=False,
                          color=colors)

    # Format labels
    labels = [item.get_text() for item in ax.get_yticklabels()]
    labels = ['* ' + item if item in signif else item for item in labels]
    labels = [item.replace('_', ' ') for item in labels]
    ax.set_yticklabels(labels)
    ax.tick_params(axis='y', which='major', labelsize=14)

    # Add grid and adjust ticks
    ax.grid(True, which='major', axis='both', alpha=0.5)
    xticks = np.arange(0, round(max(pvalues['pvalue'])) + 1)
    ax.set_xticks(xticks)

    # Adjust axis properties
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    ax.tick_params(length=0, axis='x')
    ax.tick_params(length=0, axis='y')

    # Set labels with proper font - 
    ax.set_xlabel('$-\log_{10}\\ P$', fontfamily='Times New Roman',fontsize=16)
    ax.set_ylabel('')

    # Add Bonferroni cutoff line (solid)
    ax.vlines(x=bonferroni_threshold,
              ymin=min(ax.get_yticks()) - 0.5,
              ymax=max(ax.get_yticks()) + 0.5,
              color='#444444',
              linestyle='-',
              linewidth=1)

    # Add uncorrected p-value cutoff line (dashed)
    ax.vlines(x=uncorrected_threshold,
              ymin=min(ax.get_yticks()) - 0.5,
              ymax=max(ax.get_yticks()) + 0.5,
              color='#444444',
              linestyle='--',
              linewidth=1)

    # Save figure with appropriate format
    fig.savefig(f_plot, bbox_inches='tight', dpi=300)


def heatmap_triangle(dataframe, axes):
    """Create a heatmap of the lower triangle of a pairwise correlation
    matrix of all pairs of columns in the given dataframe. The heatmap
    triangle is rotated 45 degrees clockwise and drawn on the given axes.

    Args:
        dataframe (pandas.DataFrame): Input data
        axes (matplotlib.axes.Axes): Matplotlib axes object
    """
    N = dataframe.shape[1]
    D = dataframe.corr(method='pearson')

    # UPGMA clustering
    Z = sch.linkage(D, method='average')
    R = sch.dendrogram(Z, no_plot=True)
    cluster_order = R['leaves']
    D = D.iloc[cluster_order, cluster_order]

    # Create lower triangle matrix
    C = np.tril(D)
    C = np.ma.masked_array(C, C == 0)
    for i in range(N):
        C[i, i] = 0

    # Rotation transformation
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = np.array([[0.5, 1], [0.5, -1]])
    A = np.dot(A, t)

    # Set up colormap
    cmap = plt.cm.RdBu_r
    norm = mp.colors.BoundaryNorm(np.linspace(-1, 1, 14), cmap.N)

    axes.set_xticks([])
    axes.set_yticks([])

    # Create heatmap
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)
    caxes = plt.pcolormesh(X, Y, np.flipud(C), axes=axes, cmap=cmap, norm=norm)

    axes.set_xlim(right=0)

    # Add colorbar
    cb = plt.colorbar(caxes, ax=axes, orientation='horizontal', shrink=0.7,
                     fraction=0.05, pad=-0.05, ticks=np.linspace(-1, 1, 5),
                     use_gridspec=True)
    cb.set_label("$\mathrm{Pearson's}\\ r$", fontfamily='Times New Roman',fontsize=16)

    return caxes, D.index


def shorten(x, n=48):
    """Shorten long strings by adding ellipsis in the middle."""
    if len(x) > n:
        return x[:n//2] + '...' + x[-n//2:]
    return x


def flatten(lists):
    """Flatten a list of lists into a single list."""
    return list(itertools.chain.from_iterable(lists))


def find_gene_matrix(filename):
    """Find the gene matrix filename in args.txt."""
    with open(filename) as f:
        for line in f:
            m = re.search('--gene-matrix\s+([^\n]+)', line)
            if m:
                return m.groups()[0].rstrip()
    return None


if __name__ == '__main__':
    main()