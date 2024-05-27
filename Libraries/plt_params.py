import matplotlib.pyplot as plt
#import matplotlib as plt

size = 20

plt.rcParams.update({
    # taglia della figura
    'figure.figsize': (8,6),

    # font
    'text.usetex': True,
    'font.family': 'cmr10',
    'mathtext.fontset': 'cm',
    'axes.formatter.use_mathtext': True,
    'figure.titlesize' : size,

    # assi
    'axes.labelsize': size,
    'axes.titlesize': size,
    'axes.linewidth': 2,
    'axes.axisbelow': 'True',
    'axes.grid': True,
    'grid.alpha': 0.5,

    # ticks
    'xtick.labelsize': size,
    'ytick.labelsize': size,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.minor.width': 0.6,
    'ytick.minor.width': 0.6,    

    # legenda
    'legend.fontsize': size,

    # plots
    #'lines.linewidth': 1.6,
    #'scatter.marker': 'o',
    #'lines.markersize': 8,
    'scatter.edgecolors': 'black'
})