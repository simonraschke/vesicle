import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np



    
def makesize(width,ratio=None, scale=1):
    fig_width_pt = width  # Get this from LaTeX using \the\XXXXwidth
    inches_per_pt = 1.0/72.27                                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0
    if ratio == None:
        ratio = golden_mean
    #golden_mean = 0.54  if twocolumn == True else 0.9
    fig_width = fig_width_pt*inches_per_pt*scale                    # width in inches
    fig_height = fig_width*ratio                              # height in inches
    return [fig_width,fig_height]



def setStyle():
    # print(plt.rcParams.keys())
    plt.style.use('seaborn-paper')
    plt.rc('figure', figsize=makesize(246,ratio=0.8))
    plt.rc('figure', titleweight='light')
    plt.rc('axes', labelweight='light')
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble='\\usepackage[utf8x]{inputenc}, \\usepackage[T1]{fontenc}, \\usepackage[c]{esvect}, \\usepackage{amsfonts}, \\usepackage{amssymb}, \\usepackage{amsmath},')
    plt.rc('font', weight='light', family='serif', size=10)
    plt.rc('lines', linewidth=1)
    plt.rc('lines', markersize=3)
    plt.rc('lines', markeredgewidth=0.5)
    plt.figure().tight_layout(pad=0)
