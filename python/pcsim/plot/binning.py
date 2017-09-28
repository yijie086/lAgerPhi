import numpy as np

def make_binning(vmin, vmax, nbins, use_log=False):
    '''Create a (nbins, edges) tuple between vmin and vmax.

    parameters:
        * vmin: minimum value
        * vmax: maximum value
        * nbins: number of bins
        * use_log: do we want a logarithmic binning?
    '''
    vmin, vmax = float(vmin), float(vmax)
    if vmin > vmax:
        vmin, vmax = vmax, vmin
    step = (vmax-vmin) / nbins
    ## note: arange does not include the last value, so add another step for or
    ##       final bin edge
    if use_log:
        vmin = np.log(vmin)
        vmax = np.log(vmax)
        step = (vmax-vmin) / nbins
    edges = []
    for i in xrange(0, nbins + 1):
        edges.append(vmin + step * i)
    if use_log:
        edges = np.exp(edges)
    else:
        edges = np.array(edges)
    return (nbins, edges)

def concat_binning(*bintuples):
    edges = np.sort(np.concatenate([b[1][:-1] for b in bintuples]))
    last_edge = np.sort(np.concatenate([b[1] for b in bintuples]))[-1]
    edges = np.append(edges, [last_edge])
    return (len(edges)-1, edges)
