import numpy as np

def make_binning(vmin, vmax, nbins, use_log=False):
    '''Create a (nbins, binning) tuple between vmin and vmax.

    parameters:
        * vmin: minimum value
        * vmax: maximum value
        * nbins: number of bins
        * use_log: do we want a logarithmic binning?
    '''
    vmin, vmax = float(vmin), float(vmax)
    if vmin > vmax:
        vmin, vmax = vmax, vmin
    bins = np.arange(vmin, vmax, (vmax - vmin) / (nbins+1))
    if use_log:
        bins = np.exp(np.arange(np.log(vmin), np.log(vmax),
                        np.log(vmax/vmin)/(nbins+1)))
    return (nbins, bins)

