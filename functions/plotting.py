import numpy as np
import sys
from scipy.stats import norm, chi2
sys.path.insert(1, '../dicts')
from dicts.my_dicts import cov_scale

def plot_ellipse(cov_plot, q=None, nsig=None, **kwargs):
    """
    Returns the erguments to print a `matplotlib.patches` `Ellipse` calculating the
    eighenvalues

    Parameters
    ----------
    cov : (2, 2) array
        Covariance matrix.
    q : float, optional
        Confidence level, should be in (0, 1)
    nsig : int, optional
        Confidence level in unit of standard deviations. 
        E.g. 1 stands for 68.3% and 2 stands for 95.4%.

    Returns
    -------
    width, height, rotation :
         The lengths of two axises and the rotation angle in degree
    for the ellipse.
    """

    if q is not None:
        q = np.asarray(q)
    elif nsig is not None:
        q = 2 * norm.cdf(nsig) - 1
    else:
        raise ValueError('One of `q` and `nsig` should be specified.')
    r2 = chi2.ppf(q, 2)

    # def eigsorted(cov_plot):
    #     vals, vecs = np.linalg.eigh(cov_plot)
    #     order = vals.argsort()[::-1]
    #     return vals[order], vecs[:,order]
    # vals, vecs = eigsorted(cov_plot)
    # width, height = 2 * nsig * np.sqrt(vals)
    # rotation = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    val, vec = np.linalg.eigh(cov_plot)
    width, height = 2 * np.sqrt(np.abs(val[:, None] * r2))
    rotation = np.degrees(np.arctan2(*vec[::-1, 0]))

    return width, height, rotation


def my_ellipse(cov, q=None):
    """
    Returns the erguments to print a `matplotlib.patches` `Ellipse` using the analytical forulae
    Arguments:
    - `cov`: the covariance matrix
    - `q`: the `int` key of the confidence level dictionary `cov_scale`
    Returns:
    - width
    - height
    - angle
    """
    a, c = cov[0, 0], cov[1, 1]
    b = np.mean([cov[0, 1], cov[1, 0]])
    l1 = (a+c)/2 + np.sqrt( ((a-c)/2)**2 + b**2 )
    l2 = (a+c)/2 - np.sqrt( ((a-c)/2)**2 + b**2 )
    theta = np.degrees(np.arctan2(l1-a, b))

    if q is None:
        return 2*l1, 2*l2, theta
    else:
        s = cov_scale[int(q)]
        return 2*np.sqrt(np.abs(l1)*s), 2*np.sqrt(np.abs(l2)*s), theta