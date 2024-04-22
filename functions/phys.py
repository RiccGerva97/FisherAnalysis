import numpy as np

def PacMan(x, d = np.array((0, 0, 1000.)) ):
    """Returns an array x with the 3rd component in the interval [0; d]"""
    if 0 <= x[2] <= d[2]: return x
    elif x[2] > d[2]: return PacMan(x - d)
    elif x[2] < 0: return PacMan(x + d)

def HubblePar(z, O_m = 0.3175, H_0 = 67.11, O_k = 0, O_de = -1, w = -1):
    """Evaluates the Hubble parameter depending on redshift `z` and cosmological parameters:
       - `O_m` matter density;
       - `O_k` curvature density;
       - `O_de` dark matter density;
       - `w` time exponent of dark energy."""
    if H_0 < 1.: H_0 *= 100.
    if O_de == -1: O_de = 1 - (O_m + O_k)
    H2 = H_0**2 * (O_m * (1+z)**3 +\
                  O_k * (1+z)**2 +\
                  O_de* (1+z)**(3*(1+w))  )
    return np.sqrt(H2)
