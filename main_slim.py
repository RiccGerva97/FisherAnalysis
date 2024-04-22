import numpy as np
import sys, os, h5py
from matplotlib import rcParams
from scipy.stats import norm, chi2

sys.path.insert(1, './dicts')
sys.path.append('./functions')
from dicts.my_dicts import order_folders, cosmological_pars, order_dimension, VarCosmoPar
from functions.parsers import cosmo_parser
from functions.math import Hartlap
from functions.pcg import randind

# name of the files typology
TYPE_FILE = '3D_cubes'

root = f"/home/{TYPE_FILE}_files/"
files_to_read = os.listdir(root)

# if you want to select a specified number of fiducial realizations
max_fids = 10000
# order of the neutrino mass derivate
order_derivate = 4

## reading estimator coeffs.
# number of estimator coeffs
num_coeff = 100
every_mean = np.zeros((len(order_folders), num_coeff))    # define empty total array
errbar = np.zeros((len(order_folders), num_coeff))        # define error bar array

for FC in files_to_read:
    every_coeff = []

    # identify cosmology model
    cosmo = cosmo_parser(FC)
    
    # `index`: give always the same intdex to the models
    index = order_folders[cosmo]
    
    # read .hdf5 file
    with h5py.File(root+FC, 'r') as fh:
        for k in fh.keys():
            dat = np.array(fh[k])
            every_coeff.append(dat)
    every_mean[index] = np.mean(every_coeff, axis=0)
    errbar[index] = np.std(every_coeff, axis=0)


# ------------- DERIVATES -------------------------------------------------------------------------------

n_pars = len(cosmological_pars)
derivates = np.zeros((n_pars, num_coeff))

# in evaluating the derivates using Quijote, you must take into account two facts:
#   1) 'Mnu' has 3 options of derivate (3 orders)
#   2) 'Ob' has been in(de)cremented both once and twice; here I only have used the double variation
#        a) if you want/have to use only the single variation simply comment the last 'elif' and modify
#            the first 'if'
#        b) if you want/have to use both, I simply haven't take into account this possibility XD

for i in cosmological_pars:
    if "Mnu" not in i and "Ob" not in i:
        ind = order_dimension[i]
        derivates[ind]    = (every_mean[order_folders[i+"_p"]] - every_mean[order_folders[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
    elif "Mnu" in i:
        # # >>> ORDER 1 >>>
        if order_derivate == 1:
            print('first')
            derivates[order_dimension['Mnu']]    = (every_mean[order_folders["Mnu_p"]] - every_mean[order_folders["zeldovich"]]) / (0.1)
        # >>> ORDER 2 >>>
        elif order_derivate == 2:
            print('second')
            derivates[order_dimension['Mnu']] = (-every_mean[order_folders['Mnu_pp']] + 4 * every_mean[order_folders['Mnu_p']] - 3 * every_mean[order_folders['zeldovich']]) / (2 * 0.1)
        # # >>> ORDER 4 >>>
        elif order_derivate == 4:
            print('third')
            derivates[order_dimension['Mnu']] = (every_mean[order_folders['Mnu_ppp']] - 12 * every_mean[order_folders['Mnu_pp']] + 32 * every_mean[order_folders["Mnu_p"]] - 21 * every_mean[order_folders['zeldovich']]) / (12 * 0.1)
        else:
            assert False, 'Error in Mnu derivate: must choose order in [1; 2; 4]'
    elif "Ob" in i:
        derivates[order_dimension['Ob']]    = (every_mean[order_folders[i+"2_p"]] - every_mean[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])

# ---------------- COEFFICIENTS COVARIANCE MATRIX ----------------------------------------------------------------------------------------

# initializing empty fiducial arrays
fiducial = []

## Read WST files
with h5py.File(root+'fiducial_est_real.hdf5', 'r') as fh:
    for k in fh.keys():
        dat = np.array(fh[k])
        fiducial.append(dat)

# if you choose a minor number of fiducial realizations, here the diluition is executed
if max_fids < len(fiducial):
    new_indices = randind(len(fiducial, max_fids))
    fiducial = fiducial[new_indices,:]

# transpose arrays
fiducial = np.array(fiducial).transpose()

# covariance matices
cov = np.cov(fiducial)
    
# invese matrices with Hartlap formalism
Hart = Hartlap(cov, max_fids)

# initialize empty FIMs
Fish = np.zeros((n_pars,n_pars))
# create FIMs (approx form)
for a in range(6): 
    for b in range(6):
        Fish[a, b] = np.sum(np.dot(derivates[a], np.dot(Hart, np.transpose(derivates[b]))))

# get the inverse FIM aka (co)variances
inverse = np.linalg.inv(Fish)

# write results on file
with h5py.File('./inv_FIM_est.hdf5', 'w') as wf:
    dataset = wf.create_dataset('inv_FIM_est', data = inverse)
