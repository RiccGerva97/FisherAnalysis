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

# if `True` rsd have been calculated
bool = False

root_Pk  = f"/home/Pk_{TYPE_FILE}_files/"
files_to_read_pk = os.listdir(root_Pk)

# if you want to select a specified number of fiducial realizations
max_fids = 10000
# order of the neutrino mass derivate
order_derivate = 4

## reading Pk values
# it's assumed you calculate three multipoles, if not comment/delete the unnecessary lines...
# -> AND modify the `3*` factor in the rsd derivate of power spectrum

# number of Pk coeffs
num_pk = 96
every_mean_m_pk = np.zeros((len(order_folders), num_pk)) # monopole
every_mean_q_pk = np.zeros((len(order_folders), num_pk)) # quadrupole
every_mean_h_pk = np.zeros((len(order_folders), num_pk)) # hexadecapole
sn = np.zeros(len(order_folders))                        # shot noise for veriation realizations
k = []                                                   # expected to be the same for all realizations
errbar_m_pk = np.zeros((len(order_folders), num_pk))
errbar_q_pk = np.zeros((len(order_folders), num_pk))

for FC in range(len(files_to_read_pk)):
    every_power_m, every_power_q, every_power_h = [], [], []
    every_sn = []
    currfile = files_to_read_pk[FC]
    cosmo = cosmo_parser(files_to_read_pk[FC])
    assert cosmo in files_to_read_pk[FC]
    index = order_folders[cosmo]
    with h5py.File(root_Pk+currfile, 'r') as f:
        for key in f.keys(): #tqdm(f.keys()):
            dat = np.array(f[key])
            if   'mono' in key: every_power_m.append(dat)
            elif 'quad' in key: every_power_q.append(dat)
            elif 'hexa' in key: every_power_h.append(dat)
            elif 'k_coord' in key: k = dat
            elif 'sn'    in key: every_sn.append(dat)
    every_mean_m_pk[index] = np.mean(every_power_m, axis=0)
    every_mean_q_pk[index] = np.mean(every_power_q, axis=0)
    every_mean_h_pk[index] = np.mean(every_power_h, axis=0)
    sn = np.mean(every_sn, axis=0)
    errbar_m_pk[index] = np.std(every_power_m, axis=0)
    errbar_q_pk[index] = np.std(every_power_q, axis=0)

# remove shot noise
every_mean_m_pk = every_mean_m_pk - sn[:, None]

# every_mean_poles = np.zeros((len(order_folders), 3*num_pk))
# for i in order_folders:
#     ind = order_folders[i]
#     every_mean_poles[ind] = np.concatenate((every_mean_m_pk[ind], every_mean_q_pk[ind], every_mean_h_pk[ind]))



# ------------- DERIVATES -------------------------------------------------------------------------------

n_pars = len(cosmological_pars)

derivates_pk  = np.zeros((n_pars, num_pk))
# derivates_poles_pk  = np.zeros((n_pars, 3*num_pk))    # <- MODIFY if haven't 3 multip.

# in evaluating the derivates using Quijote, you must take into account two facts:
#   1) 'Mnu' has 3 options of derivate (3 orders)
#   2) 'Ob' has been in(de)cremented both once and twice; here I only have used the double variation
#        a) if you want/have to use only the single variation simply comment the last 'elif' and modify
#            the first 'if'
#        b) if you want/have to use both, I simply haven't take into account this possibility XD 

for i in cosmological_pars:
    if "Mnu" not in i and "Ob" not in i:
        ind = order_dimension[i]
        derivates_pk[ind] = (every_mean_m_pk[order_folders[i+"_p"]] - every_mean_m_pk[order_folders[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
        # derivates_poles_pk[ind]  = (every_mean_poles[order_dimension[i+"_p"]] - every_mean_poles[order_dimension[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
    elif "Mnu" in i:
        # # >>> ORDER 1 >>>
        if order_derivate == 1:
            derivates_pk[order_dimension['Mnu']] = (every_mean_m_pk[order_folders["Mnu_p"]] - every_mean_m_pk[order_folders["zeldovich"]]) / (0.1)
            # derivates_poles_pk[order_dimension['Mnu']] = (every_mean_poles[order_folders["Mnu_p"]] - every_mean_poles[order_folders["zeldovich"]]) / (0.1)
        # >>> ORDER 2 >>>
        elif order_derivate == 2:
            derivates_pk[order_dimension['Mnu']] = (-every_mean_m_pk[order_folders['Mnu_pp']] + 4 * every_mean_m_pk[order_folders["Mnu_p"]] - 3 * every_mean_m_pk[order_folders['zeldovich']]) / (2 * 0.1)
            # derivates_poles_pk[order_dimension['Mnu']] = (-every_mean_poles[order_folders['Mnu_pp']] + 4 * every_mean_poles[order_folders["Mnu_p"]] - 3 * every_mean_poles[order_folders['zeldovich']]) / (2 * 0.1)
        # # >>> ORDER 4 >>>
        elif order_derivate == 4:
            derivates_pk[order_dimension['Mnu']] = (every_mean_m_pk[order_folders['Mnu_ppp']] - 12 * every_mean_m_pk[order_folders['Mnu_pp']] + 32 * every_mean_m_pk[order_folders["Mnu_p"]] - 21 * every_mean_m_pk[order_folders['zeldovich']]) / (12 * 0.1)
            # derivates_poles_pk[order_dimension['Mnu']] = (every_mean_poles[order_folders['Mnu_ppp']] - 12 * every_mean_poles[order_folders['Mnu_pp']] + 32 * every_mean_poles[order_folders["Mnu_p"]] - 21 * every_mean_poles[order_folders['zeldovich']]) / (12 * 0.1)
        else:
            assert False, 'Error in Mnu derivate: must choose order in [1; 2; 4]'
    elif "Ob" in i:
        derivates_pk[order_dimension['Ob']] = (every_mean_m_pk[order_folders[i+"2_p"]] - every_mean_m_pk[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])
        # derivates_poles_pk[order_dimension['Ob']] = (every_mean_poles[order_folders[i+"2_p"]] - every_mean_poles[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])



# ---------------- COEFFICIENTS COVARIANCE MATRIX ----------------------------------------------------------------------------------------

# initializing empty fiducial arrays
fiducial_m_pk = []
fiducial_q_pk = []
fiducial_h_pk = []
k, f_sn = [], []

## Read Pk files
with h5py.File(root_Pk+'fiducial_pk_real.hdf5', 'r') as f:
    for k in f.keys():
        dat = np.array(f[k])
        if   'mono'  in k: fiducial_m_pk.append(dat)
        elif 'quad'  in k: fiducial_q_pk.append(dat)
        elif 'hexa'  in k: fiducial_h_pk.append(dat)
        elif 'k_coord' in k: k.append(dat)
        elif 'sn'    in k: f_sn.append(dat)
        else: assert False, f'keyword {k} not found in `if` code'

# if you choose a minor number of fiducial realizations, here the diluition is executed
if max_fids < len(fiducial_m_pk):
    new_indices = randind(len(fiducial_m_pk, max_fids))
    fiducial_m_pk = fiducial_m_pk[new_indices,:] 
    fiducial_q_pk = fiducial_q_pk[new_indices,:]
    fiducial_h_pk = fiducial_h_pk[new_indices,:] 
    k = k[new_indices,:] 
    f_sn = f_sn[new_indices,:]

# remove shot noise
if not  np.shape(sn) == ():
    fiducial_m_pk = fiducial_m_pk - f_sn[:, None]

# fiducial_poles = np.concatenate((fiducial_m_pk-sn, fiducial_q_pk, fiducial_h_pk))

# transpose arrays
fiducial_m_pk = np.array(fiducial_m_pk).transpose()
# fiducial_poles = np.array(fiducial_poles).transpose()

# covariance matices
cov_pk = np.cov(fiducial_m_pk)
# cov_poles = np.cov(fiducial_poles)
    
# invese matrices with Hartlap formalism
Hart_pk  = Hartlap(cov_pk, max_fids)
# Hart_poles  = Hartlap(cov_poles, max_fids)

# initialize empty FIMs
Fish_pk = np.zeros((n_pars,n_pars))
# Fish_poles = np.zeros((n_pars,n_pars))
# create FIMs (approx form)
for a in range(6): 
    for b in range(6):
        Fish_pk[a, b] = np.sum(np.dot(derivates_pk[a], np.dot(Hart_pk, np.transpose(derivates_pk[b]))))
        # Fish_poles[a, b] = np.sum(np.dot(derivates_poles_pk[a], np.dot(Hart_poles, np.transpose(derivates_poles_pk[b]))))

# get the inverse FIM aka (co)variances
inverse_pk = np.linalg.inv(Fish_pk)
# inverse_poles = np.linalg.inv(Fish_poles)

# write results on file
with h5py.File('./inv_FIM_est.hdf5', 'a') as wf:
    dataset_pk_real = wf.create_dataset('inv_pk_real', data = inverse_pk)
    # dataset_poles_real = wf.create_dataset('inv_poles_real', data = inverse_poles)