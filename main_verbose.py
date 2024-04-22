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

print(f'Progress:  |-------|  {round(0/21*100, 2)}%', end='\r')

# name of the files typology
TYPE_FILE = '3D_cubes'

# if `True` rsd have been calculated
bool_rsd = True

root_WST = f"/home/WST_{TYPE_FILE}_files/"
root_Pk  = f"/home/Pk_{TYPE_FILE}_files/"
files_to_read_wst = os.listdir(root_WST)
files_to_read_pk = os.listdir(root_Pk)

# if you want to select a specified number of fiducial realizations
max_fids = 10000
# order of the neutrino mass derivate
order_derivate = 4

## reading WST coeffs.
# number of wst coeffs
num_wst = 76
every_mean_wst = np.zeros((len(order_folders), num_wst))    # define empty total array
errbar_wst = np.zeros((len(order_folders), num_wst))        # define error bar array
if bool_rsd:
    every_mean_wst_rsd = np.zeros((len(order_folders), num_wst))
    errbar_wst_rsd = np.zeros((len(order_folders), num_wst))

for FC in range(len(files_to_read_wst)):
    every_wst = []

    # identify cosmology model
    cosmo = cosmo_parser(files_to_read_wst[FC])
    
    # `index`: give always the same index to the models
    index = order_folders[cosmo]
    
    # read .hdf5 file
    with h5py.File(root_WST+files_to_read_wst[FC], 'r') as fh:
        for key in fh.keys():
            dat = np.array(fh[key])
            every_wst.append(dat)
    
    if 'real' in files_to_read_wst[FC]:
        every_mean_wst[index] = np.mean(every_wst, axis=0)
        errbar_wst[index] = np.std(every_wst, axis=0)
    elif 'rsd' in files_to_read_wst[FC]:
        every_mean_wst_rsd[index] = np.mean(every_wst, axis=0)
        errbar_wst_rsd[index] = np.std(every_wst, axis=0)
    else:
        assert False, 'Undefined space in '+files_to_read_wst[FC]

print(f'Progress:  |/------|  {round(1/21*100, 2)}%', end='\r')

## reading Pk values
# it's assumed you calculate three multipoles, if not comment/delete the unnecessary lines...
# -> AND modify the `3*` factor in the rsd derivate of power spectrum

# number of Pk coeffs
num_pk = 96
every_mean_m_pk = np.zeros((len(order_folders), num_pk)) # monopole
every_mean_q_pk = np.zeros((len(order_folders), num_pk)) # quadrupole
every_mean_h_pk = np.zeros((len(order_folders), num_pk)) # hexadecapole
k = []                                                   # expected to be the same for all realizations
sn = np.zeros(len(order_folders))          # shot noise for veriation realizations
errbar_m_pk = np.zeros((len(order_folders), num_pk))
errbar_q_pk = np.zeros((len(order_folders), num_pk))
if bool_rsd:
    every_mean_m_pk_rsd = np.zeros((len(order_folders), num_pk))
    every_mean_q_pk_rsd = np.zeros((len(order_folders), num_pk))
    every_mean_h_pk_rsd = np.zeros((len(order_folders), num_pk))
    # k_rsd = []
    errbar_m_pk_rsd = np.zeros((len(order_folders), num_pk))
    errbar_q_pk_rsd = np.zeros((len(order_folders), num_pk))
    sn_rsd = np.zeros(len(order_folders))

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
    if 'real' in currfile:
        every_mean_m_pk[index] = np.mean(every_power_m, axis=0)
        every_mean_q_pk[index] = np.mean(every_power_q, axis=0)
        every_mean_h_pk[index] = np.mean(every_power_h, axis=0)
        sn = np.mean(every_sn, axis=0)
        errbar_m_pk[index] = np.std(every_power_m, axis=0)
        errbar_q_pk[index] = np.std(every_power_q, axis=0)
    elif 'rsd' in currfile:
        every_mean_m_pk_rsd[index] = np.mean(every_power_m, axis=0)
        every_mean_q_pk_rsd[index] = np.mean(every_power_q, axis=0)
        every_mean_h_pk_rsd[index] = np.mean(every_power_h, axis=0)
        sn_rsd = np.mean(every_sn, axis=0)
        errbar_m_pk_rsd[index] = np.std(every_power_m, axis=0)
        errbar_q_pk_rsd[index] = np.std(every_power_q, axis=0)
    else: assert False, f'Undefined space in {currfile}'

print(f'Progress:  |\------|  {round(2/21*100, 2)}%', end='\r')

every_mean_poles_rsd = np.zeros((len(order_folders), 3*num_pk))
if bool_rsd:
    for i in order_folders:
        ind = order_folders[i]
        every_mean_poles_rsd[ind] = np.concatenate((every_mean_m_pk_rsd[ind], every_mean_q_pk_rsd[ind], every_mean_h_pk_rsd[ind]))

# remove shot noise
if not  np.shape(sn) == ():
    every_mean_m_pk = every_mean_m_pk - sn[:, None]
    if bool_rsd:
        every_mean_m_pk_rsd = every_mean_m_pk_rsd - sn_rsd[:, None]

print(f'Progress:  |+------|  {round(3/21*100, 2)}%', end='\r')

# ------------- DERIVATES -------------------------------------------------------------------------------

n_pars = len(cosmological_pars)

derivates_wst = np.zeros((n_pars, num_wst))
derivates_pk  = np.zeros((n_pars, num_pk))
# derivates_pk  = np.zeros((n_pars, 3*num_pk))
if bool_rsd:
    derivates_wst_rsd = np.zeros((n_pars, num_wst))
    derivates_pk_rsd  = np.zeros((n_pars, 3*num_pk))    # <- MODIFY if haven't 3 multip.

# in evaluating the derivates using Quijote, you must take into account two facts:
#   1) 'Mnu' has 3 options of derivate (3 orders)
#   2) 'Ob' has been in(de)cremented both once and twice; here I only have used the double variation
#        a) if you want/have to use only the single variation simply comment the last 'elif' and modify
#            the first 'if'
#        b) if you want/have to use both, I simply haven't take into account this possibility XD 

for i in cosmological_pars:
    if "Mnu" not in i and "Ob" not in i:
        ind = order_dimension[i]
        derivates_wst[ind]    = (every_mean_wst[order_folders[i+"_p"]] - every_mean_wst[order_folders[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
        derivates_pk[ind] = (every_mean_m_pk[order_folders[i+"_p"]] - every_mean_m_pk[order_folders[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
        if bool_rsd:
            derivates_wst_rsd[ind] = (every_mean_wst_rsd[order_folders[i+"_p"]] - every_mean_wst_rsd[order_folders[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
            derivates_pk_rsd[ind]  = (every_mean_poles_rsd[order_folders[i+"_p"]] - every_mean_poles_rsd[order_folders[i+"_m"]]) / (2 * VarCosmoPar['d_'+i] )
    elif "Mnu" in i:
        # # >>> ORDER 1 >>>
        if order_derivate == 1:
            derivates_wst[order_dimension['Mnu']]    = (every_mean_wst[order_folders["Mnu_p"]] - every_mean_wst[order_folders["zeldovich"]]) / (0.1)
            derivates_pk[order_dimension['Mnu']] = (every_mean_m_pk[order_folders["Mnu_p"]] - every_mean_m_pk[order_folders["zeldovich"]]) / (0.1)
            if bool_rsd:
                derivates_wst_rsd[order_dimension['Mnu']]    = (every_mean_wst_rsd[order_folders["Mnu_p"]] - every_mean_wst_rsd[order_folders["zeldovich"]]) / (0.1)
                derivates_pk_rsd[order_dimension['Mnu']] = (every_mean_poles_rsd[order_folders["Mnu_p"]] - every_mean_poles_rsd[order_folders["zeldovich"]]) / (0.1)
        # >>> ORDER 2 >>>
        elif order_derivate == 2:
            derivates_wst[order_dimension['Mnu']] = (-every_mean_wst[order_folders['Mnu_pp']] + 4 * every_mean_wst[order_folders['Mnu_p']] - 3 * every_mean_wst[order_folders['zeldovich']]) / (2 * 0.1)
            derivates_pk[order_dimension['Mnu']] = (-every_mean_m_pk[order_folders['Mnu_pp']] + 4 * every_mean_m_pk[order_folders["Mnu_p"]] - 3 * every_mean_m_pk[order_folders['zeldovich']]) / (2 * 0.1)
            if bool_rsd:
                derivates_wst_rsd[order_dimension['Mnu']] = (-every_mean_wst_rsd[order_folders['Mnu_pp']] + 4 * every_mean_wst_rsd[order_folders['Mnu_p']] - 3 * every_mean_wst_rsd[order_folders['zeldovich']]) / (2 * 0.1)
                derivates_pk_rsd[order_dimension['Mnu']] = (-every_mean_poles_rsd[order_folders['Mnu_pp']] + 4 * every_mean_poles_rsd[order_folders["Mnu_p"]] - 3 * every_mean_poles_rsd[order_folders['zeldovich']]) / (2 * 0.1)
        # # >>> ORDER 4 >>>
        elif order_derivate == 4:
            derivates_wst[order_dimension['Mnu']] = (every_mean_wst[order_folders['Mnu_ppp']] - 12 * every_mean_wst[order_folders['Mnu_pp']] + 32 * every_mean_wst[order_folders["Mnu_p"]] - 21 * every_mean_wst[order_folders['zeldovich']]) / (12 * 0.1)
            derivates_pk[order_dimension['Mnu']] = (every_mean_m_pk[order_folders['Mnu_ppp']] - 12 * every_mean_m_pk[order_folders['Mnu_pp']] + 32 * every_mean_m_pk[order_folders["Mnu_p"]] - 21 * every_mean_m_pk[order_folders['zeldovich']]) / (12 * 0.1)
            if bool_rsd:
                derivates_wst_rsd[order_dimension['Mnu']] = (every_mean_wst_rsd[order_folders['Mnu_ppp']] - 12 * every_mean_wst_rsd[order_folders['Mnu_pp']] + 32 * every_mean_wst_rsd[order_folders["Mnu_p"]] - 21 * every_mean_wst_rsd[order_folders['zeldovich']]) / (12 * 0.1)
                derivates_pk_rsd[order_dimension['Mnu']] = (every_mean_poles_rsd[order_folders['Mnu_ppp']] - 12 * every_mean_poles_rsd[order_folders['Mnu_pp']] + 32 * every_mean_poles_rsd[order_folders["Mnu_p"]] - 21 * every_mean_poles_rsd[order_folders['zeldovich']]) / (12 * 0.1)
        else:
            assert False, 'Error in Mnu derivate: must choose order in [1; 2; 4]'
    elif "Ob" in i:
        derivates_wst[order_dimension['Ob']]    = (every_mean_wst[order_folders[i+"2_p"]] - every_mean_wst[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])
        derivates_pk[order_dimension['Ob']] = (every_mean_m_pk[order_folders[i+"2_p"]] - every_mean_m_pk[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])
        if bool_rsd:
            derivates_wst_rsd[order_dimension['Ob']]    = (every_mean_wst_rsd[order_folders[i+"2_p"]] - every_mean_wst_rsd[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])
            derivates_pk_rsd[order_dimension['Ob']] = (every_mean_poles_rsd[order_folders[i+"2_p"]] - every_mean_poles_rsd[order_folders[i+"2_m"]]) / (2*VarCosmoPar['d_'+i+"2"])

print(f'Progress:  |+/-----|  {round(4/21*100, 2)}%', end='\r')

# ---------------- COEFFICIENTS COVARIANCE MATRIX ----------------------------------------------------------------------------------------

# initializing empty fiducial arrays
fiducial_wst = []
fiducial_m_pk = []
fiducial_q_pk = []
fiducial_h_pk = []
k, sn_f = [], []

## Read WST files
with h5py.File(root_WST+'fiducial_wst_real.hdf5', 'r') as fh:
    for key in fh.keys():
        dat = np.array(fh[key])
        fiducial_wst.append(dat)

print(f'Progress:  |+\-----|  {round(5/21*100, 2)}%', end='\r')

## Read Pk files
#with h5py.File(root_Pk+'fiducial_pk_real.hdf5', 'r') as f:
with h5py.File(root_Pk+'fiducial_pk_nbk_real.hdf5', 'r') as f:
    for key in f.keys():
        dat = np.array(f[key])
        if   'mono'  in key: fiducial_m_pk.append(dat)
        elif 'quad'  in key: fiducial_q_pk.append(dat)
        elif 'hexa'  in key: fiducial_h_pk.append(dat)
        # elif 'k_coord' in key: k.append(dat)
        elif 'nbk_k' in key: k.append(dat)
        elif 'sn'    in key: sn_f.append(dat)
        else: assert False, fr'keyword {key} not found in `if` code'

print(f'Progress:  |++-----|  {round(6/21*100, 2)}%', end='\r')

# if you choose a minor number of fiducial realizations, here the diluition is executed
if max_fids < len(fiducial_m_pk):
    new_indices = randind(len(fiducial_m_pk, max_fids))
    fiducial_wst = fiducial_wst[new_indices,:]
    fiducial_m_pk = fiducial_m_pk[new_indices,:]
    fiducial_q_pk = fiducial_q_pk[new_indices,:]
    fiducial_h_pk = fiducial_h_pk[new_indices,:]
    k = k[new_indices,:]
    if not  np.shape(sn_f) == ():
        sn_f = sn_f[new_indices,:]

print(f'Progress:  |++/----|  {round(7/21*100, 2)}%', end='\r')

# remove shot noise
if not  np.shape(sn) == ():
    fiducial_m_pk = fiducial_m_pk - sn_f[:, None]

# transpose arrays
fiducial_wst = np.array(fiducial_wst).transpose()
fiducial_m_pk = np.array(fiducial_m_pk).transpose()

print(f'Progress:  |++\----|  {round(8/21*100, 2)}%', end='\r')

# covariance matices
cov_wst = np.cov(fiducial_wst)
print(f'Progress:  |++*----|  {round(8.5/21*100, 2)}%', end='\r')
cov_pk = np.cov(fiducial_m_pk)
    
print(f'Progress:  |+++----|  {round(9/21*100, 2)}%', end='\r')

# invese matrices with Hartlap formalism
Hart_wst = Hartlap(cov_wst, max_fids)
print(f'Progress:  |+++^---|  {round(9.5/21*100, 2)}%', end='\r')
Hart_pk  = Hartlap(cov_pk, max_fids)

print(f'Progress:  |+++/---|  {round(10/21*100, 2)}%', end='\r')

# initialize empty FIMs
Fish_wst, Fish_pk = np.zeros((n_pars,n_pars)), np.zeros((n_pars,n_pars))
# create FIMs (approx form)
for a in range(6): 
    for b in range(6):
        Fish_wst[a, b] = np.sum(np.dot(derivates_wst[a], np.dot(Hart_wst, np.transpose(derivates_wst[b]))))
        Fish_pk[a, b] = np.sum(np.dot(derivates_pk[a], np.dot(Hart_pk, np.transpose(derivates_pk[b]))))

print(f'Progress:  |+++\---|  {round(11/21*100, 2)}%', end='\r')

# get the inverse FIM aka (co)variances
inverse_wst = np.linalg.inv(Fish_wst)
print(f'Progress:  |+++*---|  {round(11.5/21*100, 2)}%', end='\r')
inverse_pk = np.linalg.inv(Fish_pk)

print(f'Progress:  |++++---|  {round(12/21*100, 2)}%', end='\r')

# write results on file
with h5py.File('./inv_FIM_est.hdf5', 'a') as wf:
    dataset_wst_real = wf.create_dataset('inv_wst_real', data = inverse_wst)
    dataset_pk_real = wf.create_dataset('inv_pk_real', data = inverse_pk)

print(f'Progress:  |++++/--|  {round(13/21*100, 2)}%', end='\r')

# the same as before, but about redshift space
if bool_rsd:
    fiducial_wst_rsd = []
    fiducial_m_pk_rsd = []
    fiducial_q_pk_rsd = []
    fiducial_h_pk_rsd = []
    k_rsd, sn_f_rsd = [], []

    with h5py.File(root_WST+'fiducial_wst_rsd.hdf5', 'r') as fh:
        for key in fh.keys():
            dat = np.array(fh[key])
            fiducial_wst_rsd.append(dat)

    print(f'Progress:  |++++\--|  {round(14/21*100, 2)}%', end='\r')

    with h5py.File(root_Pk+'fiducial_pk_nbk_rsd.hdf5', 'r') as f:
        for key in f.keys():
            dat = np.array(f[key])
            if   'mono' in key: fiducial_m_pk_rsd.append(dat)
            elif 'quad' in key: fiducial_q_pk_rsd.append(dat)
            elif 'hexa' in key: fiducial_h_pk_rsd.append(dat)
            elif 'nbk_k' in key: k.append(dat)
            elif 'sn'   in key: sn_f_rsd.append(dat)
            else: assert False, f'keyword {k} not found in `if` code'

    print(f'Progress:  |+++++--|  {round(15/21*100, 2)}%', end='\r')

    if max_fids < len(fiducial_wst_rsd):
        new_indices = randind(len(fiducial_wst_rsd, max_fids))
        print(f'Progress:  |+++++^-|  {round(15.5/21*100, 2)}%', end='\r')
        fiducial_wst = fiducial_wst[new_indices,:]
        fiducial_m_pk_rsd = fiducial_m_pk_rsd[new_indices] 
        fiducial_q_pk_rsd = fiducial_q_pk_rsd[new_indices]
        fiducial_h_pk_rsd = fiducial_h_pk_rsd[new_indices] 
        k_rsd = k_rsd[new_indices]
        if not  np.shape(sn_f_rsd) == ():
            sn_f_rsd = sn_f_rsd[new_indices]

    print(f'Progress:  |+++++/-|  {round(16/21*100, 2)}%', end='\r')

    # remove shot noise
    if not  np.shape(sn) == ():
        fiducial_m_pk_rsd = fiducial_m_pk_rsd - sn_f_rsd[:, None]

    fiducial_poles_rsd = np.concatenate((fiducial_m_pk_rsd, fiducial_q_pk_rsd, fiducial_h_pk_rsd), axis=1)

    # transpose arrays
    fiducial_wst_rsd = np.array(fiducial_wst_rsd).transpose()
    fiducial_poles_rsd = np.array(fiducial_poles_rsd).transpose()

    cov_wst_rsd = np.cov(fiducial_wst_rsd)
    print(f'Progress:  |+++++!-|  {round(16.5/21*100, 2)}%', end='\r')
    cov_pk_rsd = np.cov(fiducial_poles_rsd)

    print(f'Progress:  |+++++\-|  {round(17/21*100, 2)}%', end='\r')

    Hart_wst_rsd = Hartlap(cov_wst_rsd, max_fids)
    print(f'Progress:  |+++++*-|  {round(17.5/21*100, 2)}%', end='\r')
    Hart_pk_rsd = Hartlap(cov_pk_rsd, max_fids)

    print(f'Progress:  |++++++-|  {round(18/21*100, 2)}%', end='\r')

    Fish_wst_rsd, Fish_pk_rsd = np.zeros((n_pars,n_pars)), np.zeros((n_pars,n_pars))

    for a in range(6): 
        for b in range(6):
            Fish_wst_rsd[a, b] = np.sum(np.dot(derivates_wst_rsd[a], np.dot(Hart_wst_rsd, np.transpose(derivates_wst_rsd[b]))))
            Fish_pk_rsd[a, b] = np.sum(np.dot(derivates_pk_rsd[a], np.dot(Hart_pk_rsd, np.transpose(derivates_pk_rsd[b]))))

    print(f'Progress:  |++++++/|  {round(19/21*100, 2)}%', end='\r')

    inverse_wst_rsd = np.linalg.inv(Fish_wst_rsd)
    print(f'Progress:  |++++++!|  {round(19.5/21*100, 2)}%', end='\r')
    inverse_pk_rsd = np.linalg.inv(Fish_pk_rsd)

    print(f'Progress:  |++++++\|  {round(20/21*100, 2)}%', end='\r')


    # write results on file
    with h5py.File('./inv_FIM_est.hdf5', 'a') as wf:
        dataset_wst_rsd = wf.create_dataset('inv_wst_rsd', data = inverse_wst_rsd)
        dataset_pk_rsd  = wf.create_dataset('inv_pk_rsd',  data = inverse_pk_rsd)

    print(f'Progress:  |+++++++|  {round(21/21*100, 2)}%\nEND\n')