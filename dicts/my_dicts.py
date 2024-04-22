COSMOPAR = {
#                   | Om   | Ob   |   h   |  n_s  | s_8 | Mnu | w |

    'fiducial' :    [0.3175, 0.049, 0.6711, 0.9624, 0.834, 0, -1],
    'zeldovich':    [0.3175, 0.049, 0.6711, 0.9624, 0.834, 0, -1],
    
    'Mnu_p' :       [0.3175, 0.049, 0.6711, 0.9624, 0.834, 0.1, -1],
    'Mnu_pp' :      [0.3175, 0.049, 0.6711, 0.9624, 0.834, 0.2, -1],
    'Mnu_ppp' :     [0.3175, 0.049, 0.6711, 0.9624, 0.834, 0.4, -1],
    
    'h_m' :         [0.3175, 0.049, 0.6511, 0.9624, 0.834, 0, -1],
    'h_p' :         [0.3175, 0.049, 0.6911, 0.9624, 0.834, 0, -1],
    
    'ns_m' :        [0.3175, 0.049, 0.6711, 0.9424, 0.834, 0, -1],
    'ns_p' :        [0.3175, 0.049, 0.6711, 0.9824, 0.834, 0, -1],
    
    'Ob_m' :        [0.3175, 0.048, 0.6711, 0.9624, 0.834, 0, -1],
    'Ob_p' :        [0.3175, 0.050, 0.6711, 0.9624, 0.834, 0, -1],
    'Ob2_m' :       [0.3175, 0.047, 0.6711, 0.9624, 0.834, 0, -1],
    'Ob2_p' :       [0.3175, 0.051, 0.6711, 0.9624, 0.834, 0, -1],
    
    'Om_m' :        [0.3075, 0.049, 0.6711, 0.9624, 0.834, 0, -1],
    'Om_p' :        [0.3275, 0.049, 0.6711, 0.9624, 0.834, 0, -1],
    
    's8_m' :        [0.3175, 0.049, 0.6711, 0.9624, 0.819, 0, -1],
    's8_p' :        [0.3175, 0.049, 0.6711, 0.9624, 0.849, 0, -1]
}
"""Dictionary for the 15 cosmology models used"""

folders_ordered = {
    0 :'fiducial'  ,    
    1 :'h_m'       , 
    2 :'h_p'       ,    
    3 :'Mnu_p'     , 
    4 :'Mnu_pp'    ,
    5 :'Mnu_ppp'   ,
    6 :'ns_m'      , 
    7 :'ns_p'      ,    
    8 :'Ob2_m'     , 
    9 :'Ob2_p'     ,    
    10 :'Om_m'      , 
    11 :'Om_p'      ,
    12 :'s8_m'      , 
    13 :'s8_p'      ,
    14 :'zeldovich' ,
}
"""Dictionary to associate a number 0-14 to a cosmology model"""

order_folders = {
    'fiducial'  : 0,
    'h_m'       : 1,  'h_p'       : 2,
    'Mnu_p'     : 3,  'Mnu_pp'    : 4, 'Mnu_ppp'   : 5, 
    'ns_m'      : 6,  'ns_p'      : 7,
    'Ob2_m'     : 8,  'Ob2_p'     : 9,
    'Om_m'      : 10, 'Om_p'      : 11,
    's8_m'      : 12, 's8_p'      : 13,
    'zeldovich' : 14
}
"""Dictionary to associate a cosmology model to a number 0-14"""

order_dimension = {
    'Om'  : 0, 'Om ' : 0,
    'Ob'  : 1, 'Ob ' : 1, 'Ob2' : 1,
    'h'   : 4, 'h  ' : 4, # 4
    'ns'  : 3, 'ns ' : 3,
    's8'  : 2, 's8 ' : 2, # 2
    'Mnu' : 5
}
"""Dictionary to associate a parameter to a number 0-5"""

cosmological_pars = {
    'Om'  : 0, 'Ob'  : 1, 'h'   : 4,
    'ns'  : 3, 's8'  : 2, 'Mnu' : 5
}
"""Slim version of `order_dimension`"""


VarCosmoPar = {
    'd_h'  : 0.02, 'd_ns' : 0.02,  'd_Ob' : 0.001, 'd_Ob2': 0.002,
    'd_Om' : 0.01, 'd_s8' : 0.015
}
"""Dictionary that gives the variation used to create the variant models"""

fiducial_vals = {
    'Ob'  : 0.049, 'Ob2' : 0.049,
    'Om'  : 0.3175,
    'h'   : 0.6711,
    'n_s' : 0.9624, 'ns'  : 0.9624,
    's_8' : 0.834,  's8'  : 0.834,
    'Mnu' : 0
}
"""Gives the fiducial value of the parameters"""

cov_scale ={
    0: 1,
    1: 1.52,
    2: 2.48,
    3: 3.44
}
"""Gives the alpha value in the ellipse/confidence level calculus"""