#!/usr/bin/env python
# coding: utf-8

# Import all what we need.
import os
import math
import glob
import requests
import numpy as np
import pandas as pd
from tqdm import tqdm
from io import StringIO


# # Part 1: HITRAN Online Information.
# Get the names of molecules, iso-slugs and isotopoluge datasets from the api__urls.txt
# which saved the URLs with molecule, iso-slug and isotopologue.
# Combine them with '/' for reading files from folders more convenient later.
molecule = []
iso_slug = []
isotopologue = []
path_mol_iso = []
for line in open('./data/url/api__urls.txt'):
    molecule.append(line.split('/')[-4])
    iso_slug.append(line.split('/')[-3])
    isotopologue.append(line.split('/')[-2])
    path_mol_iso.append(line.split('/')[-4] + '/' + line.split('/')[-3]
                        + '/' + line.split('/')[-2])

molecule_list = list(set(molecule))
molecule_list.sort(key=molecule.index)
iso_slug_list = list(set(iso_slug))
iso_slug_list.sort(key=iso_slug.index)
isotopologue_list = list(set(isotopologue))
isotopologue_list.sort(key=isotopologue.index)
path_mol_iso_list = list(set(path_mol_iso))
path_mol_iso_list.sort(key=path_mol_iso.index)

# Convert the iso-slug names into the ones which are shown in the table of
# HITRAN online website. It will help us to get their corresponding molecule
# numbers, isotopologue numbers and fractional abundances. 
# The HITRAN online URL is: https://hitran.org/docs/iso-meta/.
unc_formula = pd.DataFrame(eval(str(iso_slug_list).replace('1H','H')
                                .replace('-','').replace('_p','+')))
unc_formula.columns = ['exomol formula']

# CO2 information for HITRAN format.
path_mol_iso_list = ['CO2/12C-16O2/UCL-4000']
hitran_online = pd.DataFrame()
hitran_online['molecule ID'] = ['2']
hitran_online['isotopologue ID'] = ['1']
hitran_online['exomol formula'] = ['12C16O2']
hitran_online['fractional abundance'] = ['0.984204']
hitran_online['Q(296K)'] = ['286.09']

path_mol_iso = path_mol_iso_list[0]
M_mol_iso = 'M of ' + path_mol_iso.replace('/','__')
molecule_id = int(hitran_online['molecule ID'][0])
isotopologue_id = int(hitran_online['isotopologue ID'][0])
fractional_abundance = float(hitran_online['fractional abundance'][0])

read_path = './data/www.exomol.com/db/'

# Create a folder for saving result files.
# If the folder exists, save files directory,otherwise, create it.
result_path = './data/result/csv/'
# Determine whether the folder exists or not.
if os.path.exists(result_path):
    pass
else:
    # Create the folder.
    os.makedirs(result_path, exist_ok=True)


# Part 2: Process Data.
# 2.1 Read States File.
# Consider column names of states file with def files.
path_mol_iso = path_mol_iso_list[0]
states_col_name = (['ID','energy','g_tot','J_tot','Unc']
                   + ['Gtot','e/f','n1','n2','l2','n3','W',
                      'm1','m2','l','m3','r','v1','v2','v3'])

# Read compressed states file in chunks directly.
# Extract rows of states file whose uncertainty indices are small than 0.001
# to be the considered states file.
s_df = dict()
states_df = pd.DataFrame()
states_filenames = glob.glob(read_path + path_mol_iso + '/' + path_mol_iso.split('/')[1]
                             + '__' + path_mol_iso.split('/')[2] + '.states.bz2')

for states_filename in states_filenames:
    s_df[states_filename] = pd.read_csv(states_filename, compression='bz2', sep='\s+',
                                        header=None, names=states_col_name,
                                        chunksize=100_000_000, iterator=True,
                                        low_memory=False)
    for chunk in s_df[states_filename]:
        states_df = states_df.append(chunk)
        
# Extract rows of states file whose uncertainty indices are small than 0.001.
unc_states_df = states_df[states_df['Unc'] < float(0.001)]

# 2.2 Read Partition Function File.
pf_col_name = ['T', 'Q']
pf_path = read_path + path_mol_iso + '/' + path_mol_iso.replace('/','__') + '_pf.csv'
pf_url = ('http://www.exomol.com/db/' + path_mol_iso + '/'
          + path_mol_iso.split('/')[1] + '__' + path_mol_iso.split('/')[2] + '.pf')   
response = requests.get(pf_url)
content = response.text  
pf_data = pd.read_csv(StringIO(content), sep='\s+', names=pf_col_name, header=None, engine='python')
pf_data.to_csv(pf_path, header=False)
pf_df = pd.read_csv(pf_path, header=None, names=pf_col_name)
# Partition function defined as sum over states at standard 296K
Q = pf_df.iloc[296-1]['Q']

# 2.3 Read Transitions Files.
# Extract rows of transitionos files whose upper states ID and
# lower states ID are all in considered states file.
def extract_unc_trans(trans_df):
    upper_id = trans_df['i'].values
    lower_id = trans_df['f'].values
    state_id = unc_states_df['ID'].values
    
    # Extract the same upper states ID from states_df
    unc_trans_i_df = pd.DataFrame()
    for id_1 in tqdm(state_id):
        unc_trans_i_df = unc_trans_i_df.append(trans_df[trans_df['i'].isin([id_1])])
        
    # Extract the same lower states ID from states_df
    unc_trans_df = pd.DataFrame()
    for id_2 in tqdm(state_id):
        unc_trans_df = unc_trans_df.append(unc_trans_i_df[unc_trans_i_df['f'].isin([id_2])])
        
    return unc_trans_df

# 2.4 Calculating.
# HITRAN Parameters for calculating.
T = 296                      # Reference temperature k
h = 6.62607015e-34               # Planck's const (J s)
c = 299792458                  # Velocity of light (m s-1)
kB = 1.380649e-23                # Boltzmann's const (J K-1)
c2 = h * c * 100 / kB             # Second radiation constant (cm K)
pi_c_8 = 1 / (8 * np.pi * c * 100)    # 8 * pi * c (cm-1 s)
c2_T = c2 / T                  # c2 / T  (cm)

# Process data for CSV format.
def calculate_csv(unc_states_df, unc_trans_df):
    unc_upper_id = unc_trans_df['i'].values
    unc_lower_id = unc_trans_df['f'].values 
    state_id = unc_states_df['ID']
    unc_trans_num = unc_trans_df['i'].count()

    wavenumber = []                      # Vacuum wavenumber (cm−1)
    intensity = pd.DataFrame()               # Intensities (cm-1/molecule cm-2) at standard 296K         
    A_coefficient = []                    # Einstein A-coefficient
    lower_state_energy = pd.DataFrame()         # lower state energy
    uncertainty = []                      # Uncertainty indices
    weight_upper_state = pd.DataFrame()         # Statistical weight of upper state
    weight_lower_state = pd.DataFrame()         # Statistical weight of lower state
    upper_global_quanta = []                # Upper-state 'global' quanta
    lower_global_quanta = []                # Lower-state 'global' quanta
    upper_local_quanta = []                 # Upper-state 'local' quanta
    lower_local_quanta = []                 # Lower-state 'local' quanta

    for i in tqdm(range(unc_trans_num)):
        id_i = unc_upper_id[i]
        id_f = unc_lower_id[i]
        A = unc_trans_df['A_if'].values[i]                       # Einstein-A coefficient (s−1)
        g_i = unc_states_df[state_id.isin([id_i])]['g_tot'].values       # Total degeneracy of upper state
        g_f = unc_states_df[state_id.isin([id_f])]['g_tot'].values       # Total degeneracy of lower state
        E_i = unc_states_df[state_id.isin([id_i])]['energy'].values      # Upper state energy
        E_f = unc_states_df[state_id.isin([id_f])]['energy'].values      # Lower state energy
        unc_i = unc_states_df[state_id.isin([id_i])]['Unc'].values       # Uncertainty indices of upper state
        unc_f = unc_states_df[state_id.isin([id_f])]['Unc'].values       # Uncertainty indices of lower state
        
        u_m1 = unc_states_df[state_id.isin([id_i])]['m1'].values[0]
        u_m2 = unc_states_df[state_id.isin([id_i])]['m2'].values[0]
        u_m3 = unc_states_df[state_id.isin([id_i])]['m3'].values[0]
        u_l = unc_states_df[state_id.isin([id_i])]['l'].values[0]
        u_r = unc_states_df[state_id.isin([id_i])]['r'].values[0]
        u_J = unc_states_df[state_id.isin([id_i])]['J_tot'].values[0]
        u_ef = unc_states_df[state_id.isin([id_i])]['e/f'].values[0]

        l_m1 = unc_states_df[state_id.isin([id_f])]['m1'].values[0]
        l_m2 = unc_states_df[state_id.isin([id_f])]['m2'].values[0]
        l_m3 = unc_states_df[state_id.isin([id_f])]['m3'].values[0]
        l_l = unc_states_df[state_id.isin([id_f])]['l'].values[0]
        l_r = unc_states_df[state_id.isin([id_f])]['r'].values[0]
        l_J = unc_states_df[state_id.isin([id_f])]['J_tot'].values[0]
        l_ef = unc_states_df[state_id.isin([id_f])]['e/f'].values[0]

        V_i = '     %2d%2d%2d%2d%2d' % (u_m1,u_m2,u_l,u_m3,u_r) + ','    # Upper-state 'global' quanta
        V_f = '     %2d%2d%2d%2d%2d' % (l_m1,l_m2,l_l,l_m3,l_r) + ','    # Lower-state 'global' quanta
        Q_i = '          %3d%2s' % (u_J,u_ef) + ','               # Upper-state 'local' quanta
        Q_f = '          %3d%2s' % (l_J,l_ef) + ','               # Lower-state 'local' quanta

        unc = math.sqrt(unc_i ** 2 + unc_f ** 2)                   # Uncertainty idices
        v = float(abs(E_i - E_f))                             # Vacuum wavenumber (cm−1)
        S = g_i * A * np.exp(- c2_T * E_f) * (1 - np.exp(- c2_T * v)) * pi_c_8 / (v ** 2) / Q    # Intensities

        wavenumber.append(v)
        intensity = intensity.append(pd.DataFrame(S))
        A_coefficient.append(A)
        lower_state_energy = lower_state_energy.append(pd.DataFrame(E_f))
        uncertainty.append(unc)
        weight_upper_state = weight_upper_state.append(pd.DataFrame(g_i))
        weight_lower_state = weight_lower_state.append(pd.DataFrame(g_f))
        upper_global_quanta += V_i.split(',')
        lower_global_quanta += V_f.split(',')
        upper_local_quanta += Q_i.split(',')
        lower_local_quanta += Q_f.split(',')
    
    iso_csv_df = pd.DataFrame()
    iso_csv_df['v'] = wavenumber               # Vacuum wavenumber (cm−1)
    iso_csv_df['S'] = intensity.values           # Intensities (cm-1/molecule cm-2) at standard 296K  
    iso_csv_df['A'] = A_coefficient             # Einstein A-coefficient
    iso_csv_df['E_f'] = lower_state_energy.values    # Lower state energy
    iso_csv_df['Ierr'] = uncertainty             # Uncertainty indices
    iso_csv_df['g_i'] = weight_upper_state.values    # Statistical weight of upper state
    iso_csv_df['g_f'] = weight_lower_state.values    # Statistical weight of lower state
    iso_csv_df[M_mol_iso] = molecule_id           # Molecule number
    iso_csv_df['I'] = isotopologue_id            # Isotopologue number
    iso_csv_df['gm_a'] = np.nan                # Air-broadened half-width
    iso_csv_df['gm_s'] = np.nan                # Self-broadened half-width
    iso_csv_df['n_a'] = np.nan                 # Temperature-dependence exponent for gamma_air
    iso_csv_df['dt_a'] = np.nan                # Air pressure-induced line shift
    iso_csv_df['V_i'] = list(filter(None, upper_global_quanta))      # Upper-state 'global' quanta
    iso_csv_df['V_f'] = list(filter(None, lower_global_quanta))      # Lower-state 'global' quanta
    iso_csv_df['Q_i'] = list(filter(None, upper_local_quanta))       # Upper-state 'local' quanta
    iso_csv_df['Q_f'] = list(filter(None, lower_local_quanta))       # Lower-state 'local' quanta                  
    iso_csv_df['Iref'] = np.nan                 # Reference indices
    iso_csv_df['*'] = np.nan                   # Flag

    return iso_csv_df
    

# Part 3: Save as CSV Format.
trans_col_name = ['i', 'f', 'A_if']
t_df = dict()
species_csv_df = pd.DataFrame()
trans_filenames = glob.glob(read_path + path_mol_iso + '/' + '*trans.bz2')
for trans_filename in tqdm(trans_filenames):
    t_df[trans_filename] = pd.read_csv(trans_filename, compression='bz2', sep='\s+',
                                       usecols=[0,1,2], header=None, names=trans_col_name,
                                       chunksize=100_000_000, iterator=True, low_memory=False)
    # Set an empty DataFrame to avoid meeting empty considered transitions files.
    iso_csv_df = pd.DataFrame()
    for trans_df in t_df[trans_filename]:
        unc_trans_df = extract_unc_trans(trans_df)
        if len(unc_trans_df) != 0:
            iso_csv_df = calculate_csv(unc_states_df, unc_trans_df)
            
    species_csv_df = species_csv_df.append(iso_csv_df)

order = [M_mol_iso, 'I', 'v', 'S', 'A', 'gm_a', 'gm_s', 'E_f', 'n_a', 'dt_a',
         'V_i', 'V_f', 'Q_i', 'Q_f', 'Ierr', 'Iref', '*', 'g_i', 'g_f']
species_csv_df = species_csv_df[order]
# Sort by increasing wavenumber
species_csv_df = species_csv_df.sort_values(['v'], ascending = True).reset_index(drop=True)
# Save into a CSV file with column names which contain molecule name.
species_csv_df.to_csv(result_path + path_mol_iso.replace('/','__') + '.csv', header=True, index=False) 
