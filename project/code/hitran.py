


import os
import math
import glob
import requests
import numpy as np
import pandas as pd
from io import StringIO





molecule = []
iso_slug = []
isotopologue = []
path_mol_iso = []
for line in open('./data/url/api__urls.txt'):
    molecule.append(line.split('/')[-4])
    iso_slug.append(line.split('/')[-3])
    isotopologue.append(line.split('/')[-2])
    path_mol_iso.append(line.split('/')[-4] + '/' + line.split('/')[-3] + '/' + line.split('/')[-2])

molecule_list = list(set(molecule))
molecule_list.sort(key=molecule.index)

iso_slug_list = list(set(iso_slug))
iso_slug_list.sort(key=iso_slug.index)

isotopologue_list = list(set(isotopologue))
isotopologue_list.sort(key=isotopologue.index)

path_mol_iso_list = list(set(path_mol_iso))
path_mol_iso_list.sort(key=path_mol_iso.index)



unc_formula = pd.DataFrame(eval(str(iso_slug_list).replace('1H','H').replace('-','').replace('_p','+')))
unc_formula.columns = ['exomol formula']




hitran_online = pd.DataFrame()
hitran_online['exomol formula'] = unc_formula['exomol formula']
hitran_online['molecule ID'] = ['50','26','51','2','1','52','11','53','53','53','53','53']
hitran_online['isotopologue ID'] = ['1','1','1','1','1','1','1','1','2','3','4','5']
hitran_online['fractional abundance'] = ['1','0.977599','1','0.984204','0.997317','1','0.995872','0.2','0.2','0.2','0.2','0.2']


def convert_uncertainty_code(HITRAN_df):
    HITRAN_num = HITRAN_df['Ierr'].count()
    uncertainty_code = 0
    Ierr = []
    for i in range(HITRAN_num):
        uncertainty = HITRAN_df['Ierr'].values[i]
        uncertainty_value = float(uncertainty)
        if (1 <= uncertainty_value):
            uncertainty_code = '{:>1}'.format(0) + '40000'
        elif (0.1 <= uncertainty_value < 1):
            uncertainty_code = '{:>1}'.format(1) + '40000'
        elif (0.01 <= uncertainty_value < 0.1):
            uncertainty_code = '{:>1}'.format(2) + '40000'
        elif (0.001 <= uncertainty_value < 0.01):
            uncertainty_code = '{:>1}'.format(3) + '40000'
        elif (0.0001 <= uncertainty_value < 0.001):
            uncertainty_code = '{:>1}'.format(4) + '40000'
        elif (0.00001 <= uncertainty_value < 0.0001):
            uncertainty_code = '{:>1}'.format(5) + '40000'
        elif (0.000001 <= uncertainty_value < 0.00001):
            uncertainty_code = '{:>1}'.format(6) + '40000'
        elif (0.0000001 <= uncertainty_value < 0.000001):
            uncertainty_code = '{:>1}'.format(7) + '40000'
        elif (uncertainty_value < 0.0000001):
            uncertainty_code = '{:>1}'.format(8) + '40000'
        Ierr.append(uncertainty_code)
    return Ierr




def convert_csv_to_HITRAN(csv_df):
    HITRAN_df = csv_df[csv_df.S > 1.0E-30]
    Ierr = convert_uncertainty_code(HITRAN_df)
    
    HITRAN_df['M'] = HITRAN_df.M.map('{:_>2}'.format)
    HITRAN_df['I'] = HITRAN_df.I.map('{:>1}'.format)
    HITRAN_df['v'] = HITRAN_df.v
    HITRAN_df['S'] = HITRAN_df.S * fractional_abundance
    HITRAN_df['S'] = HITRAN_df.S.map('{:_>10.3E}'.format)
    HITRAN_df['A'] = HITRAN_df.A.map('{:_>10.3E}'.format)
    HITRAN_df['gm_a'] = '_' * 5
    HITRAN_df['gm_s'] = '_' * 5
    HITRAN_df['E_f'] = HITRAN_df.E_f.map('{:_>10.4F}'.format)
    HITRAN_df['n_a'] = '_' * 4
    HITRAN_df['dt_a'] = '_' * 8
    HITRAN_df['V_i'] = HITRAN_df.V_i
    HITRAN_df['V_f'] = HITRAN_df.V_f
    HITRAN_df['Q_i'] = HITRAN_df.Q_i
    HITRAN_df['Q_f'] = HITRAN_df.Q_f
    HITRAN_df['Ierr'] = Ierr
    HITRAN_df['Iref'] = '_' * 12
    HITRAN_df['*'] = '_'
    HITRAN_df['g_i'] = HITRAN_df.g_i.map('{:_>7.1F}'.format)
    HITRAN_df['g_f'] = HITRAN_df.g_f.map('{:_>7.1F}'.format)

    return HITRAN_df




col_name = ['M', 'I', 'v', 'S', 'A', 'gamma_air', 'gamma_self',
            'E_f', 'n_air', 'delta_air', 'V_i', 'V_f', 'Q_i', 'Q_f',
            'Ierr', 'Iref', '*', 'g_i', 'g_f']

df = dict()
one_csv_df = pd.DataFrame()
HITRAN_df = pd.DataFrame()
csv_filenames = glob.glob('./data/result/' + '*csv')
for csv_filename in csv_filenames:
    df[csv_filename] = pd.read_csv(csv_filename, header=None, names=col_name,
                                   chunksize=100_000_000, iterator=True, low_memory=False)
    formula = csv_filename.replace('_p','+').split('_')[2].replace('1H','H').replace('-','')
    fractional_abundance = float(hitran_online[hitran_online['exomol formula'].isin([formula])]['fractional abundance'].values)
    for chunk in df[csv_filename]:
        one_csv_df = one_csv_df.append(chunk)
        csv_chunk = chunk
        csv_df = csv_chunk[~csv_chunk['I'].isin(['I'])]
        csv_df[['M', 'I']] = csv_df[['M', 'I']].astype(int)
        csv_df[['v', 'S', 'A', 'E_f','Ierr', 'g_i', 'g_f']] = csv_df[['v', 'S', 'A', 'E_f','Ierr', 'g_i', 'g_f']].astype(float)
        hitran_df = convert_csv_to_HITRAN(csv_df)
        HITRAN_df = HITRAN_df.append(hitran_df)


# Create a folder for saving result files.
# If the folder exists, save files directory,otherwise, create it.
one_file_path = './data/result/one_file/'
if os.path.exists(one_file_path):                   # Determine whether the folder exists or not.
    pass
else:
    os.makedirs(one_file_path, exist_ok=True)       # Create the folder.
    
save_csv = one_csv_df.to_csv(one_file_path + 'csv_format.csv', header=True, index=False)



HITRAN_df = HITRAN_df.sort_values(['v'], ascending = True).reset_index(drop=True)
HITRAN_df['v'] = HITRAN_df['v'].map('{:_>12.6F}'.format)




HITRAN_df['V_i'] = HITRAN_df['V_i'].str.replace("'","upper")
HITRAN_df['V_i'] = HITRAN_df['V_i'].str.replace('"','lower')
HITRAN_df['V_f'] = HITRAN_df['V_f'].str.replace("'","upper")
HITRAN_df['V_f'] = HITRAN_df['V_f'].str.replace('"','lower')
HITRAN_df['Q_i'] = HITRAN_df['Q_i'].str.replace("'","upper")
HITRAN_df['Q_i'] = HITRAN_df['Q_i'].str.replace('"','lower')
HITRAN_df['Q_f'] = HITRAN_df['Q_f'].str.replace("'","upper")
HITRAN_df['Q_f'] = HITRAN_df['Q_f'].str.replace('"','lower')





order = ['M', 'I', 'v', 'S', 'A', 'gm_a', 'gm_s', 'E_f', 'n_a', 'dt_a',
         'V_i', 'V_f', 'Q_i', 'Q_f', 'Ierr', 'Iref', '*', 'g_i', 'g_f']
HITRAN_df = HITRAN_df[order]



HITRAN_df.to_csv(one_file_path + 'demo_hitran.txt', header=None, index=False)




def read_txt_in_chunks(path, chunk_size=1024*1024):
    file = open(path, 'r')
    while True:
        chunk_data = file.read(chunk_size)
        if not chunk_data:
            break
        yield chunk_data

HITRAN_path = one_file_path + 'demo_hitran.txt'
with open(one_file_path + 'HITRAN_format.txt', 'w') as save_file:
    for chunk in read_txt_in_chunks(HITRAN_path):
        string = str(chunk).replace(',','').replace('_',' ').replace("upper","'").replace('lower','"')
        save_file.write(string)

