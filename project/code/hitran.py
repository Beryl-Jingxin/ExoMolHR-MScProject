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


# Part 1: HITRAN Online Information.
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
    path_mol_iso.append(line.split('/')[-4] + '/' + line.split('/')[-3] + '/' + line.split('/')[-2])

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
unc_formula = pd.DataFrame(eval(str(iso_slug_list).replace('1H','H').replace('-','').replace('_p','+')))
unc_formula.columns = ['exomol formula']

# Information for calculations to obtain the HITRAN format data.
hitran_online = pd.DataFrame()
hitran_online['exomol formula'] = unc_formula['exomol formula']
hitran_online['molecule ID'] = ['50','26','51','2','1','52','11','53','53','53','53','53']
hitran_online['isotopologue ID'] = ['1','1','1','1','1','1','1','1','2','3','4','5']
hitran_online['fractional abundance'] = ['1','0.977599','1','0.984204','0.997317','1','0.995872','1','1','1','1','1']


# # Part 2: Process Data to Satisfy HITRAN Format
# Convert uncertainty values which we calculated in each molecule codes into uncertainty code.
# See information from https://hitran.org/docs/uncertainties/.
def convert_uncertainty_code(HITRAN_df):
    HITRAN_num = HITRAN_df['Ierr'].count()
    uncertainty_code = 0
    Ierr = []
    for i in range(HITRAN_num):
        uncertainty = HITRAN_df['Ierr'].values[i]
        uncertainty_value = float(uncertainty)
        if (0.1 <= uncertainty_value < 1):
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
        elif (0.00000001 <= uncertainty_value < 0.0000001):
            uncertainty_code = '{:>1}'.format(8) + '40000'
        elif (uncertainty_value < 0.00000001):
            uncertainty_code = '{:>1}'.format(9) + '40000'
        else:
            uncertainty_code = '{:>1}'.format(0) + '40000'
        Ierr.append(uncertainty_code)
    return Ierr


# Convert CSV format into HITRAN format. All intensities less than 1.0E-30 can be ignored.
# We only extract those rows whose intensity is larger than 1.0E-30.
# To save data as HITRAN format, we just use _ to instead of blanks and save as a demo result.
# Then we will use this demo result to convert it into HITRAN format result.
def convert_csv_to_HITRAN(csv_df):
    HITRAN_df = csv_df[csv_df['S'] > 1.0E-30]
    Ierr = convert_uncertainty_code(HITRAN_df)

    HITRAN_df['M'] = HITRAN_df['M'].map('{:_>2}'.format)
    HITRAN_df['I'] = HITRAN_df['I'].map('{:>1}'.format)
    HITRAN_df['v'] = HITRAN_df['v']
    HITRAN_df['S'] = HITRAN_df['S'] * fractional_abundance
    HITRAN_df['S'] = HITRAN_df['S'].map('{:_>10.3E}'.format)
    HITRAN_df['A'] = HITRAN_df['A'].map('{:_>10.3E}'.format)
    HITRAN_df['gm_a'] = '_' * 5
    HITRAN_df['gm_s'] = '_' * 5
    HITRAN_df['E_f'] = HITRAN_df['E_f'].map('{:_>10.4F}'.format)
    HITRAN_df['n_a'] = '_' * 4
    HITRAN_df['dt_a'] = '_' * 8
    HITRAN_df['V_i'] = HITRAN_df['V_i']
    HITRAN_df['V_f'] = HITRAN_df['V_f']
    HITRAN_df['Q_i'] = HITRAN_df['Q_i']
    HITRAN_df['Q_f'] = HITRAN_df['Q_f']
    HITRAN_df['Ierr'] = Ierr
    HITRAN_df['Iref'] = '_' * 12
    HITRAN_df['*'] = '_'
    HITRAN_df['g_i'] = HITRAN_df['g_i'].map('{:_>7.1F}'.format)
    HITRAN_df['g_f'] = HITRAN_df['g_f'].map('{:_>7.1F}'.format)

    return HITRAN_df


# Read data in chunks.
def read_txt_in_chunks(path, chunk_size=1024*1024):
    file = open(path, 'r')
    while True:
        chunk_data = file.read(chunk_size)
        if not chunk_data:
            break
        yield chunk_data


def save_hitran(df, demo_path, save_path):
    '''
    Sort HITRAN format data by increasing wavenumbers and then convert wavenumbers format
    to be similar as other columns which is using _ instead of blank.
    Since if there is " in a value of DataFrame, then when we save data into a text file,
    there will be two more " at the begining and the end of this value.
    To avoid this problem, we replace ' to be upper and replace " to be lower.
    We will then convert upper and lower back into ' and " later when we write data
    into a HITRAN format result text file.
    Change the order to Change the column order to satisfy the HITRAN format.
    Save a demo file as a text file for convering into HITRAN format result text file later.
    Read the demo file and replace string upper and lower back to ' and ".
    Then write data into a text file. After all we obtain the HITRAN format result text file.
    '''
    # Sort HITRAN format data by increasing wavenumbers.
    df = df.sort_values(['v'], ascending = True).reset_index(drop=True)
    # Convert wavenumbers format to be similar as other columns which is using _ instead of blank.
    df['v'] = df['v'].map('{:_>12.6F}'.format)
    # To avoid the changes of ' and " when converting from csv to txt.
    df['V_i'] = df['V_i'].str.replace("'","upper")
    df['V_i'] = df['V_i'].str.replace('"','lower')
    df['V_f'] = df['V_f'].str.replace("'","upper")
    df['V_f'] = df['V_f'].str.replace('"','lower')
    df['Q_i'] = df['Q_i'].str.replace("'","upper")
    df['Q_i'] = df['Q_i'].str.replace('"','lower')
    df['Q_f'] = df['Q_f'].str.replace("'","upper")
    df['Q_f'] = df['Q_f'].str.replace('"','lower')
    # Change the column order to satisfy the HITRAN format.
    order = ['M', 'I', 'v', 'S', 'A', 'gm_a', 'gm_s', 'E_f', 'n_a', 'dt_a',
             'V_i', 'V_f', 'Q_i', 'Q_f', 'Ierr', 'Iref', '*', 'g_i', 'g_f']
    df = df[order]
    # Save a demo file for converting into HITRAN format.
    df.to_csv(demo_path, header=None, index=False)

    with open(save_path, 'w') as save_file:
        for chunk in read_txt_in_chunks(demo_path):
            # Replace back to ' and " to satisfy the HITRAN format.
            string = str(chunk).replace(',','').replace('_',' ').replace("upper","'").replace('lower','"')
            save_file.write(string)


# A folder for save demo files, the files which are named with demo are not the results.
# Create a folder for saving demo files.
# If the folder exists, save files directory,otherwise, create it.
demo_file_path = './data/result/demo/'
if os.path.exists(demo_file_path):            # Determine whether the folder exists or not.
    pass
else:
    os.makedirs(demo_file_path, exist_ok=True)   # Create the folder.


# A folder for save HITRAN format files, the files which are named as
# molecule__iso-slug__isotopologue_HITRAN.txt
# Create a folder for saving HITRAN format files.
# If the folder exists, save files directory,otherwise, create it
hitran_file_path = './data/result/hitran/'
if os.path.exists(hitran_file_path):           # Determine whether the folder exists or not.
    pass
else:
    os.makedirs(hitran_file_path, exist_ok=True)  # Create the folder.


'''
Save HITRAN format files

Read all molecule CSV format results.
Concatenate these CSV files into a large CSV file.
In this CSV file, data are sorted by wavenumbers and grouped by different
molecule__iso-slug__isotopologue names.

Convert this large CSV format result into HITRAN format.
Thre are columns names in CSV format result DataFrame,
however, we need to calculate intensity * fractional abundance.
Therefore, we remove those rows which have column names.
We just use the column isotopologue number (I),
extract those rows whose isotopologue numbers are not the string I.
Then we reset the type of the none empty column values.
After these, we convert CSV format into HITRAN format.

Save HITRAN format files in a loop and save the filenames as
molecule__iso-slug__isotopologue_HITRAN.txt
'''

# Concatenate HITRAN dataframes for saving as one final HITRAN format text with the whole data.
col_name = ['M', 'I', 'v', 'S', 'A', 'gamma_air', 'gamma_self',
            'E_f', 'n_air', 'delta_air', 'V_i', 'V_f', 'Q_i', 'Q_f',
            'Ierr', 'Iref', '*', 'g_i', 'g_f']

df = dict()
one_csv_df = pd.DataFrame()
HITRAN_df = pd.DataFrame()
csv_filenames = glob.glob('./data/result/csv/' + '*csv')
for csv_filename in csv_filenames:
    df[csv_filename] = pd.read_csv(csv_filename, header=None, names=col_name,
                                   chunksize=100_000_000, iterator=True, low_memory=False)
    formula = csv_filename.replace('_p','+').split('_')[2].replace('1H','H').replace('-','')
    fractional_abundance = float(hitran_online[hitran_online['exomol formula'].isin([formula])]['fractional abundance'].values)
    for chunk in df[csv_filename]:
        # Concatenate CSV files.
        one_csv_df = one_csv_df.append(chunk)

        # For converting HITRAN format.
        csv_chunk = chunk
        # Remove the rows which has column names.
        csv_df = csv_chunk[~csv_chunk['I'].isin(['I'])]
        # Reset type of each column values.
        csv_df[['M', 'I']] = csv_df[['M', 'I']].astype(int)
        csv_df[['v', 'S', 'A', 'E_f','Ierr', 'g_i', 'g_f']] = csv_df[['v', 'S', 'A', 'E_f','Ierr', 'g_i', 'g_f']].astype(float)
        # Convert CSV format into HITRAN format.
        hitran_df = convert_csv_to_HITRAN(csv_df)

        # Save as HITRAN format per species
        hitran_dfs = hitran_df
        save_hitran_filename = csv_filename.replace('\\','/').split('/')[-1].split('.')[0]
        demo_paths = demo_file_path + save_hitran_filename + '_demo_hitran.txt'
        save_paths = hitran_file_path + save_hitran_filename + '_hitran.txt'
        save_hitrans = save_hitran(hitran_df, demo_paths, save_paths)

        # For saving all data into one HITRAN format file.
        HITRAN_df = HITRAN_df.append(hitran_df)


# 2.1 Save CSV Format Result
# This csv format result includes whole data (all molecules) and is sorted by increasing wavenumbers.
# Save the large concatenated CSV format result file.
# Create a folder for saving result files.
# If the folder exists, save files directory,otherwise, create it.
two_files_path = './data/result/two_files/'
if os.path.exists(two_files_path):               # Determine whether the folder exists or not.
    pass
else:
    os.makedirs(two_files_path, exist_ok=True)   # Create the folder.

save_csv = one_csv_df.to_csv(two_files_path + 'CSV_format.csv', header=True, index=False)


# 2.2 Save HITRAN Format Result
# This HITRAN format result includes whole data (all molecules) and is sorted by increasing wavenumbers.
demo_path = demo_file_path + 'demo_hitran.txt'
save_path = two_files_path + 'HITRAN_format.txt'
save_HITRAN = save_hitran(HITRAN_df, demo_path, save_path)
