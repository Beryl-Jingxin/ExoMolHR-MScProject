# Scientific Computing Individual Research Project

Scientific Computing Individual Research Project for PHAS0077 in UCL during the 2019/2020 academic season


## Project Title

The title of this project is:
A relational database to store measured high-resolution spectra.


## Purpose

The aim of this project is to extract molecule line lists of high-resolution spectra from ExoMol website and convert ExoMol format data into HITRAN format data.


## Environment

```bash
physics cluster `prima`
Python 3.6.8
Linux
```


## Install Dependancies

If the packages used in codes are not installed, please use the following command in terminal.

For example: install pandas

```bash
[jzhang@prima ~]$ python3 -V
Python 3.6.8
[jzhang@prima ~]$ which python
/user/bin/python
[jzhang@prima ~]$ python3.6 -m pip install --user pandas
```

## Usage Instructions

Download codes.

```bash
git clone https://github.com/Beryl-Jingxin/PHAS0077.git
```

Part of this project build instructions:

```
├── project
│     └── code
│          ├── exomol_master.py   
│          ├── AlH.py                    
│          ├── C2.py                        
│          ├── C2H2.py                  
│          ├── CO2.py                   
│          ├── H2O.py
│          ├── H3O_p.py                    
│          ├── NH3.py                        
│          ├── TiO.py                                  
│          ├── hitran.py
│          ├── ipynb_with_sample_output
│          │     ├── exomol_master.ipynb
│          │     ├── AlH.ipynb
│          │     ├── C2.ipynb
│          │     ├── C2.ipynb
│          │     ├── C2H2.ipynb
│          │     ├── CO2.ipynb
│          │     ├── H2O.ipynb
│          │     ├── H3O_p.ipynb
│          │     ├── NH3.ipynb
│          │     ├── TiO.ipynb
│          │     └── hitran.ipynb
│          └── data
│                └── url
│                     └── api__urls.txt
└── README.md
```

Run the codes in Linux:

```bash
chmod u+x xxx.py
python3 xxx.py
```

To run the codes of project, there are three steps:
Access in code folder:

###  Step 1    Download Files

```bash
[jzhang@prima code]$ chmod u+x exomol_master.py
[jzhang@prima code]$ python3 exomol_master.py             // Get API URLs and download def files
[jzhang@prima code]$ cd data
[jzhang@prima data]$ wget -b -r -i ./url/api__urls.txt    // Download states and trans files in backgrounder
```

###  Step 2    Process Data

For example, if we process the data of NH3:

```bash
[jzhang@prima data]$ cd ..
[jzhang@prima code]$ chmod u+x NH3.py
[jzhang@prima code]$ python3 NH3.py
```

Or use `nohup` to run program in backgrounder.

```bash
[jzhang@prima data]$ mkdir log
[jzhang@prima data]$ cd ..
[jzhang@prima code]$ chmod u+x NH3.py
[jzhang@prima code]$ nohup python3 NH3.py > ./data/log/NH3.log 2>&1 &
```

Check the course of NH3.py (this command can be typed anywhere in Linux):

```bash
[jzhang@prima ~]$ ps -aux|grep NH3.py    // NH3.py can be replaced as its PID
```

###  Step 3    Save Data

Save data as CSV and HITRAN format. This process is fast so that do not need to run this code in backgrounder.

```bash
[jzhang@prima code]$ chmod u+x hitran.py
[jzhang@prima code]$ python3 hitran.py
```

Compress the results:

```bash
[jzhang@prima data]$ tar -jcvf result.tar.bz2 result
```

## Final structure

Part of the final structure:

```
├── project
│     └── code
│          ├── exomol_master.py   
│          ├── AlH.py                    
│          ├── C2.py                        
│          ├── C2H2.py                  
│          ├── CO2.py                   
│          ├── H2O.py
│          ├── H3O_p.py                    
│          ├── NH3.py                        
│          ├── TiO.py                                  
│          ├── hitran.py
│          ├── ipynb_with_sample_output
│          │     ├── exomol_master.ipynb
│          │     ├── AlH.ipynb
│          │     ├── C2.ipynb
│          │     ├── C2H2.ipynb
│          │     ├── CO2.ipynb
│          │     ├── H2O.ipynb
│          │     ├── H3O_p.ipynb
│          │     ├── NH3.ipynb
│          │     ├── TiO.ipynb
│          │     └── hitran.ipynb
│          └── data
│                ├── def
│                │    ├── AlCl_(27Al)(35Cl)_27Al-35Cl__MoLLIST.def
│                │    ├── AlCl_(27Al)(37Cl)_27Al-37Cl__MoLLIST.def
│                │    ├── AlF_(27Al)(19F)_27Al-19F__MoLLIST.def
│                │    └── AlH_(26Al)(1H)_26Al-1H__AlHambra.def
│                ├── url
│                │    └── api__urls.txt
│                ├── wget-log
│                ├── www.exomol.com
│                │    └── db
│                │         ├── NH3    
│                │         │    └── 14N-1H3
│                │         │          └── CoYuTe
│                │         │                ├── 14N-1H3__CoYuTe__00000-00100.trans.bz2
│                │         │                ├── 14N-1H3__CoYuTe__00100-00200.trans.bz2
│                │         │                ├── 14N-1H3__CoYuTe.states.bz2
│                │         │                └── NH3__14N-1H3__CoYuTe_pf.csv
│                │         ├── H2O    
│                │         │    └── 1H2-16O
│                │         │          └── POKAZATEL
│                │         │                ├── 1H2-16O__POKAZATEL__00000-00100.trans.bz2
│                │         │                ├── 1H2-16O__POKAZATEL__00100-00200.trans.bz2
│                │         │                ├── 1H2-16O__POKAZATEL.states.bz2
│                │         │                ├── 1H2-16O__POKAZATEL_v1.states.bz2
│                │         │                └── H2O__1H2-16O__POKAZATEL_pf.csv
│                │         └── TiO    
│                │              ├── 46Ti-16O
│                │              │      └── Toto
│                │              │           ├── 46Ti-16O__Toto.trans.bz2
│                │              │           ├── 46Ti-16O__Toto.states.bz2
│                │              │           └── TiO__46Ti-16O__Toto_pf.csv
│                │              └── 48Ti-16O
│                │                    └── Toto
│                │                          ├── 48Ti-16O__Toto.trans.bz2
│                │                          ├── 48Ti-16O__Toto.states.bz2
│                │                          └── TiO__48Ti-16O__Toto_pf.csv
│                ├── log
│                │    ├── AlH.log
│                │    ├── C2.log
│                │    ├── C2H2.log
│                │    ├── CO2.log
│                │    ├── H2O.log
│                │    ├── H3O_p.log
│                │    ├── NH3.log
│                │    └── TiO.log
│                ├── result
│                │    ├── AlH__27Al-1H__AlHambra.csv
│                │    ├── C2__12C2__8states.csv
│                │    ├── C2H2__12C2-1H2__aCeTY.csv
│                │    ├── CO2__12C-16O2__UCL-4000.csv
│                │    ├── H2O__1H2-16O__POKAZATEL.csv
│                │    ├── H3O_p__1H3-16O_p__eXeL.csv
│                │    ├── NH3__14N-1H3__CoYuTe.csv
│                │    ├── TiO__48Ti-16O__Toto.csv
│                │    └── one_file
│                │           ├── csv_format.csv        // CSV format result
│                │           ├── demo_hitran.txt       // Not the result
│                │           └── HITRAN_format.txt     // HITRAN format result
│                └── result.tar.bz2     
└── README.md
```
