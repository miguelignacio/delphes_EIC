import matplotlib as mpl
import uproot
import matplotlib.pyplot as plt
import scipy
import numpy as np
import math
import pandas as pd
import seaborn as sns
import mplhep as hep
#import zfit
import inspect
import sys
import argparse 
import glob

from concurrent.futures import ThreadPoolExecutor

plt.style.use(hep.style.ATLAS)

plt.rcParams.update({'font.sans-serif': "Arial",
                     'font.family': "sans-serif",
                     'font.size': 30,
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': 'Arial',
                     })


import EICAnalysisTools as eat


# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", type=str,
                    help="Main input subfolder")

args = parser.parse_args()


# Grab all the CSV files
csvfiles = glob.glob(f"{args.input}/*/tagging_study.csv")

print(csvfiles)


dataframes = []

for csvfile in csvfiles:
    dataframes.append(pd.read_csv(csvfile))


data = pd.concat(dataframes)

data = data.groupby(["Variation"], as_index=False).sum()

print(data.head(10))

# Add the Punzi DOM

DiscoverySignificance = 5.0
xsection = eat.TotalXSection('CC_DIS_e10_p275_CT18NNLO')
lumi = 100 #/fb
n_gen = len(dataframes)*2e5 # 200,000 collisions generated per file

allLight = float(data[data["Variation"] == "all" ]["Light"])
allCharm = float(data[data["Variation"] == "all" ]["Charm"])


data["LightEff"] = data["Light"]/allLight
data["Light_100fb"] = xsection*lumi*data["Light"]/n_gen
data["CharmEff"] = data["Charm"]/allCharm
data["Charm_100fb"] = xsection*lumi*data["Charm"]/n_gen




data["PunziFOM"] = (data["Charm"]/allCharm)/(DiscoverySignificance/2.0 + np.sqrt(data["Light_100fb"]))
data["CharmErrFOM"] =  data["Charm_100fb"]/np.sqrt(data["Charm_100fb"] + 2.0*data["Light_100fb"])

pd.set_option('display.max_rows', data.shape[0]+1)

print(data)


# Find the row with the max(PunziFOM)

#FOM_choice = "PunziFOM"
FOM_choice = "CharmErrFOM"

# Remove that "all" row before optimization printing
data = data[data["Variation"] != "all"]

optimal_row = data[FOM_choice].idxmax()

print(data.iloc[optimal_row])
print(f"Number of generated events: {n_gen}")

# Hold light jet efficiency ~constant at some value and find the best operating point there.

#print(data[ np.abs(data["LightEff"] - 2.0e-5)/1.0e-5 < 0.5 ])

#lighteff_target = 4.00e-5
target = float(data.iloc[optimal_row][FOM_choice])
tolerance = 2

print(f"Optimization of charm efficiency holding {FOM_choice} at ~{target}")
print(data[ np.abs(data[FOM_choice] - target) < tolerance ].head(20))
optimal_row = data[ np.abs(data[FOM_choice] - target) < tolerance ][FOM_choice].idxmax()

print(data.iloc[optimal_row])

