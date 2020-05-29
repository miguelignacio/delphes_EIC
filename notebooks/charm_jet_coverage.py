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

parser.add_argument("-n", "--input", type=str,
                    help="Directory containing input files")
parser.add_argument("-x", "--xvar", type=str, default='jet_p',
                    help="jet_pt, jet_p, etc.")

args = parser.parse_args()

pt_name = "Jet.PT"
eta_name = "Jet.Eta"
flavor_name = "Jet.Flavor"

if args.xvar == "genjet_p":
    pt_name = "GenJet.PT"
    eta_name = "GenJet.Eta"
    flavor_name = "GenJet.Flavor"
    
branchlist=[pt_name, eta_name, flavor_name]

print("Loading data...")

df = eat.UprootLoad([f"../{args.input}/*/out.root"], "Delphes", branches=branchlist)

#df = df[:10000]

n_gen = len(df)
print(f"n_gen = {n_gen}")


# Code by D. Higinbotham
#========================
# Example Contour Plot

# Bins in angle and momentum

jet_pt = np.concatenate(df[pt_name].to_numpy()).ravel()
jet_eta = np.concatenate(df[eta_name].to_numpy()).ravel()
#jet_tag = np.concatenate(df['GenJet.BTag'].to_numpy()).ravel()
jet_flavor = np.concatenate(df[flavor_name].to_numpy()).ravel()
jet_theta = 2*np.arctan(np.exp(-jet_eta))
jet_p = jet_pt*np.cosh(jet_eta)

#jet_tagged = (jet_tag == 1)
charm_flavor = ( jet_flavor == 4 )



angles = np.radians(np.linspace(0, 180, 90))

mom = np.linspace(0,100,10)
xlabel = "Jet $p_T$ [GeV]"
xvals = jet_pt[charm_flavor]
thetavals = jet_theta[charm_flavor]

if args.xvar == "jet_pt":
    mom=np.linspace(0,50,10)
    xlabel = "Charm Jet $p_T$ [GeV]"
    xvals = jet_pt[charm_flavor]
elif args.xvar == "jet_p":
    mom=np.linspace(0,80,16)
    xlabel = "Charm Jet Momentum [GeV]"
    xvals = jet_p[charm_flavor]
elif args.xvar == "genjet_p":
    mom=np.linspace(0,80,16)
    xlabel = "Generator-Level Charm Jet Momentum [GeV]"
    xvals = jet_p[charm_flavor]
else:
    print("Unknown x variable")
    sys.exit()


values, thetaedges, redges = np.histogram2d(thetavals, xvals, bins=[angles, mom])
r, theta = np.meshgrid( redges[:-1], thetaedges[:-1])


# Make the plot
fig,ax = plt.subplots(subplot_kw=dict(projection='polar'),dpi=300)
#ax.contourf(theta, r, values,levels=18)
ax.contourf(theta, r, values)
list=[0,np.pi/6,np.pi/3,np.pi/2,4*np.pi/6,5*np.pi/6,np.pi]
ax.set_xticks(list)
ax.set_thetamin(0)
ax.set_thetamax(180)
plt.xlabel(xlabel,labelpad=-40,fontsize=18)
plt.title(f"CC-DIS, 10GeVx275GeV, $Q^2>100\\mathrm{{GeV^2}}$", fontsize=24)
plt.text(3.45/2,np.max(mom)+12.5,'Degrees',fontsize=18,multialignment='right')
ax.tick_params(axis='x', labelsize=12 )
ax.tick_params(axis='y', labelsize=12 )
plt.tight_layout()

plt.savefig(f"charm_jet_coverage_{args.xvar}_{args.input}.png")
plt.savefig(f"charm_jet_coverage_{args.xvar}_{args.input}.pdf")
