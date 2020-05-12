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

from concurrent.futures import ThreadPoolExecutor

plt.style.use(hep.style.ATLAS)

plt.rcParams.update({'font.sans-serif': "Arial",
                     'font.family': "sans-serif",
                     'font.size': 30,
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': 'Arial',
                     })


datasets = {}

datacache = uproot.cache.ArrayCache(100*1024)   # 100 kB

def UprootLoad(filestring):
    print(f"::UprootLoad({filestring}")
    executor = ThreadPoolExecutor(8)

    df_list = []
    data = uproot.pandas.iterate([filestring], 
                                 "Delphes", 
                                 branches=["GenJet.PT", "GenJet.Eta", "GenJet.Flavor"],
                                 flatten=False,
                                 executor=executor)
    for dataframe in data:
        df_list.append(dataframe)
    
    df = pd.concat(df_list)

    return df

print("Loading data...")

df_ct18nnlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18NNLO/[0-3]/out.root")
df_ct18nlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18NLO/[0-3]/out.root")
df_ct18Annlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18ANNLO/[0-3]/out.root")
df_ct18Anlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18ANLO/[0-3]/out.root")


print("Done loading data for the study!")


# Units
global u_fb, u_mb, xsection
u_fb = 1
u_mb = 1e12*u_fb
xsection = 2.408e-08*u_mb

def DifferentialXS(df, which='all'):
    global u_fb, u_mb, xsection, n_gen
    n_gen = len(df)
    print(f"n_gen = {n_gen}")
    jet_pt = np.concatenate(df['GenJet.PT'].to_numpy()).ravel()
    jet_eta = np.concatenate(df['GenJet.Eta'].to_numpy()).ravel()
    jet_flavor = np.concatenate(df['GenJet.Flavor'].to_numpy()).ravel()

    #jet_basics = (0.01 < jet_eta) & (jet_eta < 0.9)
    jet_basics = (-5.0 < jet_eta) & (jet_eta < 5.0)
    
    all_flavor = ( jet_flavor > -999.0 ) & (jet_basics)
    light_flavor = ( jet_flavor < 4 ) & (jet_basics)
    charm_flavor = ( jet_flavor == 4 ) & (jet_basics)

    selection = all_flavor
    if which == 'charm':
        selection = charm_flavor


    (counts, bins) = np.histogram(jet_pt[ selection ], range=(10,50), 
                                  bins=(10,12.5,15,20,25,30,40,50))
                                  #bins=40)

    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths/2

    errors = np.sqrt(counts)
    rel_errors = errors/counts


    # convert counts to dsigma/dpT * 100/fb
    dsigma_dPT = counts * xsection * 100*u_fb**(-1) / n_gen / bin_widths
    dsigma_dPT_errors = rel_errors * dsigma_dPT
    
    return (bin_centers, bin_widths, dsigma_dPT, dsigma_dPT_errors)

fig, ax = plt.subplots(figsize=(12,12))
plt.axis('off')


gridspec = fig.add_gridspec(ncols=1, nrows=2, width_ratios=[
    1], height_ratios=[3, 1])


# Log plot of differential cross-sections
ax1 = fig.add_subplot(gridspec[0, 0])

ax1.grid(which='both', axis='both')

(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nnlo, 'all')
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='mediumblue')

(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18Annlo, 'all')
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18ANNLO)', marker='s', ms=10, ls='none', linewidth=2, color='saddlebrown')

(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nlo, 'all')
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NLO)', marker='v', ms=10, ls='none', linewidth=2, color='magenta')

(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18Anlo, 'all')
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18ANLO)', marker='D', ms=10, ls='none', linewidth=2, color='coral')

charm_ct18nnlo = DifferentialXS(df_ct18nnlo, 'charm')
(bins, bin_widths, xs, xs_error) = charm_ct18nnlo
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, fillstyle='none', color='mediumblue')

charm_ct18Annlo = DifferentialXS(df_ct18Annlo, 'charm')
(bins, bin_widths, xs, xs_error) = charm_ct18Annlo
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18ANNLO)', marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')

charm_ct18nlo = DifferentialXS(df_ct18nlo, 'charm')
(bins, bin_widths, xs, xs_error) = charm_ct18nlo
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NLO)', marker='v', ms=10, ls='none', linewidth=2, fillstyle='none', color='magenta')

charm_ct18Anlo = DifferentialXS(df_ct18Anlo, 'charm')
(bins, bin_widths, xs, xs_error) = charm_ct18Anlo
ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18ANLO)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='coral')


plt.ylabel('$d\\sigma/dp_T \\times 100\\mathrm{fb^{-1}}$ [$\\mathrm{GeV^{-1}}$]')
plt.xlabel('Generated Jet $p_T$ [GeV]')
plt.yscale('log')

plt.title("CC-DIS, 18x275GeV, $Q^2>100\\mathrm{GeV^2}$", fontsize=20)

ax1.set_ylim([10, 2e5])

ax1.legend(fontsize=18)


# ratio plot for charm production
ax2 = fig.add_subplot(gridspec[1, 0])

ax2.grid(which='both', axis='both')

ratio_charm_ct18nnlo_ct18nnlo = charm_ct18nnlo[2]/charm_ct18nnlo[2]
ratio_charm_ct18Annlo_ct18nnlo = charm_ct18Annlo[2]/charm_ct18nnlo[2]
ratio_charm_ct18nlo_ct18nnlo = charm_ct18nlo[2]/charm_ct18nnlo[2]
ratio_charm_ct18Anlo_ct18nnlo = charm_ct18Anlo[2]/charm_ct18nnlo[2]

ratioerr_charm_ct18nnlo_ct18nnlo = ratio_charm_ct18nnlo_ct18nnlo*np.sqrt((charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
ratioerr_charm_ct18Annlo_ct18nnlo = ratio_charm_ct18Annlo_ct18nnlo*np.sqrt((charm_ct18Annlo[3]/charm_ct18Annlo[2])**2 + (charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
ratioerr_charm_ct18nlo_ct18nnlo = ratio_charm_ct18nlo_ct18nnlo*np.sqrt((charm_ct18nlo[3]/charm_ct18nlo[2])**2 + (charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
ratioerr_charm_ct18Anlo_ct18nnlo = ratio_charm_ct18Anlo_ct18nnlo*np.sqrt((charm_ct18Anlo[3]/charm_ct18Anlo[2])**2 + (charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)

print(ratio_charm_ct18Annlo_ct18nnlo)
print(ratio_charm_ct18nlo_ct18nnlo)
print(ratio_charm_ct18Anlo_ct18nnlo)

ax2.axhline(1.0, color='k', linestyle='-', linewidth=2)
#ax2.errorbar(bins, ratio_charm_ct18nnlo_ct18nnlo, xerr = bin_widths/2, yerr=ratioerr_charm_ct18nnlo_ct18nnlo, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')
errorboxes = [mpl.patches.Rectangle((x - xe, y - ye), 2*xe, 2*ye)
              for x, y, xe, ye in zip(bins, ratio_charm_ct18nnlo_ct18nnlo, (bin_widths/2), ratioerr_charm_ct18nnlo_ct18nnlo)]

# Create patch collection with specified colour/alpha
pc = mpl.collections.PatchCollection(errorboxes, facecolor='k', alpha=0.25)

# Add collection to axes
ax2.add_collection(pc)

ax2.errorbar(bins, ratio_charm_ct18Annlo_ct18nnlo, xerr = bin_widths/2, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')
ax2.errorbar(bins, ratio_charm_ct18nlo_ct18nnlo, xerr = bin_widths/2, marker='v', ms=10, ls='none', linewidth=2, fillstyle='none', color='magenta')
ax2.errorbar(bins, ratio_charm_ct18Anlo_ct18nnlo, xerr = bin_widths/2, marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='coral')

ax2.set_ylim([0.5, 1.5])
plt.ylabel('Ratio to\nCT18NNLO', fontsize=18)

plt.tight_layout()

plt.savefig("jet_differential_xs.png")
