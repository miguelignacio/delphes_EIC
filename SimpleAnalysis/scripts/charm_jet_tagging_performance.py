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
import awkward

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

parser.add_argument("-d", "--dir", type=str,
                    help="Directory containing input files")
parser.add_argument("-i", "--input", type=str,
                    help="Main input subfolder")
parser.add_argument("-x", "--xvar", type=str, default='pt',
                    help="pt, eta, bjorken_x")

args = parser.parse_args()



branchlist=["*"]

print("Loading data...")

df = eat.UprootLoad([f"{args.dir}/{args.input}/0/out.root"], "tree", branches=branchlist)


def DrawDiffTagEfficiencyPlot(df, draw_config={}):
    xvar=draw_config['xvar']
    xrange=draw_config['xrange']
    xbins=draw_config['xbins']
    ylimits=draw_config['ylimits']
    xlimits=draw_config['xlimits']
    yunits=draw_config['yunits']
    xunits=draw_config['xunits']


    fig, ax = plt.subplots(figsize=(12,8))
    plt.axis('off')
    
    
    gridspec = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
    
    
    # Log plot of differential cross-sections
    ax1 = fig.add_subplot(gridspec[0, 0])
    
    ax1.grid(which='both', axis='both')
    
    (bins, bin_widths, eff, eff_error) = eat.DifferentialTaggingEfficiency(df, x=xvar, xrange=xrange, xbins=xbins, which='light')
    ax1.errorbar(bins, eff, xerr = bin_widths/2, yerr=eff_error, label='light jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='red')
    
    
    (bins, bin_widths, eff, eff_error) = eat.DifferentialTaggingEfficiency(df, x=xvar, xrange=xrange, xbins=xbins, which='charm')
    ax1.errorbar(bins, eff, xerr = bin_widths/2, yerr=eff_error, label='charm jets (CT18NNLO)', marker='D', ms=10, ls='none', linewidth=2, color='blue')
    
    
    
    xvar_symbol="p_T"
    if xvar.find('Eta') != -1:
        xvar_symbol='\\eta'
        
    plt.ylabel(f'$\\varepsilon_{{tag}}$ {yunits}')
    plt.xlabel(f'Reconstructed Jet ${xvar_symbol}$ {xunits}')
    
    beamconfig = "10x275"

    plt.title(f"CC-DIS, {beamconfig}GeV, $Q^2>100\\mathrm{{GeV^2}}$", fontsize=20)
    
    ax1.set_ylim(ylimits)
    ax1.set_xlim(xlimits)
    
    plt.xlabel(draw_config['xlabel'])

    if xvar == 'bjorken_x' or xvar == 'JB_x':
        plt.xscale('log')

    ax1.legend(fontsize=18)
    
    plt.yscale('log')
    
    plt.tight_layout()
    
    plt.savefig(f"jet_differential_tagefficiency_{xvar}_{args.input}.png")
    plt.savefig(f"jet_differential_tagefficiency_{xvar}_{args.input}.pdf")


def DrawDiffTagYieldPlot(df, draw_config={}, process='CC_DIS_e10_p275_CT18NNLO'):
    xvar=draw_config['xvar']
    xrange=draw_config['xrange']
    xbins=draw_config['xbins']
    ylimits=draw_config['ylimits']
    xlimits=draw_config['xlimits']
    yunits=draw_config['yunits']
    xunits=draw_config['xunits']
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.axis('off')
    
    
    gridspec = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
    
    
    # Log plot of differential cross-sections
    ax1 = fig.add_subplot(gridspec[0, 0])
    
    ax1.grid(which='both', axis='both')
    
    (bins, bin_widths, xs, xs_error) = eat.DifferentialTaggingYield(df, x=xvar, xrange=xrange, xbins=xbins, which='light', process=process, target_lumi=100)
    ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='light jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='red')
    
    
    charm_ct18nnlo = eat.DifferentialTaggingYield(df, x=xvar, xrange=xrange, xbins=xbins, which='charm', process=process, target_lumi=100)
    (bins, bin_widths, xs, xs_error) = charm_ct18nnlo
    ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='blue')
    
    
    
    xvar_symbol="p_T"
    if xvar.find('Eta') != -1:
        xvar_symbol='\\eta'
        
    plt.ylabel(f'$\\varepsilon_{{tag}} \\times \\sigma_{{\\mathrm{{CC-DIS}}}} \\times 100\\mathrm{{fb^{{-1}}}}$')
    plt.xlabel(f'Reconstructed Jet ${xvar_symbol}$ {xunits}')
    plt.yscale('log')
    
    beamconfig = "10x275"
    process_name = 'CC-DIS'

    if process.find('NC') != -1:
        process_name = 'NC-DIS'
    plt.title(f"{process_name}, {beamconfig}GeV, $Q^2>100\\mathrm{{GeV^2}}$", fontsize=20)
    
    ax1.set_ylim(ylimits)
    ax1.set_xlim(xlimits)

    plt.xlabel(draw_config['xlabel'])

    if xvar == 'bjorken_x' or xvar == 'JB_x':
        plt.xscale('log')


    ax1.legend(fontsize=18)
    
    
    plt.tight_layout()
    
    plt.savefig(f"jet_differential_tagyield_{xvar}_{args.input}.png")
    plt.savefig(f"jet_differential_tagyield_{xvar}_{args.input}.pdf")


# Inputs
draw_config={}

if args.xvar == 'eta':
    draw_config['xvar'] = 'jet_eta'
    draw_config['xrange'] = [-4.0, 4.0]
    draw_config['xbins'] = [-4.0, -3.5, -2.5, -1.0, 1.0, 2.5, 3.5, 4.0]
    draw_config['ylimits'] = [1e-5, 1]
    draw_config['xlimits'] = [-5,5]
    draw_config['yunits'] = ''
    draw_config['xunits'] = ''
    draw_config['xlabel'] = 'Reconstructed Jet $\eta$'
elif args.xvar == 'pt':
    draw_config['xvar'] = 'jet_pt'
    draw_config['xrange'] = [10,50]
    draw_config['xbins'] = [10,12.5,15,20,25,35,60]
    draw_config['ylimits'] = [1e-5, 1]
    draw_config['xlimits'] = [0,60]
    draw_config['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
    draw_config['xunits'] = '[GeV]'
    draw_config['xlabel'] = 'Reconstructed Jet $p_T$ [GeV]'
elif args.xvar == 'bjorken_x':
    draw_config['xvar'] = 'bjorken_x'
    draw_config['xrange'] = [5e-3,1]
    draw_config['xbins'] = np.concatenate([np.arange(1e-2,1e-1,3.33333e-2),np.arange(0.1,0.5,0.2),[0.5,1.0]])
    #draw_config['xbins'] = [5e-3, 1e-2, 1e-1, 1]
    draw_config['ylimits'] = [1e-5, 1]
    draw_config['xlimits'] = [5e-3,1]
    draw_config['yunits'] = ''
    draw_config['xunits'] = ''
    draw_config['xlabel'] = 'Bjorken x'
elif args.xvar == 'JB_x':
    draw_config['xvar'] = 'JB_x'
    draw_config['xrange'] = [5e-3,1]
    draw_config['xbins'] = np.concatenate([np.arange(1e-2,1e-1,3.33333e-2),np.arange(0.1,0.5,0.2),[0.5,1.0]])
    #draw_config['xbins'] = [5e-3, 1e-2, 1e-1, 1]
    draw_config['ylimits'] = [1e-5, 1]
    draw_config['xlimits'] = [5e-3,1]
    draw_config['yunits'] = ''
    draw_config['xunits'] = ''
    draw_config['xlabel'] = 'Reconstructed $x_{JB}$'


DrawDiffTagEfficiencyPlot(df,draw_config)

draw_config={}
if args.xvar == 'eta':
    draw_config['xvar'] = 'jet_eta'
    draw_config['xrange'] = [-4.0, 4.0]
    draw_config['xbins'] = [-4.0, -3.5, -2.5, -1.0, 1.0, 2.5, 3.5, 4.0]
    draw_config['ylimits'] = [1, 1e3]
    draw_config['xlimits'] = [-5,5]
    draw_config['yunits'] = ''
    draw_config['xunits'] = ''
    draw_config['xlabel'] = 'Reconstructed Jet $\eta$'
elif args.xvar == 'pt':
    draw_config['xvar'] = 'jet_pt'
    draw_config['xrange'] = [10,50]
    draw_config['xbins'] = [10,12.5,15,20,25,35,60]
    draw_config['ylimits'] = [1, 1e4]
    draw_config['xlimits'] = [0,60]
    draw_config['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
    draw_config['xunits'] = '[GeV]'
    draw_config['xlabel'] = 'Reconstructed Jet $p_T$ [GeV]'
elif args.xvar == 'bjorken_x':
    draw_config['xvar'] = 'bjorken_x'
    draw_config['xrange'] = [1e-2,1]
    draw_config['xbins'] = np.concatenate([np.arange(1e-2,1e-1,3.33333e-2),np.arange(0.1,0.5,0.2),[0.5,1.0]])
    draw_config['ylimits'] = [1e-1, 1e4]
    draw_config['xlimits'] = [1e-2,1]
    draw_config['xunits'] = ''
    draw_config['yunits'] = ''
    draw_config['xlabel'] = 'Bjorken x'
elif args.xvar == 'JB_x':
    draw_config['xvar'] = 'JB_x'
    draw_config['xrange'] = [1e-2,1]
    draw_config['xbins'] = np.concatenate([np.arange(1e-2,1e-1,3.33333e-2),np.arange(0.1,0.5,0.2),[0.5,1.0]])
    draw_config['ylimits'] = [1e-1, 1e4]
    draw_config['xlimits'] = [1e-2,1]
    draw_config['yunits'] = ''
    draw_config['xunits'] = ''
    draw_config['xlabel'] = 'Reconstructed $x_{JB}$'
    draw_config['xunits'] = ''

DrawDiffTagYieldPlot(df, draw_config)


# Finally, project the data statistical uncertainties to 100/fb of EIC data
# Overlay the PDF range variations with these predicted stat. errors

df_20rs2 = eat.UprootLoad([f"{args.dir}/CC_DIS_e10_p275_lha_20Rs2/*/out.root"], "tree", branches=branchlist)
df_21rs2 = eat.UprootLoad([f"{args.dir}/CC_DIS_e10_p275_lha_21Rs2/*/out.root"], "tree", branches=branchlist)

xvar=draw_config['xvar']
xrange=draw_config['xrange']
xbins=draw_config['xbins']
ylimits=draw_config['ylimits']
xlimits=draw_config['xlimits']
yunits=draw_config['yunits']
xunits=draw_config['xunits']

charm_ct18nnlo = eat.DifferentialTaggingYield(df, x=xvar, xrange=xrange, xbins=xbins, which='charm', process='CC_DIS_e10_p100_CT18NNLO')
charm_ct18nnlo_20rs2 = eat.DifferentialTaggingYield(df_20rs2, x=xvar, xrange=xrange, xbins=xbins, which='charm', process='CC_DIS_e10_p100_CT1820Rs2')
charm_ct18nnlo_21rs2 = eat.DifferentialTaggingYield(df_21rs2, x=xvar, xrange=xrange, xbins=xbins, which='charm', process='CC_DIS_e10_p100_CT1821Rs2')


# Calculate data statistical errors for 100/fb
#(bins, bin_widths, xs, xs_error) 

print(charm_ct18nnlo[2])
print(charm_ct18nnlo[3])



N_20 = charm_ct18nnlo_20rs2[2]
errN_20 = np.zeros(len(N_20))

for index in range(len(errN_20)):
    errN_20[index] = np.sqrt(N_20[index])

R_N_20 = np.ones(len(N_20)) + errN_20/N_20

N_21 = charm_ct18nnlo_21rs2[2]

diff_20_21 = N_21 - N_20

R_N_diff = np.ones(len(N_20)) + diff_20_21/N_20


fig, ax = plt.subplots(figsize=(12,8))
plt.axis('off')


gridspec = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])

bins = charm_ct18nnlo[0]
bin_widths = charm_ct18nnlo[1]

one_line = np.ones(len(bins))

# Log plot of differential cross-sections
ax1 = fig.add_subplot(gridspec[0, 0])

ax1.grid(which='both', axis='both')

R_N_20_label = "Stat. Uncertainty [CT18NNLO, $R_s=2s/(\overline{u}+\overline{d})=0.325$ (suppressed)]"

ax1.axhline(1.0, color='k', linestyle='-', linewidth=2)
errorboxes = [mpl.patches.Rectangle((x - xe, y - ye), 2*xe, 2*ye, facecolor='#7A6E67', alpha=0.35, 
                                    label=R_N_20_label)
              for x, y, xe, ye in zip(bins, one_line, (bin_widths/2), errN_20/N_20)]

# Create patch collection with specified colour/alpha
pc = mpl.collections.PatchCollection(errorboxes, facecolor='#7A6E67', alpha=0.35, label=R_N_20_label)

# Add collection to axes
ax1.add_collection(pc)

enhanced = ax1.errorbar(bins, R_N_diff, xerr = bin_widths/2, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='#003066', label='CT18ZNNLO with enhanced strangeness, $R_s=2s/(\overline{u}+\overline{d})=0.863$')

#ax1.fill_between(bins, pct_diff_ct18nnlo_20rs2, pct_diff_ct18nnlo_21rs2, color='#2D6CC0', alpha=0.35)

ax1.set_ylim([0.0, 2.00])
ax1.set_xlim(xlimits)
plt.ylabel('')

plt.xlabel(draw_config['xlabel'])

if xvar == 'bjorken_x' or xvar == 'JB_x':
    plt.xscale('log')
    
ax1.legend(handles=[errorboxes[0],enhanced], fontsize=18)


plt.tight_layout()

plt.savefig(f"jet_differential_tagyield_100fb_{xvar}_{args.input}.png")
plt.savefig(f"jet_differential_tagyield_100fb_{xvar}_{args.input}.pdf")
