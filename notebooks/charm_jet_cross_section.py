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

import EICAnalysisTools as eat


# Units
global u_fb, u_mb, xsection
u_fb = 1
u_mb = 1e12*u_fb
#xsection = 2.408e-08*u_mb

## dataset dictionary (cross-sections, etc.)
#global xs_map
#xs_map = {}
#xs_map['CC_DIS_e18_p275_B15'] = 2.408e-08*u_mb # from PYTHIA8
#xs_map['NC_DIS_e18_p275_B15'] = 3.671e-06*u_mb # from PYTHIA8


print("Loading data...")

branchlist = ["Jet.PT", "Jet.Eta", "Jet.Flavor", "Jet.BTag", "GenJet.PT", "GenJet.Eta", "GenJet.Flavor"]


df_ct18nnlo = eat.UprootLoad([f"../CC_DIS_e10_p275_B15_dR5_maxIP3mm_trkpt10_22sigmin_CT18NNLO/*/out.root"], "Delphes",
                             branches=branchlist)
#df_ct18nnlo_20rs2 = eat.UprootLoad([f"../CC_DIS_e10_p275_B15_dR5_maxIP3mm_trkpt10_22sigmin_lha_20Rs2/*/out.root"], "Delphes",
#                             branches=branchlist)
#df_ct18nnlo_21rs2 = eat.UprootLoad([f"../CC_DIS_e10_p275_B15_dR5_maxIP3mm_trkpt10_22sigmin_lha_21Rs2/*/out.root"], "Delphes",
#                             branches=branchlist)


#df_ct18nnlo = df_ct18nnlo[:1000]
#df_ct18nnlo_20rs2 = df_ct18nnlo_20rs2[:1000]
#df_ct18nnlo_21rs2 = df_ct18nnlo_21rs2[:1000]


#df_ct18nlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18NLO/*/out.root")
#df_ct18Annlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18ANNLO/*/out.root")
#df_ct18Anlo = UprootLoad(f"../CC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18ANLO/*/out.root")

#df_ct18nnlo_nc = UprootLoad(f"../NC_DIS_e18_p275_B15_dR5_maxIP3mm_CT18NNLO/*/out.root")


print("Done loading data for the study!")



def DifferentialXS(df, x='GenJet.PT', xrange=[10,50], xbins=[10,12.5,15,20,25,30,40,50], which='all', process='CC_DIS_e18_p275'):
    global u_fb, u_mb, n_gen
    n_gen = len(df)
    print(f"n_gen = {n_gen}")
    jet_pt = np.concatenate(df['GenJet.PT'].to_numpy()).ravel()
    jet_eta = np.concatenate(df['GenJet.Eta'].to_numpy()).ravel()
    jet_flavor = np.concatenate(df['GenJet.Flavor'].to_numpy()).ravel()

    jet_x = np.concatenate(df[x].to_numpy()).ravel()


    #jet_basics = (0.01 < jet_eta) & (jet_eta < 0.9)
    jet_basics = (-5.0 < jet_eta) & (jet_eta < 5.0)
    
    all_flavor = ( jet_flavor > -999.0 ) & (jet_basics)
    light_flavor = ( jet_flavor < 4 ) & (jet_basics)
    charm_flavor = ( jet_flavor == 4 ) & (jet_basics)

    selection = all_flavor
    if which == 'charm':
        selection = charm_flavor


    (counts, bins) = np.histogram(jet_x[ selection ], range=xrange, 
                                  bins=xbins)
                                  #bins=40)

    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths/2

    errors = np.sqrt(counts)
    rel_errors = errors/counts


    # convert counts to dsigma/dpT * 100/fb
    dsigma_dx = counts * eat.TotalXSection(process)  * 100*u_fb**(-1) / n_gen / bin_widths
    dsigma_dx_errors = rel_errors * dsigma_dx
    
    return (bin_centers, bin_widths, dsigma_dx, dsigma_dx_errors)





draw_config_eta={}

draw_config_eta['xvar'] = 'GenJet.Eta'
draw_config_eta['xrange'] = [-4.0, 4.0]
draw_config_eta['xbins'] = [-4.0, -3.5, -2.5, -1.0, 1.0, 2.5, 3.5, 4.0]
draw_config_eta['ylimits'] = [10, 1e7]
draw_config_eta['xlimits'] = [-5,5]
draw_config_eta['yunits'] = ''
draw_config_eta['xunits'] = ''

draw_config_pt={}
draw_config_pt['xvar'] = 'GenJet.PT'
draw_config_pt['xrange'] = [10,50]
draw_config_pt['xbins'] = [10,12.5,15,20,25,30,40,50]
draw_config_pt['ylimits'] = [10, 2e5]
draw_config_pt['xlimits'] = [0,60]
draw_config_pt['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
draw_config_pt['xunits'] = '[GeV]'


def DrawDiffXSPlot(draw_config={}, process='CC_DIS_e18_p275'):
    xvar=draw_config['xvar']
    xrange=draw_config['xrange']
    xbins=draw_config['xbins']
    ylimits=draw_config['ylimits']
    xlimits=draw_config['xlimits']
    yunits=draw_config['yunits']
    xunits=draw_config['xunits']
    
    fig, ax = plt.subplots(figsize=(12,12))
    plt.axis('off')
    
    
    
    if process.find('CC') != -1:
        gridspec = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
        
        
        # Log plot of differential cross-sections
        ax1 = fig.add_subplot(gridspec[0, 0])
        
        ax1.grid(which='both', axis='both')
        (bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nnlo, x=xvar, xrange=xrange, xbins=xbins, which='all', process='CC_DIS_e10_p275_CT18NNLO')
        ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='red')

        #(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nnlo_20rs2, x=xvar, xrange=xrange, xbins=xbins, which='all', process='CC_DIS_e10_p275_CT1820Rs2')
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NNLO 20Rs2)', marker='o', ms=10, ls='none', linewidth=2, color='saddlebrown')

        #(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nnlo_21rs2, x=xvar, xrange=xrange, xbins=xbins, which='all', process='CC_DIS_e10_p275_CT1821Rs2')
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NNLO 21Rs2)', marker='o', ms=10, ls='none', linewidth=2, color='magenta')
        
        #(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18Annlo, x=xvar, xrange=xrange, xbins=xbins, which='all')
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18ANNLO)', marker='s', ms=10, ls='none', linewidth=2, color='saddlebrown')
        
        #(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nlo, x=xvar, xrange=xrange, xbins=xbins, which='all')
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NLO)', marker='v', ms=10, ls='none', linewidth=2, color='magenta')
        
        #(bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18Anlo, x=xvar, xrange=xrange, xbins=xbins, which='all')
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18ANLO)', marker='D', ms=10, ls='none', linewidth=2, color='coral')
        
        charm_ct18nnlo = DifferentialXS(df_ct18nnlo, x=xvar, xrange=xrange, xbins=xbins, which='charm')
        (bins, bin_widths, xs, xs_error) = charm_ct18nnlo
        ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='blue')

        #charm_ct18nnlo_20rs2 = DifferentialXS(df_ct18nnlo_20rs2, x=xvar, xrange=xrange, xbins=xbins, which='charm')
        #(bins, bin_widths, xs, xs_error) = charm_ct18nnlo_20rs2
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO 20Rs2)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')

        #charm_ct18nnlo_21rs2 = DifferentialXS(df_ct18nnlo_21rs2, x=xvar, xrange=xrange, xbins=xbins, which='charm')
        #(bins, bin_widths, xs, xs_error) = charm_ct18nnlo_21rs2
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO 21rs2)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='magenta')
        
        #charm_ct18Annlo = DifferentialXS(df_ct18Annlo, x=xvar, xrange=xrange, xbins=xbins, which='charm')
        #(bins, bin_widths, xs, xs_error) = charm_ct18Annlo
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18ANNLO)', marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')
        
        #charm_ct18nlo = DifferentialXS(df_ct18nlo, x=xvar, xrange=xrange, xbins=xbins, which='charm')
        #(bins, bin_widths, xs, xs_error) = charm_ct18nlo
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NLO)', marker='v', ms=10, ls='none', linewidth=2, fillstyle='none', color='magenta')
        
        #charm_ct18Anlo = DifferentialXS(df_ct18Anlo, x=xvar, xrange=xrange, xbins=xbins, which='charm')
        #(bins, bin_widths, xs, xs_error) = charm_ct18Anlo
        #ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18ANLO)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='coral')
        
        
        xvar_symbol="p_T"
        if xvar.find('Eta') != -1:
            xvar_symbol='\\eta'
            
        plt.ylabel(f'$d\\sigma/d{xvar_symbol} \\times 100\\mathrm{{fb^{{-1}}}}$ {yunits}')
        plt.xlabel(f'Generated Jet ${xvar_symbol}$ [GeV]')
        plt.yscale('log')
        
        plt.title("CC-DIS, 10x275GeV, $Q^2>100\\mathrm{GeV^2}$", fontsize=20)
        
        ax1.set_ylim(ylimits)
        ax1.set_xlim(xlimits)
    
        ax1.legend(fontsize=18)
    
        
        # ratio plot for charm production
        #ax2 = fig.add_subplot(gridspec[1, 0])
        
        #ax2.grid(which='both', axis='both')
        
        #ratio_charm_ct18nnlo_ct18nnlo = charm_ct18nnlo[2]/charm_ct18nnlo[2]
        #ratio_charm_ct18nnlo20rs2_ct18nnlo = charm_ct18nnlo_20rs2[2]/charm_ct18nnlo[2]
        #ratio_charm_ct18nnlo21rs2_ct18nnlo = charm_ct18nnlo_21rs2[2]/charm_ct18nnlo[2]
        #ratio_charm_ct18Annlo_ct18nnlo = charm_ct18Annlo[2]/charm_ct18nnlo[2]
        #ratio_charm_ct18nlo_ct18nnlo = charm_ct18nlo[2]/charm_ct18nnlo[2]
        #ratio_charm_ct18Anlo_ct18nnlo = charm_ct18Anlo[2]/charm_ct18nnlo[2]
        
        #ratioerr_charm_ct18nnlo_ct18nnlo = ratio_charm_ct18nnlo_ct18nnlo*np.sqrt((charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
        #ratioerr_charm_ct18nnlo20rs2_ct18nnlo = ratio_charm_ct18nnlo20rs2_ct18nnlo*np.sqrt((charm_ct18nnlo_20rs2[3]/charm_ct18nnlo_20rs2[2])**2)
        #ratioerr_charm_ct18nnlo21rs2_ct18nnlo = ratio_charm_ct18nnlo21rs2_ct18nnlo*np.sqrt((charm_ct18nnlo_21rs2[3]/charm_ct18nnlo_21rs2[2])**2)
        #ratioerr_charm_ct18Annlo_ct18nnlo = ratio_charm_ct18Annlo_ct18nnlo*np.sqrt((charm_ct18Annlo[3]/charm_ct18Annlo[2])**2 + (charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
        #ratioerr_charm_ct18nlo_ct18nnlo = ratio_charm_ct18nlo_ct18nnlo*np.sqrt((charm_ct18nlo[3]/charm_ct18nlo[2])**2 + (charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
        #ratioerr_charm_ct18Anlo_ct18nnlo = ratio_charm_ct18Anlo_ct18nnlo*np.sqrt((charm_ct18Anlo[3]/charm_ct18Anlo[2])**2 + (charm_ct18nnlo[3]/charm_ct18nnlo[2])**2)
        
        #print(ratio_charm_ct18nnlo20rs2_ct18nnlo)
        #print(ratio_charm_ct18nnlo21rs2_ct18nnlo)
        #print(ratio_charm_ct18Annlo_ct18nnlo)
        #print(ratio_charm_ct18nlo_ct18nnlo)
        #print(ratio_charm_ct18Anlo_ct18nnlo)
        
        #ax2.axhline(1.0, color='k', linestyle='-', linewidth=2)
        #ax2.errorbar(bins, ratio_charm_ct18nnlo_ct18nnlo, xerr = bin_widths/2, yerr=ratioerr_charm_ct18nnlo_ct18nnlo, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')
        #errorboxes = [mpl.patches.Rectangle((x - xe, y - ye), 2*xe, 2*ye)
        #              for x, y, xe, ye in zip(bins, ratio_charm_ct18nnlo_ct18nnlo, (bin_widths/2), ratioerr_charm_ct18nnlo_ct18nnlo)]
        
        # Create patch collection with specified colour/alpha
        #pc = mpl.collections.PatchCollection(errorboxes, facecolor='k', alpha=0.25)
        
        # Add collection to axes
        #ax2.add_collection(pc)
        
        #ax2.errorbar(bins, ratio_charm_ct18nnlo20rs2_ct18nnlo, xerr = bin_widths/2, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')
        #ax2.errorbar(bins, ratio_charm_ct18nnlo21rs2_ct18nnlo, xerr = bin_widths/2, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='magenta')
        #ax2.errorbar(bins, ratio_charm_ct18Annlo_ct18nnlo, xerr = bin_widths/2, marker='s', ms=10, ls='none', linewidth=2, fillstyle='none', color='saddlebrown')
        #ax2.errorbar(bins, ratio_charm_ct18nlo_ct18nnlo, xerr = bin_widths/2, marker='v', ms=10, ls='none', linewidth=2, fillstyle='none', color='magenta')
        #ax2.errorbar(bins, ratio_charm_ct18Anlo_ct18nnlo, xerr = bin_widths/2, marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='coral')
        
        #ax2.set_ylim([0.25, 1.75])
        #ax2.set_xlim(xlimits)
        #plt.ylabel('Ratio to\nCT18NNLO', fontsize=18)
    elif process.find('NC') != -1:
        gridspec = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
        
        
        # Log plot of differential cross-sections
        ax1 = fig.add_subplot(gridspec[0, 0])
        
        ax1.grid(which='both', axis='both')
        (bins, bin_widths, xs, xs_error) = DifferentialXS(df_ct18nnlo_nc, x=xvar, xrange=xrange, xbins=xbins, which='all', process=process)
        ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='all jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='red')

        charm_ct18nnlo = DifferentialXS(df_ct18nnlo_nc, x=xvar, xrange=xrange, xbins=xbins, which='charm', process=process)
        (bins, bin_widths, xs, xs_error) = charm_ct18nnlo
        ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='blue')

        xvar_symbol="p_T"
        if xvar.find('Eta') != -1:
            xvar_symbol='\\eta'
            
        plt.ylabel(f'$d\\sigma/d{xvar_symbol} \\times 100\\mathrm{{fb^{{-1}}}}$ {yunits}')
        plt.xlabel(f'Generated Jet ${xvar_symbol}$ [GeV]')
        plt.yscale('log')
        
        plt.title("NC-DIS, 18x275GeV, $Q^2>100\\mathrm{GeV^2}$", fontsize=20)
        
        ax1.set_ylim(ylimits)
        ax1.set_xlim(xlimits)
    
        ax1.legend(fontsize=18)

        
    plt.tight_layout()
    
    plt.savefig(f"jet_differential_xs_{xvar}_{process}.png")
    plt.savefig(f"jet_differential_xs_{xvar}_{process}.pdf")

### Charged Current Scattering
DrawDiffXSPlot(draw_config_pt, process='CC_DIS_e10_p275')
DrawDiffXSPlot(draw_config_eta, process='CC_DIS_e10_p275')

sys.exit()

### Neutral Current Scattering
draw_config_pt={}
draw_config_pt['xvar'] = 'GenJet.PT'
draw_config_pt['xrange'] = [10,50]
draw_config_pt['xbins'] = [10,12.5,15,20,25,30,40,50]
draw_config_pt['ylimits'] = [1e2, 5e8]
draw_config_pt['xlimits'] = [0,60]
draw_config_pt['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
draw_config_pt['xunits'] = '[GeV]'
#DrawDiffXSPlot(draw_config_pt, process='NC_DIS_e18_p275')


# Flavour tagging study code
def DifferentialTaggingYield(df, x='Jet.PT', xrange=[10,50], xbins=[10,12.5,15,20,25,30,40,50], which='all', process='CC_DIS_e18_p275_B15'):
    global u_fb, u_mb, xs_map, n_gen
    n_gen = len(df)
    print(f"n_gen = {n_gen}")
    jet_pt = np.concatenate(df['Jet.PT'].to_numpy()).ravel()
    jet_eta = np.concatenate(df['Jet.Eta'].to_numpy()).ravel()
    jet_flavor = np.concatenate(df['Jet.Flavor'].to_numpy()).ravel()
    jet_tag = np.concatenate(df['Jet.BTag'].to_numpy()).ravel()

    jet_x = np.concatenate(df[x].to_numpy()).ravel()


    #jet_basics = (0.01 < jet_eta) & (jet_eta < 0.9)
    jet_basics = (jet_tag == 1)
    
    all_flavor = ( jet_flavor > -999.0 ) & (jet_basics)
    light_flavor = ( jet_flavor < 4 ) & (jet_basics)
    charm_flavor = ( jet_flavor == 4 ) & (jet_basics)

    selection = all_flavor
    if which == 'charm':
        selection = charm_flavor
    elif which == 'light':
        selection = light_flavor

    (counts, bins) = np.histogram(jet_x[ selection ], range=xrange, 
                                  bins=xbins)
                                  #bins=40)

    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths/2

    errors = np.sqrt(counts)
    rel_errors = errors/counts


    # convert counts to dsigma/dpT * 100/fb
    dsigma_dx = counts * eat.TotalXSection(process)  * 100*u_fb**(-1) / n_gen / bin_widths
    dsigma_dx_errors = rel_errors * dsigma_dx
    
    return (bin_centers, bin_widths, dsigma_dx, dsigma_dx_errors)



# Inputs

draw_config_eta={}

draw_config_eta['xvar'] = 'Jet.Eta'
draw_config_eta['xrange'] = [-4.0, 4.0]
draw_config_eta['xbins'] = [-4.0, -3.5, -2.5, -1.0, 1.0, 2.5, 3.5, 4.0]
draw_config_eta['ylimits'] = [1, 1e3]
draw_config_eta['xlimits'] = [-5,5]
draw_config_eta['yunits'] = ''
draw_config_eta['xunits'] = ''
draw_config_eta['data'] = df_ct18nnlo

draw_config_pt={}
draw_config_pt['xvar'] = 'Jet.PT'
draw_config_pt['xrange'] = [10,50]
draw_config_pt['xbins'] = [10,12.5,15,20,25,30,40,50]
draw_config_pt['ylimits'] = [1, 1e3]
draw_config_pt['xlimits'] = [0,60]
draw_config_pt['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
draw_config_pt['xunits'] = '[GeV]'
draw_config_pt['data'] = df_ct18nnlo




def DrawDiffTagYieldPlot(draw_config={}, process='CC_DIS_e18_p275'):
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
    
    (bins, bin_widths, xs, xs_error) = DifferentialTaggingYield(draw_config['data'], x=xvar, xrange=xrange, xbins=xbins, which='light', process=process)
    ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='light jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='red')
    
    
    charm_ct18nnlo = DifferentialTaggingYield(draw_config['data'], x=xvar, xrange=xrange, xbins=xbins, which='charm', process=process)
    (bins, bin_widths, xs, xs_error) = charm_ct18nnlo
    ax1.errorbar(bins, xs, xerr = bin_widths/2, yerr=xs_error, label='charm jets (CT18NNLO)', marker='D', ms=10, ls='none', linewidth=2, fillstyle='none', color='blue')
    
    
    
    xvar_symbol="p_T"
    if xvar.find('Eta') != -1:
        xvar_symbol='\\eta'
        
    plt.ylabel(f'$\\varepsilon_{{tag}} \\times d\\sigma/d{xvar_symbol} \\times 100\\mathrm{{fb^{{-1}}}}$ {yunits}')
    plt.xlabel(f'Reconstructed Jet ${xvar_symbol}$ {xunits}')
    plt.yscale('log')
    
    process_name = 'CC-DIS'
    if process.find('NC') != -1:
        process_name = 'NC-DIS'
    plt.title(f"{process_name}, 18x275GeV, $Q^2>100\\mathrm{{GeV^2}}$", fontsize=20)
    
    ax1.set_ylim(ylimits)
    ax1.set_xlim(xlimits)
    
    ax1.legend(fontsize=18)
    
    
    plt.tight_layout()
    
    plt.savefig(f"jet_differential_tagyield_{xvar}_{process}.png")
    plt.savefig(f"jet_differential_tagyield_{xvar}_{process}.pdf")


#DrawDiffTagYieldPlot(draw_config_pt)
#DrawDiffTagYieldPlot(draw_config_eta)

# Neutral Current Process
draw_config_pt={}
draw_config_pt['xvar'] = 'Jet.PT'
draw_config_pt['xrange'] = [10,50]
draw_config_pt['xbins'] = [10,12.5,15,20,25,30,40,50]
draw_config_pt['ylimits'] = [1e-1, 1e6]
draw_config_pt['xlimits'] = [0,60]
draw_config_pt['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
draw_config_pt['xunits'] = '[GeV]'
draw_config_pt['data'] = df_ct18nnlo_nc
#DrawDiffTagYieldPlot(draw_config_pt, process='NC_DIS_e18_p275_B15')


### Tagging Efficiency
# Flavour tagging study code
def DifferentialTaggingEfficiency(df, x='Jet.PT', xrange=[10,50], xbins=[10,12.5,15,20,25,30,40,50], which='all'):
    n_gen = len(df)
    print(f"n_gen = {n_gen}")
    jet_pt = np.concatenate(df['Jet.PT'].to_numpy()).ravel()
    jet_eta = np.concatenate(df['Jet.Eta'].to_numpy()).ravel()
    jet_flavor = np.concatenate(df['Jet.Flavor'].to_numpy()).ravel()
    jet_tag = np.concatenate(df['Jet.BTag'].to_numpy()).ravel()

    jet_x = np.concatenate(df[x].to_numpy()).ravel()


    jet_tagged = (jet_tag == 1)
    
    all_flavor = ( jet_flavor > -999.0 ) 
    light_flavor = ( jet_flavor < 4 ) 
    charm_flavor = ( jet_flavor == 4 ) 

    selection = all_flavor
    if which == 'charm':
        selection = charm_flavor
    elif which == 'light':
        selection = light_flavor

    (tru_counts, bins) = np.histogram(jet_x[ selection ], range=xrange, 
                                      bins=xbins)

    (tag_counts, bins) = np.histogram(jet_x[ (selection) & (jet_tagged) ], range=xrange, 
                                      bins=xbins)


    eff = tag_counts/tru_counts
    n_err = scipy.stats.binom.interval(1.0-math.exp(-1), tru_counts, eff) 

    eff_err=[np.fabs(n_err[0]/tru_counts-eff), np.fabs(n_err[1]/tru_counts-eff)]

    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths/2

    
    return (bin_centers, bin_widths, eff, eff_err)


# Inputs
draw_config_eta={}

draw_config_eta['xvar'] = 'Jet.Eta'
draw_config_eta['xrange'] = [-4.0, 4.0]
draw_config_eta['xbins'] = [-4.0, -3.5, -2.5, -1.0, 1.0, 2.5, 3.5, 4.0]
draw_config_eta['ylimits'] = [1e-5, 1]
draw_config_eta['xlimits'] = [-5,5]
draw_config_eta['yunits'] = ''
draw_config_eta['xunits'] = ''

draw_config_pt={}
draw_config_pt['xvar'] = 'Jet.PT'
draw_config_pt['xrange'] = [10,50]
draw_config_pt['xbins'] = [10,12.5,15,20,25,30,40,50]
draw_config_pt['ylimits'] = [1e-5, 1]
draw_config_pt['xlimits'] = [0,60]
draw_config_pt['yunits'] = '[$\\mathrm{GeV^{-1}}$]'
draw_config_pt['xunits'] = '[GeV]'



def DrawDiffTagEfficiencyPlot(draw_config={}):
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
    
    (bins, bin_widths, eff, eff_error) = DifferentialTaggingEfficiency(df_ct18nnlo, x=xvar, xrange=xrange, xbins=xbins, which='light')
    ax1.errorbar(bins, eff, xerr = bin_widths/2, yerr=eff_error, label='light jets (CT18NNLO)', marker='o', ms=10, ls='none', linewidth=2, color='red')
    
    
    (bins, bin_widths, eff, eff_error) = DifferentialTaggingEfficiency(df_ct18nnlo, x=xvar, xrange=xrange, xbins=xbins, which='charm')
    ax1.errorbar(bins, eff, xerr = bin_widths/2, yerr=eff_error, label='charm jets (CT18NNLO)', marker='D', ms=10, ls='none', linewidth=2, color='blue')
    
    
    
    xvar_symbol="p_T"
    if xvar.find('Eta') != -1:
        xvar_symbol='\\eta'
        
    plt.ylabel(f'$\\varepsilon_{{tag}}$ {yunits}')
    plt.xlabel(f'Reconstructed Jet ${xvar_symbol}$ {xunits}')
    
    plt.title("CC-DIS, 18x275GeV, $Q^2>100\\mathrm{GeV^2}$", fontsize=20)
    
    ax1.set_ylim(ylimits)
    ax1.set_xlim(xlimits)
    
    ax1.legend(fontsize=18)
    
    plt.yscale('log')
    
    plt.tight_layout()
    
    plt.savefig(f"jet_differential_tagefficiency_{xvar}.png")
    plt.savefig(f"jet_differential_tagefficiency_{xvar}.pdf")

#DrawDiffTagEfficiencyPlot(draw_config_pt)
#DrawDiffTagEfficiencyPlot(draw_config_eta)
