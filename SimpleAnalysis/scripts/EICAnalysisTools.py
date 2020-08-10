import sys
import pandas as pd
import numpy as np
import uproot
import scipy
import math
import glob
from uproot_methods import *
from concurrent.futures import ThreadPoolExecutor


# Define useful variables

## Units

global u_fb,u_mb,u_pb

u_fb = 1
u_mb = 1e12*u_fb
u_pb = 1e3*u_fb

global n_gen
n_gen = 10e6



def UprootLoad(files=[], treename="", branches=[]):
    global n_gen
    print(f"::UprootLoad({files}")
    executor = ThreadPoolExecutor(8)

    #df = uproot.lazyarrays(files,
    #                       treename,
    #                       branches=branches, 
    #                       entrysteps=1000)

    filenames = glob.glob(files[0])

    data = pd.DataFrame(columns=branches)

    for filename in filenames:
        tree = uproot.open(filename)[treename]
        data = pd.concat([data, tree.arrays(branches, namedecode='utf-8', outputtype=pd.DataFrame)])

    n_gen = len(data) #uproot.numentries(files, treename, total=True)

    print(f"Originally generated number of collisions: {n_gen}")

    print("Exploding the DataFrame...")
    #data = data.set_index(['bjorken_x', 'jb_x', 'jb_Q2']).apply(pd.Series.explode).reset_index()
    #data = data.apply(pd.Series.explode)
    print("Completed exploding the DataFrame")

    return data

def TotalXSection(process='CC_DIS_e18_p275'):
    global u_mb
    if process == 'CC_DIS_e18_p275':
        return 2.408e-08*u_mb
    elif process == 'CC_DIS_e10_p100_CT18NNLO':
        return 5.44044e-09*u_mb
    elif process == 'CC_DIS_e10_p275_CT18NNLO':
        return 1.47637e-08*u_mb
    elif process == 'CC_DIS_e10_p275_CT1820Rs2':
        return 1.47637e-08*u_mb
    elif process == 'CC_DIS_e10_p100_CT1820Rs2':
        return 5.44842e-09*u_mb
    elif process == 'CC_DIS_e10_p275_CT1821Rs2':
        return 1.45884e-08*u_mb
    elif process == 'CC_DIS_e10_p100_CT1821Rs2':
        return 5.36006e-09*u_mb
    elif process == 'NC_DIS_e18_p275':
        return 3.671e-06*u_mb
    
    print(f"EICAnalysisTools::TotalXSection == Unable to find requested process, {process}.")
    # Don't continue. This is bad. Nonsense will result. The user should fix this.
    sys.exit(-1)
    return 0.0


# Compute auxiliary information for a DataFrame
def DISBjorkenX(row):
    """
    pProton      = branchParticle.At(0).P4(); #these numbers 0 , 3, 5 are hardcoded in Pythia8
    pleptonIn    = branchParticle.At(3).P4();
    pleptonOut   = branchParticle.At(5).P4();
    pPhoton      = pleptonIn - pleptonOut;
  
    #Q2, W2, Bjorken x, y, nu.
    Q2 = -pPhoton.M2()
    W2 = (pProton + pPhoton).M2()
    x = Q2 / (2. * pProton.Dot(pPhoton))
    y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn))
    """

    particles_px = row['Particle.Px']
    particles_py = row['Particle.Py']
    particles_pz = row['Particle.Pz']
    particles_E = row['Particle.E']

    # Beam Particle 4-vectors
    iproton = 0
    ilepton1 = 3
    ilepton2 = 5
    
    pProton = TLorentzVector(particles_px[iproton], particles_py[iproton], particles_pz[iproton], particles_E[iproton])
    pleptonIn = TLorentzVector(particles_px[ilepton1], particles_py[ilepton1], particles_pz[ilepton1], particles_E[ilepton1])
    pleptonOut = TLorentzVector(particles_px[ilepton2], particles_py[ilepton2], particles_pz[ilepton2], particles_E[ilepton2])
    pPhoton = pleptonIn - pleptonOut

    Q2 = -pPhoton.mag2
    W2 = (pProton + pPhoton).mag2
    x = Q2 / (2. * pProton.dot(pPhoton))

    # Write this value of x for each jet, since we need to do jet-level plots
    x_array = np.ones(len(row['jet_pt'])) * x

    return x_array




###################################################################
# Differential Computation Functions
###################################################################

# Flavour tagging study code
def DifferentialTaggingYield(df, x='jet_pt', xrange=[10,50], xbins=[10,12.5,15,20,25,30,40,50], which='all', process='CC_DIS_e18_p275', target_lumi = 100, taggers="jet_sip3dtag"):
    global u_fb, n_gen
    print(f"n_gen = {n_gen}")

    jet_flavor = np.concatenate(df['jet_flavor'].to_numpy()).ravel()
    jet_tag = np.concatenate(df['jet_sip3dtag'].to_numpy()).ravel()

    jet_x = np.concatenate(df[x].to_numpy()).ravel()

    #jet_flavor = df['jet_flavor']
    #jet_tag = df['jet_sip3dtag']

    #jet_x = df[x]

    if which == 'charm':
        jet_x = jet_x[ (jet_tag == 1) & (jet_flavor == 4) ]
    elif which == 'light':
        jet_x = jet_x[ (jet_tag == 1) & ((jet_flavor < 4) | (jet_flavor == 21)) ]


    (counts, bins) = np.histogram(jet_x, range=xrange, 
                                  bins=xbins)
                                  #bins=40)

    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths/2

    errors = np.sqrt(counts)
    rel_errors = errors/counts


    # convert counts to dsigma/dpT * 100/fb
    dsigma_dx = counts * TotalXSection(process)  * target_lumi*u_fb**(-1) / n_gen #/ bin_widths
    dsigma_dx_errors = rel_errors * dsigma_dx


    print(f"========= Validation Information from DifferentialTaggingYield ========= \n \
    \n \
    Process considered:              {process} \n \
    Theory Cross-Section:            {TotalXSection(process)} \n \
    Target Integrated Luminosity:    {target_lumi} \n \
    Total Predicted Yield of Events: {np.sum(dsigma_dx)} \ \
    \n \
    ========= Validation Information from DifferentialTaggingYield =========")
          
    
    return (bin_centers, bin_widths, dsigma_dx, dsigma_dx_errors)



def DifferentialTaggingEfficiency(df, x='jet_pt', xrange=[10,50], xbins=[10,12.5,15,20,25,30,40,50], which='all', taggers="jet_sip3dtag"):

    jet_flavor = np.concatenate(df['jet_flavor'].to_numpy()).ravel()
    jet_tag = np.concatenate(df['jet_sip3dtag'].to_numpy()).ravel()

    jet_x = np.concatenate(df[x].to_numpy()).ravel()

    # too slow
    #if len(jet_x.flatten()) != len(jet_flavor.flatten()):
    #    print("Hi!")
    #    print(len(jet_x))
    #    new_jet_x = np.array([])
    #    for i in np.arange(len(jet_x)):
    #        if len(jet_flavor[i]) > 1:
    #            print(len(jet_flavor[i]), " ", np.ones(len(jet_flavor[i])), " ", np.ones(len(jet_flavor[i]))*jet_x[i])
    #        new_jet_x = np.append(new_jet_x, np.ones(len(jet_flavor[i]))*jet_x[i])
    #print(new_jet_x[:15])

    if which == 'charm':
        jet_x = jet_x[ jet_flavor == 4 ]
        jet_tag = jet_tag[ jet_flavor == 4 ]
    elif which == 'light':
        jet_x = jet_x[ ((jet_flavor < 4) | (jet_flavor == 21)) ]
        jet_tag = jet_tag[ ((jet_flavor < 4) | (jet_flavor == 21)) ]

    (tru_counts, bins) = np.histogram(jet_x, range=xrange, 
                                      bins=xbins)

    (tag_counts, bins) = np.histogram(jet_x[ jet_tag == 1 ], range=xrange, 
                                      bins=xbins)


    eff = tag_counts/tru_counts
    n_err = scipy.stats.binom.interval(1.0-math.exp(-1), tru_counts, eff) 

    eff_err=[np.fabs(n_err[0]/tru_counts-eff), np.fabs(n_err[1]/tru_counts-eff)]


    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths/2

    print(f"========= Validation Information from DifferentialTaggingEfficiency =========")
    print(eff)
    print(eff_err)
    print(f"Tagging Efficiency above X threshold for {which}:")
    for ibin in np.arange(len(bins)-1):
        #mean_eff = np.nanmean(eff[ibin:])
        #avg_errors = (eff_err[0]+eff_err[1])/2.0
        #mean_err = np.sqrt(np.sum((avg_errors[ibin:]/eff[ibin:])**2))*mean_eff
        n_tag = np.sum(tag_counts[ibin:])
        n_tru = np.sum(tru_counts[ibin:])
        eff_above_x = n_tag/n_tru
        err_above_x = np.nanmean([np.fabs(np.sum(n_err[0][ibin:])/n_tru-eff_above_x), np.fabs(np.sum(n_err[1][ibin:])/n_tru-eff_above_x)])
        #print(f"Minimum X: {bins[ibin]:.1f}   Tag. Efficiency: ({mean_eff*100:.3f} +/- {mean_err*100:.3f})%")
        print(f"Minimum X: {bins[ibin]:.1f}   Tag. Efficiency: ({eff_above_x*100:.3f} +/- {err_above_x*100:.3f})%")
    

    print(f"========= Validation Information from DifferentialTaggingEfficiency =========")

    return (bin_centers, bin_widths, eff, eff_err)


