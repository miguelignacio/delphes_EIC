######################################################################################################################                                
# EIC detector model                                                                                                                                                                                     
# Author: Miguel Arratia                                 
# email:miguel.arratia@ucr.edu
# based on parameters from EIC detector handbook v1.2 http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf
# as well as on assumptions on calorimeter granularity and tracking efficiency (not specified in handbook). 
#######################################################################################################################


#######################################
# Load any external configurations
#######################################

source customizations.tcl


#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency
 
  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger
  TrackSmearing

  ECal
  HCal
 
  Calorimeter
  EFlowMerger
  EFlowFilter
  
  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter
  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET
  
  FastJetFinder

  JetEnergyScale

  JetFlavorAssociation
  GenJetFlavorAssociation

  UniqueObjectFinder

  ScalarHT

  TrackCountingBTagging

  mRICHPID

  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
    set InputArray Delphes/stableParticles
    
    set OutputArray stableParticles
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons
    
    
    #Values taken from EIC detector handbook v1.2
    # radius of the magnetic field coverage, in m
    set Radius 0.8
    # half-length of the magnetic field coverage, in m
    set HalfLength 1.00
    
    # magnetic field
    set Bz $PARAM_BZ
}


####################################
# Common Tracking Efficiency Model
####################################

set CommonTrackingEfficiency {
    (pt <= 0.1)                                                       * (0.00) +
    (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)                     * (1.0) +
    (abs(eta) <= 1.5) * (pt > 1.0)                                    * (1.0) +
    (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (1.0) +
    (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (1.0) +
    (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (1.0) +
    (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0)                  * (1.0) +
    (abs(eta) > 3.0)                                                  * (0.00)
}


####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set EfficiencyFormula $CommonTrackingEfficiency
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  set EfficiencyFormula $CommonTrackingEfficiency

}

##############################
# Muon tracking efficiency
##############################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  set EfficiencyFormula $CommonTrackingEfficiency

}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons. 
  # Based on EIC detector handbook v1.2 

  set ResolutionFormula {                  (abs(eta) <= 1.0) * (pt > 0.1) * sqrt((5e-3)^2 + (pt*cosh(eta))^2*(5e-4)^2) +
                         (abs(eta) > 1.0 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt((1e-2)^2 + (pt*cosh(eta))^2*(5e-4)^2) +
	                 (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1) * sqrt((2e-2)^2 + (pt*cosh(eta))^2*(1e-3)^2)  }
}

###################################
# Momentum resolution for muons
###################################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}
  # Resolution is parametrized as a function of p, not pT, in the matrix.
    
  # resolution formula for charged hadrons. 
  # Based on EIC detector handbook v1.2 
  set ResolutionFormula {                  (abs(eta) <= 1.0) * (pt > 0.1) * sqrt((5e-3)^2 + (pt*cosh(eta))^2*(5e-4)^2) +
                         (abs(eta) > 1.0 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt((1e-2)^2 + (pt*cosh(eta))^2*(5e-4)^2) +
	                 (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1) * sqrt((2e-2)^2 + (pt*cosh(eta))^2*(1e-3)^2)  }
}



###################################
# Momentum resolution for electrons
###################################

##Electron needs some thinking for EIC


module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}
    # resolution formula for electrons; same as above. Needs some thinking on how to combine with EMCAL
  
  set ResolutionFormula {                  (abs(eta) <= 1.0) * (pt > 0.1) * sqrt((5e-3)^2 + (pt*cosh(eta))^2*(5e-4)^2) +
                         (abs(eta) > 1.0 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt((1e-2)^2 + (pt*cosh(eta))^2*(5e-4)^2) +
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1) * sqrt((2e-2)^2 + (pt*cosh(eta))^2*(1e-3)^2)  }
}


##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}

################################                                                                    
# Track impact parameter smearing                                                                   
################################                                                                    


module TrackSmearing TrackSmearing {
  set InputArray TrackMerger/tracks
#  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray tracks
#  set ApplyToPileUp true

  # magnetic field
  set Bz $PARAM_BZ

  set PResolutionFormula { 0.0 }
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }
  set D0ResolutionFormula "( abs(eta) <= 3.0 ) * ( pt > 0.1 ) * $PARAM_D0RES"
  
  set DZResolutionFormula "( abs(eta) <= 3.0 ) * ( pt > 0.1 ) * $PARAM_DZRES"
  
    
}


#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackSmearing/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true
  set EnergyMin 0.1  
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # Granularity is not discussed in EIC detector handbook. 


  ##BARREL
 # assume 0.0174 x 0.020 resolution in phi,eta in the barrel |eta| < 1.0
  set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi/180.0}]
  }

 #deta=0.02 units for |eta| <=1.0
  for {set i -50} {$i < 50} {incr i} {
        set eta [expr {$i * 0.02}]
        add EtaPhiBins $eta $PhiBins
  }

  ##ENDCAP
  # assume 0.020 x 0.020 resolution in phi,eta in the endcap 1.0<|eta| < 4.0
  set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi/180.0}]
  }

  #deta=0.02 units for 1.0 < |eta| <= 4.0
  #first, from -4.0 to -1.0
  for {set i 1} {$i <=151} {incr i} {
        set eta [expr {-4.02 + $i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }
  #same for 1.0 to 4.0
    for  {set i 1} {$i <=151} {incr i} {
        set eta [expr {0.98 + $i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }

 

  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  # set ResolutionFormula {resolution formula as a function of eta and energy}
  # EM - W/ScFi, granularity 2.5 cm x 2.5 cm, 12% stochastic, 2% constant terms, as suggested by Oleg Tsai.
   

   set ResolutionFormula {    (eta <= -2.0 && eta>-4.0)                          * sqrt(energy^2*0.01^2 + energy*0.02^2)+ \
                              (eta <= -1.0 && eta>-2.0 )                         * sqrt(energy^2*0.01^2 + energy*0.07^2)+ \ 
                              (eta <= 1.0  && eta> -1.0 )                        * sqrt(energy^2*0.01^2 + energy*0.10^2)+ \
                              (eta <= 4.0  &&  eta>1.0 )                         * sqrt(energy^2*0.02^2 + energy*0.12^2)} 


}


#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false

  ##Assumes noise 100 MeV per tower. 
  set EnergyMin 0.4
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # Granularity is not discussed in EIC detector handbook. 
 ## Barrel = 0.1x0.1 , following sPHENIX HCAL

    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    
    ## Coverage is -4.0, -1.0 , and +1.0 to 4.0. 
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
        
    for {set i 1} {$i <=31} {incr i} {
	set eta [expr {-4.1 + $i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }
    
    for {set i 1} {$i <=31} {incr i} {
	set eta [expr {0.9 + $i * 0.1 }]
	add EtaPhiBins $eta $PhiBins
    }


  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  ##Resolution in endcaps based on EIC detector handbook 1.2
  ## Resolution midrapidity, as per sPHENIX HCAL

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {    (eta <= -1.0 && eta>-4.0)                       * sqrt(energy^2*0.10^2 + energy*0.50^2)+
                             (eta <= 1.0 && eta>-1.0 )                       * sqrt(energy^2*0.10^2 + energy*1.00^2)+
                              (eta <= 3.0  && eta>1.0 )                       * sqrt(energy^2*0.10^2 + energy*0.50^2)+
                              (eta <= 4.0  && eta>3.0 )                       * sqrt(energy^2*0.05^2 + energy*0.45^2)}

}


#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray chargedHadrons
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}



####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray HCal/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

######################
# EFlowFilter
######################

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray eflow
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) +
                         (abs(eta) > 2.5)                                   * (0.00)}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
  set EfficiencyFormula {                                      (pt <= 1.0) * (0.00) +
                                           (abs(eta) <= 3.5) * (pt > 1.0)  * (1.0) }
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons

  set EnergyOutputArray energy
}


#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}


#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set JetPTMin 1.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}



############
# Jet finder
############

module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeTrimming 1
  set RTrim 0.4
  set PtFracTrim 0.20
  #set PtFracTrim 0.05

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 1.0}





##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets (do not apply it)
  set ScaleFormula {1.0} 
}

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 4.0
  set PartonEtaMax 4.0

}

module JetFlavorAssociation GenJetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

}



#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray JetEnergyScale/jets jets
}

############################
# b-tagging (track counting)
############################

module TrackCountingBTagging TrackCountingBTagging {
    set JetInputArray JetEnergyScale/jets
    set TrackInputArray HCal/eflowTracks
    
    set BitNumber 0
    
    # maximum distance between jet and track
    set DeltaR 0.5
    
    # minimum pt of tracks
    set TrackPtMin $PARAM_TRACKPTMIN
    #set TrackPtMin 1.0
    
    # maximum transverse impact parameter (in mm)
    set TrackIPMax $PARAM_TRACKIPMAX
    #set TrackIPMax 3
    
    # minimum ip significance for the track to be counted
    set SigMin $PARAM_TRACKSIGMIN
    #set SigMin 2.0
    set Use3D true
    # alternate setting for 2D IP (default)
    #  set SigMin 1.3
    #  set Use3D false
    
    # minimum number of tracks (high efficiency n=2, high purity n=3)
    #set Ntracks 3
    set Ntracks $PARAM_TRACKSNMIN
}

##################
# PID Systems
##################

# mRICH

module IdentificationMap mRICHPID {
  set InputArray HCal/eflowTracks
  set OutputArray tracks

  # {PID in} {PID out} {formula}
  # make sure "PID in" and "PID out" have the same charge (e.g {-13} {211} or {-321} {211})
  # {211} {-13} is equivalent to {-211} {13} (and needs to be written once only...)



 # From EIC Yellow Report PID studies, the mRICH coverage is -3.5<eta<-1 and 1<eta<2


  # --- pions ---

  add EfficiencyFormula {-11} {-11} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.99) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.98) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.96) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.92) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.87) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.82) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.77) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.73) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.69) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.65) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.62) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.60) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.58) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.57) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.56) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.55) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.54) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.53) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.53) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.50) } 

    add EfficiencyFormula {-11} {211} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.08) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.13) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.27) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.31) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.35) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.38) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.43) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.44) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.45) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.46) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.50) } 

    add EfficiencyFormula {211} {211} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.98) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.96) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.92) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.87) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.82) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.77) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.73) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.69) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.65) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.62) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.60) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.58) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.57) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.56) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.55) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.54) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.53) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.53) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.51) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.46) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.46) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.45) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.45) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.44) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.44) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.43) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.43) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.41) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.41) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.39) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.39) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.39) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.38) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.38) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.38) } 

    add EfficiencyFormula {211} {-11} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.08) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.13) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.27) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.31) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.35) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.38) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.43) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.44) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.45) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.46) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.50) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.49) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.48) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.47) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.46) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.46) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.45) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.45) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.44) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.44) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.43) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.43) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.42) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.41) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.41) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.40) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.39) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.39) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.39) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.38) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.38) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.38) } 

    add EfficiencyFormula {211} {321} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.03) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.05) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.06) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.07) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.08) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.09) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.10) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.11) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.12) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.13) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.14) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.15) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.16) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.17) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.19) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.20) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.21) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.21) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.22) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.24) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.24) } 

    add EfficiencyFormula {321} {321} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (1.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.99) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.99) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.99) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.99) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.98) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.98) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.97) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.96) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.95) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.94) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.92) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.91) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.89) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.88) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.86) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.84) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.82) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.80) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.78) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.76) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.74) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.72) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.70) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.68) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.66) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.65) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.63) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.61) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.60) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.58) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.57) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.56) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.54) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.53) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.52) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.52) } 

    add EfficiencyFormula {321} {-11} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.03) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.05) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.06) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.07) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.08) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.09) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.10) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.11) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.12) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.13) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.14) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.15) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.16) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.17) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.19) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.20) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.21) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.21) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.22) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.24) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.24) } 

    add EfficiencyFormula {321} {211} { (eta<-3.5 || (-1 < eta && eta < 1) || eta>2.0)*( 0.00 ) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 0.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (0.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 1.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (1.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 2.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (2.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 3.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (3.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 4.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (4.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 5.90) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (5.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.00) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.10) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.20) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.30) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.40) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.50) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.60) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.70) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.80) * (0.00) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 6.90) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (6.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.00) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.10) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.20) * (0.01) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.30) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.40) * (0.02) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.50) * (0.03) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.60) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.70) * (0.04) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.80) * (0.05) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 7.90) * (0.06) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (7.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.00) * (0.07) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.10) * (0.08) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.20) * (0.09) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.30) * (0.10) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.40) * (0.11) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.50) * (0.12) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.60) * (0.13) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.70) * (0.14) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.80) * (0.15) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 8.90) * (0.16) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (8.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.00) * (0.17) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.00 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.10) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.10 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.20) * (0.18) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.20 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.30) * (0.19) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.30 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.40) * (0.20) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.40 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.50) * (0.21) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.50 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.60) * (0.21) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.60 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.70) * (0.22) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.70 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.80) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.80 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 9.90) * (0.23) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (9.90 <= (pt * cosh(eta)) && (pt * cosh(eta)) < 10.00) * (0.24) + 
                                     ((-3.5 <= eta && eta <= -1.0) || (1.0 <= eta && eta <= 2.0)) * (10.00 <= (pt * cosh(eta))) * (0.24) } 


 # efficiency for other charged particles (should be always {0} {0} {formula})

    add EfficiencyFormula {0} {0}     {      (fabs(eta) < 4.0) * (0.00)     }

}




##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle

  add Branch TrackSmearing/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch mRICHPID/tracks mRICHTrack Track

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET
 
  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon
  
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}


