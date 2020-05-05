source delphes_card_EIC.tcl

#################################
# Alter the Solenoid B-Field
#################################
module ParticlePropagator ParticlePropagator {
    # magnetic field
    set Bz PARAM_BFIELD
}

