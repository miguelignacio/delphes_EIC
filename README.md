# delphes_EIC

Install Delphes3 following:
https://github.com/delphes/delphes

The detector card contains an EIC detector based on the EIC detector handbook v1.2
http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf

So far it incorporates tracking, EMCAL and HCAL but lacks implementation of PID (it can be done though, following the LHCb card example)

Magnetic field: 1.5 T,Solenoid length: 2.0 m, Tracker radius: 80 cm. 

You can run Pythia8 within Delphes. The command file shown here is suitable for DIS at EIC. 

Run generation command:
`./DelphesPythia8 cards/delphes_card_EIC.tcl examples/Pythia8/DIS.cmnd out.root`

You can see examples of analysis code in the Delphes page above

Run visualization command:
 `root -l examples/EventDisplay.C'("cards/delphes_card_EIC.tcl","out.root")'`
 
The two examples shown here are for neutral-current and charged-current event 
for beam energies of 10 GeV electron on 100 GeV proton (63 GeV center-of-mass energy). 
