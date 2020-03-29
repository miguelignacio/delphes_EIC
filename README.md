# delphes_EIC

Install Delphes3 following:
https://github.com/delphes/delphes

You can run Pythia8 within Delphes. The command shown here is suitable for DIS at EIC. 
The detector card contains an EIC detector based on the EIC detector handbook v1.2
http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf
So far it incorporates tracking, emcal and hcal but lacks implementation of PID (it can be done though, following the LHCb card example)

Run generation command:
./DelphesPythia8 cards/delphes_card_EIC.tcl examples/Pythia8/DIS.cmnd out.root

You can see examples of analysis code in the Delphes code above. 

Run visualization command:
 root -l examples/EventDisplay.C'("cards/delphes_card_EIC.tcl","out.root")'
