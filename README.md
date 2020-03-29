# delphes_EIC

Install Delphes3 following:
https://github.com/delphes/delphes

You can run Pythia8 in it, with the command that is suitable for DIS at EIC. 
The detector card contains an EIC detector based on the EIC detector handbook v1.2

Run command:
./DelphesPythia8 cards/delphes_card_EIC.tcl examples/Pythia8/DIS.cmnd out.root

Example of analysis code:


Visualize as:
 root -l examples/EventDisplay.C'("cards/delphes_card_EIC.tcl","out.root")'
