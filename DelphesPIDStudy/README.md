# Setup Your Environment

Make sure that the environment variable ```DELPHES_ROOT``` points to the top-level directory containing your local Delphes source code.

# Build

```make -j```

# Run

Using a file, ```out.root```, produced by running Delphes using the TCL file ```delphes_EICmatrixv2_3T.tcl```:

```
./DelphesPIDStudy.exe -i 'out.root' -o pidstudy.root
```

Efficiency plots for electrons, protons, kaons, and pions are now stored in pidstudy.root.
