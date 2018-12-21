Old version of J/psi and Pc Monte Carlo generator.
v1.x was used for the pentaquark proposal
v2.x is capable of simulating collider experiments, was used for first
EIC/SoLID projections

# Basic Example
1. Ensure the generator is installed and accessible from your path
2. Use the file examples/upsilon-ep-e10p100.json as configuration (included in the source tree)
3. Choose an output directory, let's call it OUTDIR `-o OUTDIR`
4. We want the run number (random seed) 1 for this run `-r 1`
5. To run, do 
```bash
pcsim -c upsilon-ep-e10p100.json -o OUTDIR -r1
```
6. That's all!

# Configuration file explanation
Note that I only documented those settings that are meant to be changed
during normal operation.
## General settings
```json
"run"   : "Run number (random seed), overwritten by the command line",
"events": "Number of events to be generated if no target luminosity given",
"lumi"  : "Desired luminosity (in fb^-1), trumps the number of events requested",
"tag"   : "Unique label attached to the output files, e.g. beam energies",
"info"  : "Optional info about this run"
```
## Electron beam ('beam') settings
```json
"beam": {
   "dir"   : "Electron beam direction (e.g. [0, 0, -1])",
   "energy": "Electron beam energy in GeV"
}
```
## Proton beam ('target') settings
```json
"target": {
   "dir"   : "Electron beam direction (e.g. [0, 0, 1])",
   "energy": "Electron beam energy in GeV"
}
```
## Secondary photon beam settings
Do not change for electro-production (defaults should always work).
## Process generator 
### Oleksii's model for VM production
```json
"process_1": {
    "type": "oleksii_2vmp",
    "vm_type": "J/psi (443) or Upsilon (553)",
    "T0": "Subtraction constant, related to binding. Vary between 0: no binding, 4: best fit and 8: strong binding"
}
```
