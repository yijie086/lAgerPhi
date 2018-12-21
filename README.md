Old version of J/psi and Pc Monte Carlo generator.
v1.x was used for the pentaquark proposal
v2.x is capable of simulating collider experiments, was used for first
EIC/SoLID projections

# Basic Example
1. Ensure the generator is installed and accessible from your path
2. Use the file examples/upsilon-ep-e10p100.json as configuration (included in the source tree)
3. Cheese and output directory, let's call it OUTDIR ```-o OUTDIR```
4. We want the run number (random seed) 1 for this run ```-r 1```
5. 4. To run, do 
```bash
pcsim -c upsilon-ep-e10p100.json -o OUTDIR -r1
```
6. That's all!