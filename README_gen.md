# Phi analysis generator

## Installation 

```
git clone git@github.com:yijie086/lAgerPhi.git
```

```
cd lAgerPhi
```

```
./install.sh
```

And also compile a tool to clean lund file:

```
g++ -std=c++17 -O2 clean_lund.cc -o clean_lund
```



## Run code for phi simulation gemc generator

```
bin/lager -c CLAS1210GeV.ep-phi.gen.json -e 5000 -r 001 -o .
```

you can change the setting in the `CLAS1210GeV.ep-phi.gen.json`

Run number and random seed `-r`

event max number `-e`

output folder `-o`

However, this code has ~5% possibility to generate `nan` value, so clean the lund file is needed, use `clean_lund.cc`, this program also change `.gemc` to `.txt` to support OSG requirements.

```
./clean_lund CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000.gemc CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000_clean.txt
```

Then one can use the `lundANA.cpp` to check the lund file's quality.

```
root -l -b -q 'lundANA.cpp("CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000_clean.txt","events.root")'
```

when you check the root file, please notic the scale of the x axis, for the exclusive variables, they are super small.

## Shell macro for phi generator

if you want to generate lund file one by one, use

```
./run_lund.sh 
```

you can change the config in the shell macro :

```
# JSON configuration file for lAgerPhi
CONFIG="CLAS1210GeV.ep-phi.gen.json"

# Number of events per run (defines the final "-1000" part)
NEV=5000

# Total number of runs to generate
NRUNS=10
```

if you want to generate lund file in parallel, use

```
./run_lund_large.sh
```

you can change the config in the shell macro:

```
# JSON configuration file for lAgerPhi
CONFIG="CLAS1210GeV.ep-phi.gen.json"

# Number of events per run (defines the final "-5000" part in filenames)
NEV=5000

# Total number of runs you want to generate in total (for launcher mode)
TOTAL_RUNS=100

# Maximum number of runs per background job
MAX_CHUNK=10
```

