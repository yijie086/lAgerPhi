# Lepton-Ion Event Generator (liege)
Generic event generator to simulate electro- and photo-production off nucleons and nuclei.

# Versions
v3.0.0 Stable release version

## Old (depecrated)
Liege used to be known as pcsim.
v2.x (called is capable of simulating collider experiments, was used for first EIC/SoLID projections)
v1.x was used for the Hall C pentaquark proposal

# Tutorial
## Setup of the liege singularity container on your system:
The default mode to run the generator is through singularity. To setup the generator on
your system, first ensure singularity is installed. Then follow these instructions:
1. Clone this repository and checkout the desired stable release (e.g. v3.0.0)
```bash
git clone https://eicweb.phy.anl.gov/monte_carlo/liege.git
cd liege && git checkout v3.0.0
```
2. Run `deploy.py` script to install the container to a prefix of your choice, e.g. `$HOME/local/opt/liege`.
```bash
./deploy.py $HOME/local/opt/liege
```
There are several flags supported by `deploy.py`. To get a full overview you can run `deploy.py -h`.
Noteable flags:
  - Print all available options: `-h`.
  - Add a custom bind path (e.g. /site) to the launcher: `-b /site`.
  - Deploy a different version (e.g. master): `-v master`.

3. If this all executed without issues you should now be able to run the generator from
   the prefix. You can either refer directly to the launcher under `<PREFIX>/bin/liege` or
   ideally add `<PREFIX>/bin` to your `$PATH` in your `bashrc` (or equivalent).

4. You can now explore the command line flags for `liege` using the `-h` flag.
```bash
liege -h # if you have <PREFIX>/bin in your $PATH
<PREFIX>/bin/liege -h # if you do not have the prefix in your $PATH
```
From now on, for simplicity it will be assumed you have `<PREFIX>/bin` in your `$PATH`. If
you do not, substitute `liege` for the explicit path `<PREFIX>/bin/liege`.

## Command-line interface

These instructions assume you are working from the singularity container, but mostly also
apply if you are using the docker container or a local build.

The main command line options for `liege` are:
  - configuration json file: `-c your_configuration.json`
  - run number (also random seed, e.g. 1): `-r 1`
  - output directory: `-o $HOME/some_output_directory/`

The generator will write 3 output files for each run into this directory. 
  - A ROOT file with the generator output
  - A log file with the terminal output from the generator run. This will contain have useful
    statistics such as total number of generated events and total integrated cross
    section.
  - A copy of the input JSON file with optional overrides from the command line added in,
    to enable easy determination of the exact configration used for a specific generated
    output file.
  
The file names will contain the run number and desired luminosity/number of events to generate, so it is
safe the use the same output directory when running a large amount of processes in
parallel.

## JSON configuration file

`liege` is mainly configured through a JSON file. You can find examples under
 `examples`. The entire configuration is placed under the `mc` key. As example, let's
 discuss `solid.ep-2gluon.json`.

### Top-level configuration
```python
  "type" : "solid",
  "tag"  : "",
  "lumi" : "1000",
  "generator" : {},
  "detector": {},
  "reconstruction: {}
```
These are the 6 main fields to setup `liege`:
1. `type`: This can be any name you want to use to identify this simulation. Here we use
   `solid` to identify this as a SoLID simulation. The first part of the output file name
   will be this type key.
2. `tag`: An additional tag to identify this simulation. Usefull to distinguish multiple
   very similar configurations, e.g. when changing beam energies.
3. `lumi` or `events`: Either the desired luminosity (in fb^-1), or the desired number of
   events. The `lumi` key is given precedence to the `events` key in case both are
   specified.
4. `generator`: The actual generator configuration, the most important component. 
5. `detector`: Optional simple geometric acceptance components barrel, spectrometer, composite
   (multiple barrels/spectrometers) or the default null/4pi (detect everything).
6. `reconstruction`: Optional requirement that certain particles were detected. Will only
   write out events that fit the reconstruction requirements. 

### Generator configuration
The main appeal of `liege` lies into the flexibility of the generator as it is split in smaller
sub-generators that can be mixed and matched together to form a very powerful generation
tool. Below you can see the generator setup for a solid-ep simulation performed for the
SoLID pCDR document.
```python
"generator" : {
  "type" : "ep-2gluon",
  "vertex" : {"type" : "linear", "range" : [ "-7.5", "7.5" ]},
  "beam" : {
    "type" : "primary",
    "particle_type" : "11",
    "dir" : [ "0", "0", "1" ],
    "energy" : "11.0"
  },
  "target" : {
    "type" : "primary",
    "particle_type" : "2212",
    "dir" : [ "0", "0", "-1" ],
    "energy" : "0.9382721"
  },
  "photon" : {"type" : "vphoton", "y_range" : [ "0.6", "1" ]},
  "process_0" : {
    "type" : "brodsky_2vmX"
    ...
  }
}
```
The generator is modular, and expects 6 keys to be present:
1. `type`: some name you want to use to refer to this combination of detector elements
2. `vertex`: your vertex definition, e.g. `linear` or `origin`
3. `beam`: lepton beam definition. For a `primary` beam it requires particle type, beam
   3-vector direction and beam energy
4. `target`: ion beam definition. Primary beam definition is identical to the lepton beam,
   but more complicated beams are also supported (e.g. nucleon in a nucleus).
5. `photon`: Real or virtual photon intensity. Here it is a virtual photon with y between 0.6 and 1. Can be set to "flat" in case you do not want to fold in a photon intensity. For some processes (e.g. DVCS/BH) you may want to work around the photon.
6. `process_0`: Your process, can be any of the supported processes. Most processes
   currently give sensible results for both electro- and photo-production.
   You can in principle use `process_0` through `process_9` but make sure not to specify
   more than a single process at a time - process mixing is not (yet) supported.

## Running an example
To run an example using the `solid.ep-2gluon.json` configuration, go into the `examples`
directory and execute `liege` (make sure to fill in your desired output directory).
```bash
cd examples
liege -c solid.ep-2gluon.json -r 1 -o <YOUR_OUTPUT_DIRECTORY>
```

Please use the issue tracker on https://eicweb.phy.anl.gov/monte_carlo/liege to file bug
reports.
