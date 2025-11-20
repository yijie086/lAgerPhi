# l/A-event Generator
This is the Argonne generic l/A-event generator (`lAger`), a flexible MC generator system
to simulate electro- and photo-production off nucleons and nuclei.

Below you can find an overview of the release versions, as well as a short tutorial and
copyright notice. If you use lAger to generate data used in a presentation or an article
in a scientific publication, please cite:

_S. Joosten, Argonne l/A-event Generator (2021), GitLab repository,
https://eicweb.phy.anl.gov/monte_carlo/lager_

# Versions
* v3.6.x Tweaked VMD formalism and better examples for DVMP at EIC
* v3.5.x New container/installer setup
* v3.4.x Add support for psi prime
* v3.3.x Adds photoproduction generators based on JPAC pomeron and baryon resonance calculations
* v3.2.x Adds HepMC3 and fermi momentum and misc fixes
* v3.1.x First stable relase of `lAger`

# Tutorial

Simple Installation
------------
1. Create a local directory that you want to work in, e.g `$HOME/lager`, and go into this
   directory.
```bash
mkdir $HOME/lager
cd $HOME/lager
```

2. Execute the following line in your terminal to setup your environment in this directory
   to install the latest stable container
```bash
curl https://eicweb.phy.anl.gov/monte_carlo/lager/-/raw/master/install.sh | bash
```

3. The launchers are installed in the `bin` directory. 
   You can start `lager` by using the installed `bin/lager` launcher script. 
   There is also a `lager-shell` launcher that opens a shell within the
   singularity container, useful for expert usage. Finally, some example configuration
   files are downloaded into the `share/lager-examples` directory.

4. You can now explore the command line flags for `lager` using the `-h` flag.
```bash
lager -h # if you have <PREFIX>/bin in your $PATH
<PREFIX>/bin/lager -h # if you do not have the prefix in your $PATH
```

## Command-line interface

These instructions assume you are working from the singularity container, but mostly also
apply if you are using the docker container or a local build.

The main command line options for `lager` are:
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

`lager` is mainly configured through a JSON file. You can find examples under
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
These are the 6 main fields to setup `lAger`:
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
The main appeal of `lager` lies into the flexibility of the generator as it is split in smaller
sub-generators that can be mixed and matched together to form a very powerful generation
tool. Below you can see the generator setup for a solid-ep simulation performed for the
SoLID pCDR document.
```python
"generator" : {
  "type" : "ep-2gluon",
  "vertex" : {"type" : "linear", "range" : [ "-7.5", "7.5" ]},
  "beam" : {
    "lepton": {
      "type" : "constant",
      "particle_type" : "e-",
      "dir" : [ "0", "0", "1" ],
      "energy" : "11.0"
    }, 
    "ion": {
      "type" : "constant",
      "particle_type" : "proton",
      "dir" : [ "0", "0", "-1" ],
      "energy" : "0.9382721"
    }
  },
  "target" : {"type": "primary"},
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
3. `beam`: Lepton and ion beam definitions. For a normal "constant" beam, it requires
   particle type, beam 3-vector and beam energy.
4. `target`: Actual target we use. In this case we use the primary ion beam as target
   (`primary`).
5. `photon`: Real or virtual photon intensity. Here it is a virtual photon with y between 0.6 and 1. Can be set to "flat" in case you do not want to fold in a photon intensity. For some processes (e.g. DVCS/BH) you may want to work around the photon.
6. `process_0`: Your process, can be any of the supported processes. Most processes
   currently give sensible results for both electro- and photo-production.
   You can in principle use `process_0` through `process_9` but make sure not to specify
   more than a single process at a time - process mixing is not (yet) supported.

## Running an example
To run an example using the `solid.ep-2gluon.json` configuration, go into the `examples`
directory and execute `lAger` (make sure to fill in your desired output directory).
```bash
cd examples
lager -c solid.ep-2gluon.json -r 1 -o <YOUR_OUTPUT_DIRECTORY>
```

Please use the issue tracker on https://eicweb.phy.anl.gov/monte_carlo/lager to file bug
reports.

# Copyright

`lAger`: General Purpose l/A-event Generator
Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>

`lAger` is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Shoftware Foundation, either version 3 of the License, or
(at your option) any later version.

`lAger` is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with `lAger`.  If not, see <https://www.gnu.org/licenses/>.



HK May 27, 2024:

To compile and have changes propagate:

Edit .cc or .hh files using a non-lager-shell terminal. 

From lager-shell, and in the directory /w/hallc-scshelf2102/jpsi-007/klest/lager2, run "cmake --build build -- install -j10"

Add to LD_LIBRARY_PATH (still in lager-shell) with: export LD_LIBRARY_PATH="/work/hallc/jpsi-007/klest/lager2/lib/:$LD_LIBRARY_PATH"

Then run for example: "./bin/lager -c jsons/hatta.ep-phi.json -o outputs -r 9999"

