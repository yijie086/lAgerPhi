#!/usr/bin/env python

import json
import argparse


def main():

    parser = argparse.ArgumentParser(description='config generator for pcsim.')
    #parser.add_argument('beam', metavar='N', type=int, nargs='+',
    #                    help='an integer for the accumulator')
    #parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                    const=sum, default=max,
    #                    help='sum the integers (default: find the max)')
    
    parser.add_argument('-i', action="store", dest="input_file", help="input json file. Which sets default values. Other commandline options override these configuration settings.")
    parser.add_argument('--e-beam-energy', action="store", dest="e_energy", help="electron beam energy in GeV")
    parser.add_argument('--ion-beam-energy', action="store", dest="ion_energy", help="ion beam energy in GeV")


    args = parser.parse_args()

    in_file = 'upsilon-ep-e10p100.json'
    #if args.input_file is not None:
    #    in_file =  args.input_file

    with open(in_file) as f:
        data = json.load(f) 
 
    if args.e_energy is not None:
        #print(str("electron beam energy is  {} ").format(args.e_energy))
        data["mc"]["generator"]["beam"]["energy"]   = float(args.e_energy)
    
    if args.ion_energy is not None:
        #print(str("electron beam energy is  {} ").format(args.e_energy))
        data["mc"]["generator"]["target"]["energy"]   = float(args.ion_energy)
    
    #data["mc"]["generator"]["beam"]["energy"]   = 8.0
    #data["mc"]["generator"]["target"]["energy"] = 120.0
    data["mc"]["lumi"] = 1.0
    
    print json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))

if __name__ == '__main__':
    main()
