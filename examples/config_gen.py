#!/usr/bin/env python

import json

# todo: add argparse

with open('upsilon-ep-e10p100.json') as f:                                                                                                                              
    data = json.load(f) 


data["mc"]["generator"]["beam"]["energy"]   = 8.0
data["mc"]["generator"]["target"]["energy"] = 120.0

print json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))

