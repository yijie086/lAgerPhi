#!/usr/bin/python -O

class classdict(dict):
    def __init__(self, **kvargs): 
        dict.__init__(self, **kvargs)
    def __getattr__(self, key): 
        ## redirect magic methods to the base
        if key.startswith('__') and key.endswith('__'):
            return dict.__getattr__(key)
        return self[key]
    def __setattr__(self, key, val): 
        self[key] = val

class autovivi_classdict(classdict):
    def __init__(self, **kvargs):
        classdict.__init__(self, **kvargs)
    def __getattr__(self, key): 
        ## redirect magic methods to the base
        if key.startswith('__') and key.endswith('__'):
            return dict.__getattr__(key)
        ## autovivification 
        if not key in self: 
            self[key] = autovivi_classdict()
        return self[key]
