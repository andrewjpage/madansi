import os

class ValidateInputArguments(object):
    
    def __init__(self, *args):
        pass
    
    def run(self, *args):
        for arg in args:
            if not os.path.isfile(arg):
                raise  ValueError("Input file does not exist. Cannot continue.")
