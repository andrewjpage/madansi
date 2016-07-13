import os

class ValidateOutputArguments(object):
    
    def __init__(self, *args):
        pass
    
    def run(self, *args):
        for arg in args:
            if  os.path.isfile(arg):
                raise  ValueError("Output file already exists. Cannot continue.")


