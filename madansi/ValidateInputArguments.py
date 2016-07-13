import os

class ValidateInputArguments(object):
    
    def __init__(self, *args):
        pass
    
    def run(self, *args):
        for arg in args:
            if not os.path.isfile(arg):
                raise  ValueError("Input file does not exist. Cannot continue.")


  #def __init__(self,inputfile,outputfile):
  #    self.inputfile = inputfile
  #    self.outputfile = outputfile
  #
  #def run(self):
  #    
  #    if not os.path.isfile(self.inputfile):
  #        raise ValueError("Input file does not exist. Cannot continue.")
  #    else:
  #        if os.path.isfile(self.outputfile):
  #            raise ValueError("Output file already exists. Cannot continue.")