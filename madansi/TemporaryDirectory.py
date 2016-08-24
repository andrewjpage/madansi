import tempfile
import shutil

class TemporaryDirectory(object):
    def create_temporary_directory(self):
        """Creates a temporary directory to store files in"""
        temp_dir = tempfile.mkdtemp()
        return temp_dir
    
    def remove_temporary_directory(self, temp_dir):
        shutil.rmtree(temp_dir)