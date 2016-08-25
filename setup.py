import glob
from setuptools import setup, find_packages

setup(
	name= 'madansi',
	version= '1.0.1',
	description= 'Reorders and orientates contigs from a fasta file using a pan genome graph.',
	packages = find_packages(),
	author = 'Joanna Tumelty',
	author_email= 'jt20@sanger.ac.uk',
	url='https://github.com/JTumelty/madansi',
	scripts=glob.glob('scripts/*'),
	test_suite='nose.collector',
	tests_require=['nose >= 1.3'],
    install_requires=[
        'networkx >= 1.11',
        'biopython'
                ],
	license='GPLv3',
	classifiers=[
	        "License :: OSI Approved :: GNU General Public License (GPLv3)",
	        "Programming Language :: Python",
	        "Development Status :: 4 - Beta",
	        "Intended Audience :: Science/Research",
	        "Topic :: Scientific/Engineering :: Bio-Informatics",
	        ],	
	)
