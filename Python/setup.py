#!/usr/bin/env python

from distutils.core import setup

setup(name='PyCIL',
      version='0.1',
      description='python scripts for CIL project on multi-spectral CT processing',
      author='Christian Kehl',
      author_email='chke@dtu.dk',
      download_url='https://gitlab.com/CKehl/MECT',
      packages=['tools'],
      package_dir={'tools' : 'src/tools'},
      py_modules=[],
      scripts=['src/FAMStoHDF5.py', 'src/HDF5toFAMS.py', 'src/HDF5toMatV7.py', 
               'src/HDF5toMHD.py', 'src/MatV7toHDF5.py', 'src/MatV7toMHD.py', 
               'src/MHDtoHDF5.py', 'src/runSegmentation_SFAMS.py']
     )
