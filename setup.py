#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='photred',
      version='1.0.0',
      description='Generic PSF Photometry Pipeline',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/PHOTRED',
      packages=['photred'],
      #scripts=['bin/photred'],
      install_requires=['numpy','astropy(>=4.0)','scipy','dlnpyutils(>=1.0.3)']
      #include_package_data=True,
)
