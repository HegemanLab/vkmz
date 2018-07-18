"""
setup.py
~~~~~~~~~~~~
:copyright: (c) 2018 Hegeman Lab, University of Minnesota
:license: MIT
"""

from setuptools import setup

setup(
    name="VKMZ",
    version="1.3.0",
    description="metabolomics formula prediction and van Krevelen diagram generation",
    author="Mark Esler",
    author_email="eslerm@umn.edu",
    url="https://github.com/HegemanLab/vkmz",
    license="MIT",
    packages=["vkmz"],
    entry_points={
      'console_scripts': [
        'vkmz = vkmz.vkmz'
      ]
    }
)
