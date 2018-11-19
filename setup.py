from distutils.core import setup

setup(name='vkmz',
      version='1.4dev',
      python_requires='>=3.6',
      description='metabolomics formula prediction and van Krevelen diagram generation',
      author='Mark Esler',
      author_email='eslerm@umn.edu',
      url='https://github.com/HegemanLab/VKMZ',
      packages=['distutils', 'distutils.command'],
     )
