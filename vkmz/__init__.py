#!/usr/bin/env/python
"""metabololite prediction and plotting tool

VKMZ predicts molecular formulas by searching a known mass-formula dictionary
for a feature observed by a mass spectrometer. Elemental ratios forpredicted-features
are calculated to create the carbon-to-oxygen and carbon-to-hydrogen axis of a
van Krevelen Diagram (VKD). VKD's are a convenient visualization tool for
briefly conveying the constituence of a complex MS mixture (e.g., untargetted
plant metabolomics). As output predicted-feature are saved to a tabular file,
an interactive VKD web page, and other optional formats.
"""

# NOTE: reviewers with Python Packaging experience
#       How do I clean my namespace so that pydoc works on all modules and so
#       that vkmz is ready for API importation (__all__)?
#       Wrapping things with if __name__ == __main__ has consequences that may
#       need to be avoided.
