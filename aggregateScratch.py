import os
import time
import extractNeededElementalData
import processElementalData
import bmrbLookup as bmrb
import argparse
from process_mzs_mzML import process_mzs as ML_process
from MzXML import MzXML
from process_mzs import process_mzs as XML_process
from flexPlot import plotVanK
import multiprocessing
from multiprocessing import Pool
from functools import partial

parser = argparse.ArgumentParser()
parser.add_argument('--load',      '-l', nargs='?', default='',     help='Load a previously generated ratio table. Set file path. Disabled by default.')
parser.add_argument('--input',     '-i', nargs='*', default='',     help='Enter full mzXML/mzML file paths. For multiple files seperate files with a space.')
parser.add_argument('--output',    '-o', action='store_true',       help='Call variable to ouput a ratio file.')
parser.add_argument('--polarity',  '-p', nargs='?', default='both', type=str, choices=['both', 'pos', 'neg'], help='Set to "pos", "neg", or "both". Default is "both". Only one plot type can be set.')
parser.add_argument('--threshold', '-t', nargs='?', default='10',   type=int, help='Set threshold as a percent integer from 0 to 100. eg. 15 for 15%%.')
parser.add_argument('--plottype', '-pt', nargs='*', default=['scatter'], choices=['scatter', 'heatmap', '3d'], help='Set to "scatter", "heatmap", or "3d". Default is "scatter".')
parser.add_argument('--multiprocessing', '-m', action='store_true', help='Call variable to use multiprocessing. One process per core.')
args = parser.parse_args()

# read load argument
vkLoad = getattr(args, "load")

# read input argument(s)
vkInput = getattr(args, "input")

# read output argument
vkOutput = getattr(args, "output")

# read polarity argument
vkPolarity = getattr(args, 'polarity')

# read threshold argument
# FLAG: a dynamic threshold function has been made, implement
vkThreshold = getattr(args, "threshold")
if 0 <= vkThreshold <= 100:
  vkThreshold = vkThreshold * 0.01
else:
  raise ValueError("The given threshold, %i, is out of bounds." % (vkThreshold))
  exit()

# read plottype argument
vkPlotTypes = getattr(args, 'plottype')

# read multiprocessing argument
vkMultiprocessing = getattr(args, "multiprocessing")

def dataParser(vkInput, vkThreshold):
  # load input files
  vkInputFiles = []
  for vkFile in vkInput:
    if vkFile.lower().endswith(('.mzml', '.mzxml')):        # verify files extension of --input
      vkInputFiles.append(vkFile)
    else:
      raise ValueError('Input file "%s" is not a mzML or mzXML file.' % (vkFile))
      exit()
  # parse input files
  vkInputMzs = [[],[]]
  for vkFile in vkInputFiles:
    if vkFile.lower().endswith('.mzxml'):
      mzXML = MzXML()
      mzXML.parse_file(vkFile)
      vkInputMzsTemp = XML_process(mzXML, threshold=vkThreshold)
      vkInputMzs[0] = vkInputMzs[0] + vkInputMzsTemp[0]
      vkInputMzs[1] = vkInputMzs[1] + vkInputMzsTemp[1]
    elif vkFile.lower().endswith('.mzxml'):   #FLAG, test this filetype
      vkInputMzs = ML_process(f, threshold=vkThreshold)
  # Removes all duplicates from both neg and pos lists
  vkInputMzs[0] = list(set(vkInputMzs[0]))
  vkInputMzs[1] = list(set(vkInputMzs[1]))
  return vkInputMzs

# vk ratio dataset builder for each polarity
# dataset is list with four elements.
# input mass to charge ratios checked against database of masses with known chemical structure
# each identified mass is saved across the four elements, at a given index, of the dataset list
# the first  element represents ratio of hydrogens in identified masses structure
# the second element represents ratio of carbons                    #FLAG VERIRFY
# the third  element represents if a nitrogen is in the strcutre    #FLAG remove element, have plotter check last element != 0
# the fourth element represents ratio of nitrogens
def buildRatios(polarity, vkInputMzs):
  lt = bmrb.getLookupTable('bmrb-db2.csv')
  # get lookup table.
  # set up. elements could be changed but would need to do some editing elsewhere.
  elements = ['C', 'H', 'O', 'N']
  if polarity == 'both':
    buildRatios('pos', vkInputMzs) 
    buildRatios('neg', vkInputMzs) 
  if polarity == 'pos' or polarity == 'neg':
    index = int(polarity == 'pos') # int(True) == 1
    if vkMultiprocessing:
      pool = Pool()
      multiprocessMzsArgs = partial(multiprocessMzs, polarity, lt, vkThreshold)
      identified = pool.map(multiprocessMzsArgs, vkInputMzs[index])
      pool.close()
      pool.join()
    else:
      identified = []
      for mz in vkInputMzs[index]: # pos2 in index 1, neg in index 0
        identified.append(bmrb.getFormulaFromMass(bmrb.adjust(mz, str(polarity)), lt, tolerance=vkThreshold)) # adjust mass and search in lookup table. Store result in list.
    identified = filter(lambda a: a != 'No Match', identified) # Filter out no matches
    identifiedElements = extractNeededElementalData.find_elements_values(elements_to_find=elements, compounds=identified) # Get elements from compounds
    identifiedRatios = processElementalData.process_elemental_data(identifiedElements) # Turn elements into ratios
    if vkOutput: 
      saveRatios(identifiedRatios, polarity)
    for type in vkPlotTypes:
      plotRatios(identifiedRatios, type)

def multiprocessMzs(polarity, lt, vkThreshold, inputMz): # recieves a single Mz
  return bmrb.getFormulaFromMass(bmrb.adjust(inputMz, str(polarity)), lt, tolerance=vkThreshold)

# write vk ratios as csv file
def saveRatios(ratios, polarity):
  try:
    filename = 'ratios-' + time.strftime("%Y%m%d%H%M%S-") + str(polarity) + '.csv'
    with open(filename, 'w') as f: 
      for ratio in ratios:
        f.writelines(str(ratio).strip('[]') + '\n')
  except ValueError:
    print('"%s" could not be saved.' % vkLoad)

# load VK ratio csv file
# FLAG: add multiple file support
def loadRatios(vkLoad):
  try:
    # read in ratios
    with open(vkLoad, 'r') as f:
      ratios = f.readlines()
    # split ratios into proper value
    for i in range(0, len(ratios)):
      ratios[i] = ratios[i].split(', ')
    # Cast to correct type
    ratios[0] = map(lambda x: float(x), ratios[0])
    ratios[1] = map(lambda x: float(x), ratios[1])
    ratios[3] = map(lambda x: float(x), ratios[3])
    for type in vkPlotTypes:
      plotRatios(ratios, type)
  except ValueError:
    print('The %s data file could not be loaded.' % vkLoad)

def plotRatios(ratios, type):
  plotVanK(ratiosList=ratios, typeOfPlot=type)

if vkLoad != '':
  loadRatios(vkLoad)
else:
  buildRatios(vkPolarity, dataParser(vkInput, vkThreshold))

