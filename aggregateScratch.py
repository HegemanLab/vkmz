import extractNeededElementalData
import processElementalData
import bmrbLookup as bmrb
import argparse
import os
from process_mzs_mzML import process_mzs as ML_process
from MzXML import MzXML
from process_mzs import process_mzs as XML_process
from flexPlot import plotVanK

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
# setup for multiprocessing
if vkMultiprocessing:
  import multiprocessing
  vkCores = multiprocessing.cpu_count()
  print(vkCores)

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
  # get lookup table.
  lt = bmrb.getLookupTable('bmrb-db2.csv')
  # set up. elements could be changed but would need to do some editing elsewhere.
  elements = ['C', 'H', 'O', 'N']
  if polarity == 'both':
    buildRatios('pos', vkInputMzs) 
    buildRatios('neg', vkInputMzs) 
  if polarity == 'pos':
    pos_comps = list()
    for mz in vkInputMzs[1]:
      pos_comps.append(bmrb.getFormulaFromMass(bmrb.adjust(mz, 'pos'), lt, tolerance=vkThreshold)) # adjust mass and search in lookup table. Store result in list.
    pos_comps = filter(lambda a: a != 'No Match', pos_comps) # Filter out no matches
    pos_elements = extractNeededElementalData.find_elements_values(elements_to_find=elements, compounds=pos_comps) # Get elements from compounds
    pos_ratios = processElementalData.process_elemental_data(pos_elements) # Turn elements into ratios
    if vkOutput:
      pos_filename = 'example-ratios-pos.csv'  # forcing file name. FLAG
      with open(pos_filename, 'w') as f: 
        for ratio in pos_ratios:
          f.writelines(str(ratio).strip('[]') + '\n')
    for vkType in vkPlotTypes:
      vkPlotter(pos_ratios, vkType)
  if polarity == 'neg':
    neg_comps = list()
    for mz in vkInputMzs[0]:
      neg_comps.append(bmrb.getFormulaFromMass(bmrb.adjust(mz, 'neg'), lt, tolerance=vkThreshold)) # adjust mass and search in lookup table. Store result in list.
    neg_comps = filter(lambda a: a != 'No Match', neg_comps) # Filter out no matches
    neg_elements = extractNeededElementalData.find_elements_values(elements_to_find=elements, compounds=neg_comps) # Get elements from compounds
    neg_ratios = processElementalData.process_elemental_data(neg_elements) # Turn elements into ratios
    if vkOutput: 
      dataOutput(neg_ratios, polarity)
    for vkType in vkPlotTypes:
      vkPlotter(neg_ratios, vkType)

# write vk ratios as csv file
def dataOutput(ratios, polarity):
  try:
    filename = 'example-ratios-' + str(polarity) + 'neg.csv'  # forcing file name. FLAG
    with open(filename, 'w') as f: 
      for ratio in ratios:
        f.writelines(str(ratio).strip('[]') + '\n')
  except ValueError:
    print('"%s" could not be saved.' % vkLoad)

# load vk ratio csv file(s)
# FLAG: add multiple file support
def loadRatios(vkLoad, vkPlotTypes):
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
    for vkType in vkPlotTypes:
      vkPlotter(ratios, vkType)
  except ValueError:
    print('The %s data file could not be loaded.' % vkLoad)

def vkPlotter(ratios, type):
  plotVanK(ratiosList=ratios, typeOfPlot=type)

if vkLoad != '':
  loadRatios(vkLoad, vkPlotTypes)
else:
  buildRatios(vkPolarity, dataParser(vkInput, vkThreshold))

