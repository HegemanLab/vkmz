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
parser.add_argument('--output',    '-o', action="store_true", default='False', help='Call variable to ouput a ratio file.')
parser.add_argument('--threshold', '-t', nargs='?', default='10',   type=int, help='Set threshold as a percent integer from 0 to 100. eg. 15 for 15%%.')
parser.add_argument('--polarity',  '-p', nargs='?', default='both', type=str, help='Set to "pos", "neg", or "both". Default is "both".')
parser.add_argument('--plottype', '-pt', nargs='*', default='[scatter]', help='Set to "scatter", "heatmap", or "3d". Default is "scatter".')
args = parser.parse_args()

# exits script if a load argument is used
vkLoad = getattr(args, "load")
if vkLoad != '':    # FLAG: not sanitation of inpute filename
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
    ratios[2] = map(lambda x: float(x), ratios[2])
  except ValueError:
    print('The %s data file could not be loaded.' % vkLoad)
  exit()    # vkLoad ~ should be ready to be piped to ploter

# code beyond this point assumes no load argument was used
# make a list from input files argument
vkInput = getattr(args, "input")
if vkInput == '':       # --input cannot be empty
  raise ValueError("No input file(s) specified. See -h for help menu.")
  exit()
else:
  vkInputFiles = []
  for vkFile in vkInput:
    if vkFile.lower().endswith(('.mzml', '.mzxml')):        # verify files extension of --input
      vkInputFiles.append(vkFile)
    else:
      raise ValueError('Input was set to "%s". It must be set to the filepath of a mzML or mzXML file.' % (vkFile))
      exit()

# set output boolean
vkOutput = getattr(args, "output")
if vkOutput == True:
  vkOutput = 'y'
else:
  vkOutput = 'n'
write_ratios = vkOutput     # original variable

# set peak threshold value
# note that a dynamic threshold function has been made and should be implemented
vkThreshold = getattr(args, "threshold")
if 0 <= vkThreshold <= 100:         # argparse has already removed non-integers
  vkThreshold = vkThreshold * 0.01
  error = vkThreshold       # original variable name
else:
  raise ValueError("The given threshold, %i, is out of bounds." % (vkThreshold))
  exit()

# read polarity
vkPolarity = getattr(args, 'polarity').lower()
if vkPolarity == 'positive' or vkPolarity == '+':
  vKPolarity = 'pos'
elif vkPolarity == 'negative' or vkPolarity == '-':
  vKPolarity = 'neg'
if vkPolarity != 'both' and vkPolarity != 'pos' and vkPolarity != 'neg':
  raise ValueError('Polarity was set to "%s". It must be set to "both", "pos", or "neg".' % (vkPolarity))
  exit()

# read plot-type
vkPlotType = getattr(args, 'plottype')
if vkPlotType == '':       # --input cannot be empty
  raise ValueError("No plot type specified. See -h for help menu.")
  exit()
else:
  vkPlotTypes = []
  for vkPlot in vkPlotType:
    #if vkPlot.lower().with(('scatter', 'heatmap', '3d')):        # verify files extension of --input
    if 'scatter' or 'heatmap' or '3d' in vkPlot.lower():        # verify files extension of --input
      vkPlotTypes.append(vkPlot.lower())
    else:
      raise ValueError('Input was set to "%s". It must be set to the filepath of a mzML or mzXML file.' % (vkFile))
      exit()

# prepare to process mzs from input files
neg_pos_mz_sets = [[], []]

# Use correct processing based on file type. Gets both pos and neg mode because this check is made anyways s o little
# cost to just gather both sets here.
for vkFile in vkInputFiles:
  if vkFile.lower().endswith('.mzxml'):
      mzXML = MzXML()
      mzXML.parse_file(vkFile)
      neg_pos_mz_sets_temp = XML_process(mzXML, threshold=vkThreshold)
      neg_pos_mz_sets[0] = neg_pos_mz_sets[0] + neg_pos_mz_sets_temp[0]
      neg_pos_mz_sets[1] = neg_pos_mz_sets[1] + neg_pos_mz_sets_temp[1]
  elif vkFile.lower().endswith('.mzxml'):
      # Use ML processing
      neg_pos_mz_sets_temp = ML_process(f, threshold=vkThreshold)
      neg_pos_mz_sets[0] = neg_pos_mz_sets[0] + neg_pos_mz_sets_temp[0]
      neg_pos_mz_sets[1] = neg_pos_mz_sets[1] + neg_pos_mz_sets_temp[1]
# Removes all duplicates from both neg and pos lists
neg_pos_mz_sets[0] = list(set(neg_pos_mz_sets[0]))
neg_pos_mz_sets[1] = list(set(neg_pos_mz_sets[1]))

mode = vkPolarity

# get lookup table.
lt = bmrb.getLookupTable('bmrb-db2.csv')

# set up. elements could be changed but would need to do some editing elsewhere.
elements = ['C', 'H', 'O', 'N']
neg_comps = list()
pos_comps = list()

# check to see if user wanted to process neg mode scans
if mode == 'neg' or mode == 'both':
  for mz in neg_pos_mz_sets[0]:
    neg_comps.append(bmrb.getFormulaFromMass(bmrb.adjust(mz, 'neg'), lt, tolerance=error)) # adjust mass and search in lookup table. Store result in list.
  neg_comps = filter(lambda a: a != 'No Match', neg_comps) # Filter out no matches
  neg_elements = extractNeededElementalData.find_elements_values(elements_to_find=elements, compounds=neg_comps) # Get elements from compounds
  neg_ratios = processElementalData.process_elemental_data(neg_elements) # Turn elements into ratios
  if write_ratios == True: 
    neg_filename = 'example-ratios-neg.csv'  # forcing file name. FLAG
    with open(neg_filename, mode='w') as f: 
      for ratio in neg_ratios:
        f.writelines(str(ratio).strip('[]') + '\n')
  plotType = vkPlotTypes
  for vkType in vkPlotTypes:
    print(vkType)
    plotVanK(ratiosList=neg_ratios, typeOfPlot=vkType)

if mode == 'pos' or mode == 'both':
  for mz in neg_pos_mz_sets[1]:
    pos_comps.append(bmrb.getFormulaFromMass(bmrb.adjust(mz, 'pos'), lt)) # adjust mass and search in lookup table. Store result in list.
  pos_comps = filter(lambda a: a != 'No Match', pos_comps) # Filter out no matches
  pos_elements = extractNeededElementalData.find_elements_values(elements_to_find=elements, compounds=pos_comps) # Get elements from compounds
  pos_ratios = processElementalData.process_elemental_data(pos_elements) # Turn elements into ratios
  if write_ratios == True: 
    pos_filename = 'example-ratios-pos.csv'  # forcing file name. FLAG
    with open(neg_filename, mode='w') as f: 
      for ratio in neg_ratios:
        f.writelines(str(ratio).strip('[]') + '\n')
  plotType = vkPlotTypes
  for vkType in vkPlotTypes:
    print(vkType)
    plotVanK(ratiosList=pos_ratios, typeOfPlot=vkType)

