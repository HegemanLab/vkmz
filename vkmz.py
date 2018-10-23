#!/usr/bin/env python

import argparse
import csv
import math
import re

parser = argparse.ArgumentParser()
inputSubparser = parser.add_subparsers(help='Select mode:', dest='mode')
parse_tsv = inputSubparser.add_parser('tsv', help='Use tabular data as input.')
parse_tsv.add_argument('--input', '-i', required=True, help='Path to tabular file. Must include columns: sample ID, mz, polarity, intensity, & retention time.')
parse_xcms = inputSubparser.add_parser('xcms', help='Use XCMS data as input.')
parse_xcms.add_argument('--data-matrix', '-xd', required=True, nargs='?', type=str, help='Path to XCMS data matrix file.')
parse_xcms.add_argument('--sample-metadata', '-xs', required=True, nargs='?', type=str, help='Path to XCMS sample metadata file.')
parse_xcms.add_argument('--variable-metadata', '-xv', required=True, nargs='?', type=str, help='Path to XCMS variable metadata file.')
for inputSubparser in [parse_tsv, parse_xcms]:
  inputSubparser.add_argument('--output',   '-o', nargs='?', type=str, required=True, help='Specify output file path.')
  inputSubparser.add_argument('--json',     '-j', action='store_true', help='Set JSON flag to save JSON output.')
  inputSubparser.add_argument('--error',    '-e', nargs='?', type=float, required=True, help='Mass error of mass spectrometer in parts-per-million.')
  inputSubparser.add_argument('--database', '-db', nargs='?', default='databases/bmrb-light.tsv', help='Select database of known formula masses.')
  inputSubparser.add_argument('--directory','-dir', nargs='?', default='', type=str, help='Define path of tool directory. Assumes relative path if unset.')
  inputSubparser.add_argument('--polarity', '-p', choices=['positive','negative'], help='Force polarity mode to positive or negative. Overrides variables in input file.')
  inputSubparser.add_argument('--neutral',  '-n', action='store_true', help='Set neutral flag if masses in input data are neutral. No mass adjustmnet will be made.')
  inputSubparser.add_argument('--alternate','-a', action='store_true', help='Set flag to keep features with multiple predictions.')
args = parser.parse_args()

class Feature(object):
  '''Feature Class'''
  def __init__(self, sample_id, polarity, mz, rt, intensity, charge = 1):
    self.sample_id = sample_id
    self.polarity = polarity
    self.mz = mz
    self.rt = rt
    self.intensity = intensity
    self.predictions = [] # list of Prediction objects
    self.charge = charge
    # CAUTION: a charge of one is assumed.

class Prediction(object):
  '''Prediction Class'''
  def __init__(self, mass, formula, delta, elementCount, hc, oc, nc):
    self.formula = formula
    self.mass = mass
    self.delta = delta
    self.elementCount = elementCount
    self.hc = hc
    self.oc = oc
    self.nc = nc

# input constants
MODE = getattr(args, "mode")
OUTPUT = getattr(args, "output")
JSON = getattr(args, "json")

def polaritySanitizer(polarity: str):
  """Sanitize input polarity values.

  Renames the, case-insensitive, values 'positive', 'pos', or '+' to 'positive'
  and 'negative', 'neg', or '-' to 'negative'.

  Errors on unrecognized polarity value.
  """
  if polarity.lower() in ['positive','pos','+']:
    polarity = 'positive'
  elif polarity.lower() in ['negative', 'neg', '-']:
    polarity = 'negative'
  else:
    raise ValueError("%s is an unknown polarity type." % polarity)
  return polarity

def readTabular(tsvFile):
  """Read a tabular file and return a list of predictions.

  Columns must be...
  """
  featureList = []
  try:
    with open(tsvFile, 'r') as f:
      next(f) # skip hearder line
      tsvData = csv.reader(f, delimiter='\t')
      for row in tsvData:
        feature = Feature(row[0], polaritySanitizer(row[1]), float(row[2]),\
                          float(row[3]), float(row[4]))
        featureList.append(feature)
  except IOError:
    print('Error while reading %s.' % tsvFile)
  return featureList

def readXcmsTabular(xcmsSampleFile, xcmsVariableFile, xcmsDataMatrixFile):
  """Read W4M's XCMS tabular files and return a list of predictions.

  The sample metadata file is read to assign polarities to samples.

  The variable metadata file is read to 
  """
  featureList = []
  # extract sample polarities
  try:
    polarity = {}
    with open(xcmsSampleFile,'r') as f:
      xcmsSampleMetadata = csv.reader(f, delimiter='\t')
      next(xcmsSampleMetadata) # skip header
      for row in xcmsSampleMetadata:
        sample = row[0]
        if POLARITY:
          polarity[sample] = POLARITY
        else:
          polarity[sample] =  polaritySanitizer(row[2])
  except IOError:
    print('Error while reading the XCMS tabular file %s.' % xcmsSampleFile)
  # extract variable mz & rt
  try:
    mz = {}
    rt = {}
    mzIndex = int()
    rtIndex = int()
    with open(xcmsVariableFile,'r') as f:
      xcmsVariableData = csv.reader(f, delimiter='\t')
      header = next(xcmsVariableData)
      mzIndex = header.index("mz")
      rtIndex = header.index("rt")
      for row in xcmsVariableData:
        mz[row[0]] = float(row[mzIndex])
        rt[row[0]] = float(row[rtIndex])
  except IOError:
    print('Error while reading the XCMS tabular file %s.' % xcmsVariableFile)
  # extract intensity and build Feature objects
  try:
    with open(xcmsDataMatrixFile,'r') as f:
      xcmsDataMatrix = csv.reader(f, delimiter='\t')
      samples = next(xcmsDataMatrix)
      # remove empty columns
      # bug in W4M XCMS input?
      # check w/ fresh input
      samples = [x for x in samples if x is not '']
      for row in xcmsDataMatrix:
        # remove empty columns
        row = [x for x in row if x is not '']
        i = 1
        while(i < len(row)):
          intensity = row[i] # keep type string for test
          if intensity not in {"NA", "#DIV/0!", '0'}:
            variable = row[0]
            sample = samples[i]
            feature = Feature(sample, polarity[sample], mz[variable],\
                              rt[variable], float(intensity))
            featureList.append(feature)
          i+=1
  except IOError:
    print('Error while reading the XCMS tabular file %s.' % xcmsDataMatrixFile)
  return featureList


POLARITY = getattr(args, "polarity")
if MODE == "tsv":
  tsvFile = getattr(args, "input")
  featureList = readTabular(tsvFile)
else: # MODE == "xcms"
  xcmsSampleFile = getattr(args, "sample_metadata")
  xcmsVariableFile = getattr(args, "variable_metadata")
  xcmsDataMatrixFile = getattr(args, "data_matrix")
  featureList = readXcmsTabular(xcmsSampleFile, xcmsVariableFile, xcmsDataMatrixFile)

# store or generate remaining constants
MASS_ERROR = getattr(args, "error")
ALTERNATE = getattr(args, "alternate")
NEUTRAL = getattr(args, "neutral")
DATABASE = getattr(args, "database")
DIRECTORY = getattr(args, "directory")
MASS = []
FORMULA = []
try:
  with open(DIRECTORY+DATABASE, 'r') as tsv:
    next(tsv)
    for row in tsv:
      mass, formula = row.split()
      MASS.append(float(mass))
      FORMULA.append(formula)
except ValueError:
  print('The %s database could not be loaded.' % DATABASE)
MAX_MASS_INDEX = len(MASS)-1

def adjust(mass, polarity, charges = 1):
  """
  Adjust a charged mass to a neutral mass.

  Mass of an electrons (1.007276) is multiplied from the charge and subtracted
  from positively charged ions and added to negatively charged ions.
  """
  # value to adjust by
  proton = 1.007276
  if polarity == 'positive':
    mass -= (proton * charges)
  else: # sanitized to negative
    mass += (proton * charges)
  return mass

def predict(mass, uncertainty, left, right):
  """Search for a matching mass within the known-mass list.

  Uses binary search to match a given mass to a known mass within a given
  uncertainty. Upon match returns known mass index.

  If no match is made returns -1.
  """
  mid = int(((right - left) / 2) + left)
  if left <= mid <= right and mid <= MAX_MASS_INDEX:
    delta = mass - MASS[mid]
    if uncertainty >= abs(delta):
      return mid
    elif uncertainty > delta:
      return predict(mass, uncertainty, left, mid-1)
    else:
      return predict(mass, uncertainty, mid+1, right)
  return -1

def predictAll(mass, uncertainty, initIndex):
  """Search for all matching masses within the known-mass list.

  Checks adjacent indexes from a given index of known-masses which are within a
  given uncertinty of a given mass.
  """
  matches = [initIndex]
  i = 0
  while initIndex+i+1 <= MAX_MASS_INDEX:
    m = initIndex+i+1
    delta = mass - MASS[m]
    if uncertainty >= abs(delta):
      matches.append(m)
      i += 1
    else:
      break
  i = 0
  while initIndex+i-1 >= 0:
    m = initIndex+i-1
    delta = float(MASS[m]) - mass
    if uncertainty >= abs(delta):
      matches.append(m)
      i -= 1
    else:
      break
  return matches

def parseFormula(formula):
  # calculate elemental ratios
  formulaList = re.findall('[A-Z][a-z]?|[0-9]+', formula)
  elementCount = {}
  hc = float()
  oc = float()
  nc = float()
  i = 0;
  while i < len(formulaList):
    # create keys
    if formulaList[i] not in elementCount:
      elementCount[formulaList[i]] = 0
    # if there is only one of this element
    if i+1 == len(formulaList) or formulaList[i+1].isalpha():
      elementCount[formulaList[i]] += 1
    else:
      elementCount[formulaList[i]] += int(formulaList[i+1])
      i+=1
    i+=1
  if 'C' in elementCount:
    if 'H' in elementCount : hc = elementCount['H']/float(elementCount['C'])
    if 'O' in elementCount : oc = elementCount['O']/float(elementCount['C'])
    if 'N' in elementCount : nc = elementCount['N']/float(elementCount['C'])
  return elementCount, hc, oc, nc


def featurePrediction(feature):
  """Make predictions for a feature.

  Reads a Feature as input and, if possible, returns it with a list of
  Prediction objects.

  By default, a feature is assumed to be charged. The observed charged mass of
  the feature is converted to a neutral mass through adjust(). The --neutral
  flag disables adjustment.

  predict() returns an index for the MASS/FORMULA lists  if a known mass is
  within the mass error uncertainty of the observed, neutral, mass.  Features
  without a prediction are thrown out.

  On match, predictAll() searches for matches adjacent to the initial match.

  By default, features with multiple predictions are thrown out unless the
  --alternate flag is set. Alternate matches are sorted by absolute delta.

  For each match an elementCount dictionary is parsed and elemental ratios are
  calculated.

  Prediction objects are made for each match and added to the features
  predictions list before returning the feature object.
  """
  if NEUTRAL:
    mass = feature.mz
  else:
    mass = adjust(feature.mz, feature.polarity, feature.charge)
  # mass error in parts per million
  uncertainty = mass * MASS_ERROR / 1e6
  initIndex = predict(mass, uncertainty, 0, MAX_MASS_INDEX) 
  if initIndex != -1:
    matches = predictAll(mass, uncertainty, initIndex)
    # remove feature if multiple predictions are made and --alternate not set
    if not ALTERNATE and len(matches) > 1:
      return
    for m in matches:
      delta = mass - MASS[m] # check with Stephen
      formula = FORMULA[m]
      elementCount, hc, oc, nc = parseFormula(formula)
      feature.predictions.append(Prediction(MASS[m], formula, delta,\
                                 elementCount, hc, oc, nc))
    # sort alternate matches by lowest absolute delta
    if not ALTERNATE and len(matches) > 1:
      feature.predictions.sort(key=lambda m: abs(m.delta))
    return(feature)
  # no prediction was made
  else:
    return
 
# write output file
def write(predictedFeatures):
  json = ''
  try: 
    with open(OUTPUT+'.tsv', 'w') as tsvFile:
      tabularHeader = "sample_id\t\
                       polarity\t\
                       mz\t\
                       rt\t\
                       intensity\t\
                       predicted_mass\t\
                       predicted_delta\t\
                       predicted_formula\t\
                       predicted_element_count\t\
                       predicted_hc\t\
                       predicted_oc\t\
                       predicted_nc\n"
      if ALTERNATE:
        tabularHeader = tabularHeader[:-1]+"\talternate_predictions\n"
      tsvFile.writelines(tabularHeader)
      for f in predictedFeatures:
        p = f.predictions[0]
        tsvRow = f.sample_id+'\t'+\
                 f.polarity+'\t'+\
                 str(f.mz)+'\t'+\
                 str(f.rt)+'\t'+\
                 str(f.intensity)+'\t'+\
                 str(p.mass)+'\t'+\
                 str(p.delta)+'\t'+\
                 p.formula+'\t'+\
                 str(p.elementCount)+'\t'+\
                 str(p.hc)+'\t'+\
                 str(p.oc)+'\t'+\
                 str(p.nc)+'\n'
        jsonElement = "{\n  \"sample_id\": \""+f.sample_id+\
                "\",\n  \"polarity\": \""+f.polarity+\
                "\",\n  \"mz\": "+str(f.mz)+\
                ",\n  \"rt\": "+str(f.rt)+\
                ",\n  \"intensity\": "+ str(f.intensity)+\
                ",\n  \"prediction\": {\n    \"mass\": "+str(p.mass)+\
                ",\n    \"delta\": \""+str(p.delta)+\
                ",\n    \"formula\": \""+p.formula+\
                "\",\n    \"element_count\": \""+str(p.elementCount)+\
                "\",\n    \"hc\": "+str(p.hc)+\
                ",\n    \"oc\": "+str(p.oc)+\
                ",\n    \"nc\": "+str(p.nc)+"\n  }\n},\n"
        if ALTERNATE and len(f.predictions) > 1:
          tsvAppend = []
          jsonAppend = str()
          for a in f.predictions[1:]:
            tsvAppend.append((a.mass, a.formula, a.delta))
            jsonAppend += "     {\n       \"mass\": "+str(a.mass)+\
                ",\n       \"delta\": \""+str(a.delta)+\
                ",\n       \"formula\": \""+a.formula+\
                "\",\n       \"element_count\": \""+str(a.elementCount)+\
                "\",\n       \"hc\": "+str(a.hc)+\
                ",\n       \"oc\": "+str(a.oc)+\
                ",\n       \"nc\": "+str(a.nc)+"\n     },\n"
          tsvRow = tsvRow[:-1]+'\t'+str(tsvAppend)+'\n'
          jsonElement = jsonElement[:-4]+",\n  \"alternate_predictions\": [\n"\
                        +jsonAppend[:-2]+"\n  ]\n},\n"
        tsvFile.writelines(tsvRow)
        json += jsonElement
    json = json[:-2] # remove final comma # [:-1] ??
  except IOError as error:
    print('IOError while writing tabular output: %s' % error.strerror)
  try: 
    with open(DIRECTORY+'d3.html','r',encoding='utf-8') as htmlTemplate,\
         open(OUTPUT+'.html','w',encoding='utf-8') as htmlFile:
      for line in htmlTemplate:
        line = re.sub('^var data.*$','var data = ['+json+']',line,flags=re.M)
        htmlFile.write(line)
  except IOError as error:
    print('IOError while writing HTML output or reading HTML template: %s'\
          % error.strerror)
  if JSON:
    try: 
      with open(OUTPUT+'.json','w') as jsonFile:
        jsonFile.write(json)
    except IOError as error:
      print('IOError while writing JSON output: %s' % error.strerror)

# main
predictedFeatures = map(featurePrediction, featureList)
predictedFeatures = [x for x in predictedFeatures if x is not None]
# sort by intensity so that D3 draws largest symbols first
predictedFeatures.sort(key=lambda x: x.intensity, reverse=True)
write(predictedFeatures)
