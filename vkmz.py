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
  inputSubparser.add_argument('--database', '-db', nargs='?', default='databases/bmrb-light.tsv', help='Define database of known fomrula-mass pairs.')
  inputSubparser.add_argument('--directory','-dir', nargs='?', default='', type=str, help='Define tool directory path. Defaults to relative path.')
  inputSubparser.add_argument('--polarity', '-p', choices=['positive','negative'], help='Force polarity mode to positive or negative. Overrides variables in input file.')
  inputSubparser.add_argument('--neutral',  '-n', action='store_true', help='Set neutral flag if masses in input data are neutral. No mass adjustmnet will be made.')
  inputSubparser.add_argument('--alternate','-a', action='store_true', help='Set flag to keep features with multiple predictions.')
  inputSubparser.add_argument('--charge','-c', action='store_true', help='Set flag if input data contains charge.')
args = parser.parse_args()

class Feature(object):
  '''Feature Class'''
  def __init__(self, sample_id, polarity, mz, rt, intensity, variable_id = False,
               charge = 1):
    self.sample_id = sample_id
    self.polarity = polarity
    self.mz = mz
    self.rt = rt
    self.intensity = intensity
    self.predictions = [] # list of Prediction objects
    if variable_id:
      self.variable_id = variable_id
    else:
      self.variable_id = sample_id+'-'+str(rt)+'-'+str(mz)
    self.charge = charge
    # CAUTION: a charge of one is assumed if not specified with --charge

class Prediction(object):
  '''Prediction Class'''
  def __init__(self, mass, formula, delta, element_count, hc, oc, nc):
    self.formula = formula
    self.mass = mass
    self.delta = delta
    self.element_count = element_count
    self.hc = hc
    self.oc = oc
    self.nc = nc

# parameter constants
MODE = getattr(args, "mode")
OUTPUT = getattr(args, "output")
JSON = getattr(args, "json")

# using charge not yet implemented
CHARGE = getattr(args, "charge")

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

def readTabular(tsv_file):
  """Read a tabular file and return a list of features.

  Reads columns named "sample_id", "polarity", "mz", "rt", "intensity", and
  optionally "charge" to create Features for each row.
  """
  feature_list = []
  try:
    with open(tsv_file, 'r') as f:
      tsvData = csv.reader(f, delimiter='\t')
      header = next(tsvData)
      sample_id_index = header.index("sample_id")
      polarity_index = header.index("polarity")
      mz_index = header.index("mz")
      rt_index = header.index("rt")
      intensity_index = header.index("intensity")
      if CHARGE:
        charge_index = header.index("charge")
      for row in tsvData:
        feature = Feature(
                    sample_id = row[sample_id_index],\
                    polarity = polaritySanitizer(row[polarity_index]),\
                    mz = float(row[mz_index]),\
                    rt = float(row[rt_index]),\
                    intensity = float(row[intensity_index]))
        if CHARGE:
          feature.charge = row[charge_index]
        feature_list.append(feature)
  except IOError:
    print('Error while reading %s.' % tsv_file)
  return feature_list

def readXcmsTabular(sample_file, variable_file, matrix_file):
  """Read W4M's XCMS tabular files and return a list of features.

  Reads sample metadata to create a dictionary of sample ids keys with,
  sanitized, polarity values.

  Reads variable metadata to create two dictionaries with variable ids as keys
  and either mass to charge or retention time as values.

  Finally, read data matrix and create all Feature objects and append to list.
  """
  feature_list = []
  # extract sample polarities
  try:
    polarity = {}
    with open(sample_file,'r') as f:
      sample_data = csv.reader(f, delimiter='\t')
      next(sample_data) # skip header
      for row in sample_data:
        sample = row[0]
        if POLARITY:
          polarity[sample] = POLARITY
        else:
          polarity[sample] = polaritySanitizer(row[2])
  except IOError:
    print('Error while reading the XCMS tabular file %s.' % sample_file)
  # extract variable mz & rt
  try:
    mz = {}
    rt = {}
    mzIndex = int()
    rtIndex = int()
    with open(variable_file,'r') as f:
      variable_data = csv.reader(f, delimiter='\t')
      header = next(variable_data)
      mzIndex = header.index("mz")
      rtIndex = header.index("rt")
      for row in variable_data:
        mz[row[0]] = float(row[mzIndex])
        rt[row[0]] = float(row[rtIndex])
  except IOError:
    print('Error while reading the XCMS tabular file %s.' % variable_file)
  # extract intensity and build Feature objects
  try:
    with open(matrix_file,'r') as f:
      matrix_data = csv.reader(f, delimiter='\t')
      samples = next(matrix_data)
      # remove empty columns
      # check w/ fresh input
      samples = [x for x in samples if x is not '']
      for row in matrix_data:
        # remove empty columns
        row = [x for x in row if x is not '']
        i = 1
        while(i < len(row)):
          intensity = row[i] # keep type string for test
          if intensity not in {"NA", "#DIV/0!", '0'}:
            variable_id = row[0]
            sample = samples[i]
            feature = Feature(
                        sample_id = sample,\
                        polarity = polarity[sample],\
                        mz = mz[variable_id],\
                        rt = rt[variable_id],\
                        intensity = float(intensity),\
                        variable_id = variable_id)
            feature_list.append(feature)
          i+=1
  except IOError:
    print('Error while reading the XCMS tabular file %s.' % matrix_file)
  return feature_list


POLARITY = getattr(args, "polarity")
if MODE == "tsv":
  tsv_file = getattr(args, "input")
  feature_list = readTabular(tsv_file)
else: # MODE == "xcms"
  sample_file = getattr(args, "sample_metadata")
  variable_file = getattr(args, "variable_metadata")
  matrix_file = getattr(args, "data_matrix")
  feature_list = readXcmsTabular(sample_file, variable_file, matrix_file)

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

def adjust(mass, polarity, charge):
  """
  Adjust a charged mass to a neutral mass.

  Mass of an electrons (1.007276) is multiplied from the charge and subtracted
  from positively charged ions and added to negatively charged ions.
  """
  # value to adjust by
  proton = 1.007276
  if polarity == 'positive':
    mass -= (proton * charge)
  else: # sanitized to negative
    mass += (proton * charge)
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

def predictAll(mass, uncertainty, init_index):
  """Search for all matching masses within the known-mass list.

  Checks adjacent indexes from a given index of known-masses which are within a
  given uncertinty of a given mass.
  """
  matches = [init_index]
  i = 0
  while init_index+i+1 <= MAX_MASS_INDEX:
    m = init_index+i+1
    delta = mass - MASS[m]
    if uncertainty >= abs(delta):
      matches.append(m)
      i += 1
    else:
      break
  i = 0
  while init_index+i-1 >= 0:
    m = init_index+i-1
    delta = float(MASS[m]) - mass
    if uncertainty >= abs(delta):
      matches.append(m)
      i -= 1
    else:
      break
  return matches

def parseFormula(formula):
  # calculate elemental ratios
  formula_list = re.findall('[A-Z][a-z]?|[0-9]+', formula)
  element_count = {}
  hc = float()
  oc = float()
  nc = float()
  i = 0;
  while i < len(formula_list):
    # create keys
    if formula_list[i] not in element_count:
      element_count[formula_list[i]] = 0
    # if there is only one of this element
    if i+1 == len(formula_list) or formula_list[i+1].isalpha():
      element_count[formula_list[i]] += 1
    else:
      element_count[formula_list[i]] += int(formula_list[i+1])
      i+=1
    i+=1
  if 'C' in element_count:
    if 'H' in element_count : hc = element_count['H']/float(element_count['C'])
    if 'O' in element_count : oc = element_count['O']/float(element_count['C'])
    if 'N' in element_count : nc = element_count['N']/float(element_count['C'])
  return element_count, hc, oc, nc


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

  For each match an element_count dictionary is parsed and elemental ratios are
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
  init_index = predict(mass, uncertainty, 0, MAX_MASS_INDEX) 
  if init_index != -1:
    matches = predictAll(mass, uncertainty, init_index)
    # remove feature if multiple predictions are made and --alternate not set
    if not ALTERNATE and len(matches) > 1:
      return
    for m in matches:
      delta = mass - MASS[m] # check with Stephen
      formula = FORMULA[m]
      element_count, hc, oc, nc = parseFormula(formula)
      feature.predictions.append(Prediction(MASS[m], formula, delta,\
                                 element_count, hc, oc, nc))
    # sort alternate matches by lowest absolute delta
    if not ALTERNATE and len(matches) > 1:
      feature.predictions.sort(key=lambda m: abs(m.delta))
    return(feature)
  # no prediction was made
  else:
    return
 
# write output file
def write(predicted_list):
  json = ''
  try: 
    with open(OUTPUT+'.tsv', 'w') as tsv_file:
      tabularHeader = "sample_id\t\
                       variable_id\t\
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
      tsv_file.writelines(tabularHeader)
      for f in predicted_list:
        p = f.predictions[0]
        tsv_row = f.sample_id+'\t'+\
                 f.variable_id+'\t'+\
                 f.polarity+'\t'+\
                 str(f.mz)+'\t'+\
                 str(f.rt)+'\t'+\
                 str(f.intensity)+'\t'+\
                 str(p.mass)+'\t'+\
                 str(p.delta)+'\t'+\
                 p.formula+'\t'+\
                 str(p.element_count)+'\t'+\
                 str(p.hc)+'\t'+\
                 str(p.oc)+'\t'+\
                 str(p.nc)+'\n'
        json_element = "{\n  \"sample_id\": \""+f.sample_id+\
                "\",\n  \"variable_id\": \""+f.variable_id+\
                "\",\n  \"polarity\": \""+f.polarity+\
                "\",\n  \"mz\": "+str(f.mz)+\
                ",\n  \"rt\": "+str(f.rt)+\
                ",\n  \"intensity\": "+ str(f.intensity)+\
                ",\n  \"prediction\": {\n    \"mass\": "+str(p.mass)+\
                ",\n    \"delta\": \""+str(p.delta)+\
                ",\n    \"formula\": \""+p.formula+\
                "\",\n    \"element_count\": \""+str(p.element_count)+\
                "\",\n    \"hc\": "+str(p.hc)+\
                ",\n    \"oc\": "+str(p.oc)+\
                ",\n    \"nc\": "+str(p.nc)+"\n  }\n},\n"
        if ALTERNATE and len(f.predictions) > 1:
          tsv_append = []
          json_append = str()
          for a in f.predictions[1:]:
            tsv_append.append((a.mass, a.formula, a.delta))
            json_append += "     {\n       \"mass\": "+str(a.mass)+\
                ",\n       \"delta\": \""+str(a.delta)+\
                ",\n       \"formula\": \""+a.formula+\
                "\",\n       \"element_count\": \""+str(a.element_count)+\
                "\",\n       \"hc\": "+str(a.hc)+\
                ",\n       \"oc\": "+str(a.oc)+\
                ",\n       \"nc\": "+str(a.nc)+"\n     },\n"
          tsv_row = tsv_row[:-1]+'\t'+str(tsv_append)+'\n'
          json_element = json_element[:-4]+",\n  \"alternate_predictions\": [\n"\
                        +json_append[:-2]+"\n  ]\n},\n"
        tsv_file.writelines(tsv_row)
        json += json_element
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
predicted_list = map(featurePrediction, feature_list)
predicted_list = [x for x in predicted_list if x is not None]
# sort by intensity so that D3 draws largest symbols first
predicted_list.sort(key=lambda x: x.intensity, reverse=True)
write(predicted_list)
