import argparse
import csv
import math
import re
 
parser = argparse.ArgumentParser()
inputSubparser = parser.add_subparsers(help='Select input type:', dest='input-type')
parse_tsv = inputSubparser.add_parser('tsv', help='Use tabular data as input.')
parse_tsv.add_argument('--input', '-i', required=True, help='Path to tabular file. Must include columns: sample ID, mz, polarity, intensity, & retention time.')
parse_xcms = inputSubparser.add_parser('xcms', help='Use XCMS data as input.')
parse_xcms.add_argument('--data-matrix', '-xd', required=True, nargs='?', type=str, help='Path to XCMS dataMatrix file.')
parse_xcms.add_argument('--sample-metadata', '-xs', required=True, nargs='?', type=str, help='Path to XCMS sampleMetadata file.')
parse_xcms.add_argument('--variable-metadata', '-xv', required=True, nargs='?', type=str, help='Path to XCMS variableMetadata file.')
for inputSubparser in [parse_tsv, parse_xcms]:
  inputSubparser.add_argument('--output',   '-o', nargs='?', type=str, required=True, help='Specify output file path.')
  inputSubparser.add_argument('--error',    '-e', nargs='?', type=float, required=True, help='Mass error of mass spectrometer in parts-per-million.')
  inputSubparser.add_argument('--database', '-db', nargs='?', default='databases/bmrb-light.tsv', help='Select database of known formula masses.')
  inputSubparser.add_argument('--directory','-dir', nargs='?', default='', type=str, help='Define path of tool directory. Assumes relative path if unset.')
  inputSubparser.add_argument('--polarity', '-p', choices=['positive','negative'], help='Force polarity mode to positive or negative. Overrides variables in input file.')
  inputSubparser.add_argument('--neutral',  '-n', action='store_true', help='Set neutral flag if masses in input data are neutral. No mass adjustmnet will be made.')
  inputSubparser.add_argument('--unique',   '-u', action='store_true', help='Set flag to remove features with multiple predictions.')
  inputSubparser.add_argument('--size',     '-s', nargs='?', default=5, type=int, help='Set maxium size of plot symbols.')
args = parser.parse_args()

# store input constants
INPUT_TYPE = getattr(args, "input-type")
POLARITY = getattr(args, "polarity")

def polaritySanitizer(sample_polarity):
  if sample_polarity.lower() in {'positive','pos','+'}:
    sample_polarity = 'positive'
  elif sample_polarity.lower() in {'negative', 'neg', '-'}:
    sample_polarity = 'negative'
  else:
    print('A sample has an unknown polarity type: %s. Polarity in the XCMS sample metadata should be set to "negative" or "positive".' % sample_polarity)
    raise ValueError
  return sample_polarity

# read input
vkInput = [] # each element is a feature from input
if INPUT_TYPE == "tsv":
  tsvFile = getattr(args, "input")
  try:
    with open(tsvFile, 'r') as f:
      next(f) # skip hearder line
      tsvData = csv.reader(f, delimiter='\t')
      for row in tsvData:
        vkInput.append([row[0],polaritySanitizer(row[1]),float(row[2]),float(row[3]),float(row[4])])
  except ValueError:
    print('The %s data file could not be read.' % tsvFile)
else: # INPUT_TYPE == "xcms"
  xcmsSampleMetadataFile = getattr(args, "sample_metadata")
  try:
    polarity = {}
    with open(xcmsSampleMetadataFile, 'r') as f:
      xcmsSampleMetadata = csv.reader(f, delimiter='\t')
      next(xcmsSampleMetadata, None) # skip header
      for row in xcmsSampleMetadata:
        sample = row[0]
        if POLARITY:
          polarity[sample] = POLARITY
        else:
          sample_polarity = polaritySanitizer(row[2])
          polarity[sample] = sample_polarity
  except ValueError:
    print('The %s data file could not be read. Check that polarity is set to "negative" or "positive"' % xcmsSampleMetadataFile)
  xcmsVariableMetadataFile = getattr(args, "variable_metadata")
  try:
    mz = {}
    rt = {}
    variable_index = {}
    mz_index = int()
    rt_index = int()
    with open(xcmsVariableMetadataFile, 'r') as f:
      xcmsVariableMetadata = csv.reader(f, delimiter='\t')
      i = 0
      for row in xcmsVariableMetadata:
        if i != 0:
          mz[row[0]] = float(row[mz_index])
          rt[row[0]] = float(row[rt_index])
        else:
          for column in row:
            variable_index[column] = i
            i += 1
          mz_index = variable_index["mz"]
          rt_index = variable_index["rt"]
  except ValueError:
    print('The %s data file could not be read.' % xcmsVariableMetadataFile)
  xcmsDataMatrixFile = getattr(args, "data_matrix")
  try:
    with open(xcmsDataMatrixFile, 'r') as f:
      xcmsDataMatrix = csv.reader(f, delimiter='\t')
      first_row = True
      for row in xcmsDataMatrix:
        if first_row:
          sample_id = row
          first_row = False
        else: 
          i = 0
          while(i < len(row)):
            if i == 0:
              i+=1
            else:
              intensity = row[i]
              if intensity not in {'NA', '#DIV/0!', '0'}:
                variable = row[0]
                sample = sample_id[i]
                # XCMS data may include empty columns
                if sample != "":
                  vkInput.append([sample, polarity[sample], mz[variable], rt[variable], float(intensity)])
            i+=1
  except ValueError:
    print('The %s data file could not be read.' % xcmsDataMatrixFile)

# generate remaining constants
OUTPUT = getattr(args, "output")
MASS_ERROR = getattr(args, "error")
UNIQUE = getattr(args, "unique")
NEUTRAL = getattr(args, "neutral")
DATABASE = getattr(args, "database")
DIRECTORY = getattr(args, "directory")
MASS = []
FORMULA = []
try:
  with open(DIRECTORY+DATABASE, 'r') as tsv:
    next(tsv) # skip first row
    for row in tsv:
      mass, formula = row.split()
      MASS.append(mass)
      FORMULA.append(formula)
except ValueError:
  print('The %s database could not be loaded.' % DATABASE)
MAX_MASS_INDEX = len(MASS)-1
MAX_SIZE = getattr(args, 'size')
SIZE_ALGORITHM = getattr(args, 'size_algorithm')

# adjust charged mass to a neutral mass
def adjust(mass, polarity):
  # value to adjust by
  proton = 1.007276
  if polarity == 'positive':
    mass -= proton
  else: # sanitized to negative
    mass += proton
  return mass

# binary search to match a neutral mass to known mass within error
def predict(mass, uncertainty, left, right):
  mid = ((right - left) / 2) + left
  if left <= mid <= right and mid <= MAX_MASS_INDEX:
    delta = float(MASS[mid]) - mass
    if uncertainty >= abs(delta):
      return mid
    elif uncertainty < delta:
      return predict(mass, uncertainty, left, mid-1)
    else:
      return predict(mass, uncertainty, mid+1, right)
  return -1

# find and rank predictions which are adjacent to the index of an intial prediction
def predictNeighbors(mass, uncertainty, prediction):
  i = 0
  neighbors = [[MASS[prediction],FORMULA[prediction],(float(MASS[prediction])-mass)],]
  while prediction+i+1 <= MAX_MASS_INDEX:
    neighbor = prediction+i+1
    delta = float(MASS[neighbor])-mass
    if uncertainty >= abs(delta):
      neighbors.append([MASS[neighbor],FORMULA[neighbor],delta])
      i += 1
    else:
      break
  i = 0
  while prediction+i-1 >= 0:
    neighbor = prediction+i-1
    delta = float(MASS[neighbor])-mass
    if uncertainty >= abs(delta):
      neighbors.append([MASS[neighbor],FORMULA[neighbor],(float(MASS[neighbor])-mass)])
      i -= 1
    else:
      break
  neighbors = sorted(neighbors, key = (lambda delta: abs(delta[2])))
  return neighbors

# predict formulas by the mass of a feature
def featurePrediction(feature):
  if NEUTRAL:
    mass = feature[2]
  else:
    mass = adjust(feature[2], feature[1]) # mz & polarity
  uncertainty = mass * MASS_ERROR / 1e6
  prediction = predict(mass, uncertainty, 0, MAX_MASS_INDEX)
  if prediction != -1: # else feature if forgotten
    predictions = predictNeighbors(mass, uncertainty, prediction)
    if UNIQUE and len(predictions) > 1:
      return
    feature.append(predictions) # feature[5]
    formula = predictions[0][1] # formula of prediction with lowest abs(delta)
    formulaList = re.findall('[A-Z][a-z]?|[0-9]+', formula)
    formulaDictionary = {'C':0,'H':0,'O':0,'N':0} # other elements are easy to add
    i = 0;
    while i < len(formulaList):
      if formulaList[i] in formulaDictionary:
        # if there is only one of this element
        if i+1 == len(formulaList) or formulaList[i+1].isalpha():
          formulaDictionary[formulaList[i]] = 1
        else:
          formulaDictionary[formulaList[i]] = formulaList[i+1]
          i+=1
      i+=1
    hc = float(formulaDictionary['H'])/float(formulaDictionary['C'])
    oc = float(formulaDictionary['O'])/float(formulaDictionary['C'])
    nc = float(formulaDictionary['N'])/float(formulaDictionary['C'])
    feature += [hc, oc, nc] # feature[6], [7], [8]
    return(feature)
 
# write output file
def write(vkData):
  json = ''
  try: 
    # write tabular file and generate json for html output
    with open(OUTPUT+'.tsv', 'w') as f: 
      f.writelines(str("sample_id\tpolarity\tmz\trt\tintensity\tpredictions\thc\toc\tnc") + '\n')
      for feature in vkData:
        f.writelines(feature[0]+'\t'+feature[1]+'\t'+str(feature[2])+'\t'+str(feature[3])+'\t'+str(feature[4])+'\t'+str(feature[5])+'\t'+str(feature[6])+'\t'+str(feature[7])+'\t'+str(feature[8])+'\n')
        json += '{sample_id:\''+str(feature[0])+'\', polarity:\''+str(feature[1])+'\', mz:'+str(feature[2])+', rt:'+str(feature[3])+', intensity:'+str(feature[4])+', predictions:'+str(feature[5])+', hc:'+str(feature[6])+', oc:'+str(feature[7])+', nc:'+str(feature[8])+'},'
    json = json[:-1] # remove final comma
    # write html
    try:
      with open(DIRECTORY+'d3.html', 'r') as template, open(OUTPUT+'.html', 'w') as f:
       for line in template:
         line = re.sub('^var data.*$', 'var data = ['+json+']', line, flags=re.M)
         f.write(line)
    except ValueError:
      print('"%s" could not be read or "%s" could not be written' % template, f)
  except ValueError:
    print('"%s" could not be saved.' % filename)

# main
vkData = map(featurePrediction, vkInput)
vkData = [x for x in vkData if x is not None]
write(vkData)
