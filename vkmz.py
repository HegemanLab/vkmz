#!/usr/bin/env python

import argparse
import csv
import math
import re
from decimal import Decimal

parser = argparse.ArgumentParser()
inputSubparser = parser.add_subparsers(help='Select mode:', dest='mode')
parse_db = inputSubparser.add_parser('db', help='Construct Database.')
parse_db.add_argument('--input', '-i', required=True, help='Path to tabular database file. The first two columns must be "mass" and "formula".')
parse_db.add_argument('--output','-o', nargs='?', type=str, required=True, help='Specify output file path.')
parse_db.add_argument('--mass','-m', nargs='?', type=str, help='Header column name for masses.')
parse_db.add_argument('--formula','-f', nargs='?', type=str, required=True, help='Header column name for formulas.')
parse_db.add_argument('--calculate','-c', action='store_true', help='Set calculate flag to calculate monoisotopic masses from molecular formulas.')
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

featureList = []
class Feature(object):
 def __init__(self, sample_id, polarity, mz, rt, intensity):
   self.sample_id = sample_id
   self.polarity = polarity
   self.mz = mz
   self.rt = rt
   self.intensity = intensity
   self.predictions = [] # list of Prediction objects

class Prediction(object):
  def __init__(self, mass, formula, delta):
    self.formula = formula
    self.elementCount = {}
    self.mass = mass
    self.delta = delta
    self.hc = float()
    self.nc = float()
    self.oc = float()

# store input constants
MODE = getattr(args, "mode")
OUTPUT = getattr(args, "output")
JSON = getattr(args, "json")

if MODE == "db":
  dbFileIn = getattr(args, "input")
  formula_column = getattr(args, "formula")
  mass_column = getattr(args, "mass")
  calculate = getattr(args, "calculate")
  #print("formula_column", bool(formula_column))
  #print("mass_column", bool(mass_column))
  #print("calculate", bool(calculate))
  if mass_column and calculate:
    raise ValueError("The mass and precssion parameters cannot both be set.")
  try:
    with open(dbFileIn, 'r') as f:
      database = []
      dbData = csv.reader(f, delimiter='\t')
      header = next(dbData) 
      formula_index = 0 if formula_column == None else header.index(formula_column)
      # test that mass_column works
      if mass_column:
        mass_index = header.index(mass_column)
        index_length = formula_index if formula_index >= mass_index else mass_index # test that this works
        for row in dbData:
          # sanity check that mass-formula indices exist
          if len(row) > index_length:
            formula = row[formula_index]
            mass = row[mass_index]
            # create tuple for each mass-formula pair only when both exist
            if mass and not mass.isspace() and formula and not formula.isspace():
              database.append((float(mass), formula.strip().replace(" ","")))
      # calculate is default
      else:
        # Dictionary of IUPAC atomic weight measurements
        # see https://github.com/HegemanLab/atomicWeightsDecimal
        elementDictionary = {"H":Decimal((0,(1,0,0,7,9,4),-5)),"He":Decimal((0,(4,0,0,2,6,0,2),-6)),"Li":Decimal((0,(6,9,4,1),-3)),"Be":Decimal((0,(9,0,1,2,1,8,2),-6)),"B":Decimal((0,(1,0,8,1,1),-3)),"C":Decimal((0,(1,2,0,1,0,7),-4)),"N":Decimal((0,(1,4,0,0,6,7),-4)),"O":Decimal((0,(1,5,9,9,9,4),-4)),"F":Decimal((0,(1,8,9,9,8,4,0,3,2),-7)),"Ne":Decimal((0,(2,0,1,7,9,7),-5)),"Na":Decimal((0,(2,2,9,8,9,7,7,0),-6)),"Mg":Decimal((0,(2,4,3,0,5,0),-4)),"Al":Decimal((0,(2,6,9,8,1,5,3,8),-6)),"Si":Decimal((0,(2,8,0,8,5,5),-4)),"P":Decimal((0,(3,0,9,7,3,7,6,1),-6)),"S":Decimal((0,(3,2,0,6,5),-3)),"Cl":Decimal((0,(3,5,4,5,3),-3)),"Ar":Decimal((0,(3,9,9,4,8),-3)),"K":Decimal((0,(3,9,0,9,8,3),-4)),"Ca":Decimal((0,(4,0,0,7,8),-3)),"Sc":Decimal((0,(4,4,9,5,5,9,1,0),-6)),"Ti":Decimal((0,(4,7,8,6,7),-3)),"V":Decimal((0,(5,0,9,4,1,5),-4)),"Cr":Decimal((0,(5,1,9,9,6,1),-4)),"Mn":Decimal((0,(5,4,9,3,8,0,4,9),-6)),"Fe":Decimal((0,(5,5,8,4,5),-3)),"Co":Decimal((0,(5,8,9,3,3,2,0,0),-6)),"Ni":Decimal((0,(5,8,6,9,3,4),-4)),"Cu":Decimal((0,(6,3,5,4,6),-3)),"Zn":Decimal((0,(6,5,4,0,9),-3)),"Ga":Decimal((0,(6,9,7,2,3),-3)),"Ge":Decimal((0,(7,2,6,4),-2)),"As":Decimal((0,(7,4,9,2,1,6,0),-5)),"Se":Decimal((0,(7,8,9,6),-2)),"Br":Decimal((0,(7,9,9,0,4),-3)),"Kr":Decimal((0,(8,3,7,9,8),-3)),"Rb":Decimal((0,(8,5,4,6,7,8),-4)),"Sr":Decimal((0,(8,7,6,2),-2)),"Y":Decimal((0,(8,8,9,0,5,8,5),-5)),"Zr":Decimal((0,(9,1,2,2,4),-3)),"Nb":Decimal((0,(9,2,9,0,6,3,8),-5)),"Mo":Decimal((0,(9,5,9,4),-2)),"Ru":Decimal((0,(1,0,1,0,7),-2)),"Rh":Decimal((0,(1,0,2,9,0,5,5,0),-5)),"Pd":Decimal((0,(1,0,6,4,2),-2)),"Ag":Decimal((0,(1,0,7,8,6,8,2),-4)),"Cd":Decimal((0,(1,1,2,4,1,1),-3)),"In":Decimal((0,(1,1,4,8,1,8),-3)),"Sn":Decimal((0,(1,1,8,7,1,0),-3)),"Sb":Decimal((0,(1,2,1,7,6,0),-3)),"Te":Decimal((0,(1,2,7,6,0),-2)),"I":Decimal((0,(1,2,6,9,0,4,4,7),-5)),"Xe":Decimal((0,(1,3,1,2,9,3),-3)),"Cs":Decimal((0,(1,3,2,9,0,5,4,5),-5)),"Ba":Decimal((0,(1,3,7,3,2,7),-3)),"La":Decimal((0,(1,3,8,9,0,5,5),-4)),"Ce":Decimal((0,(1,4,0,1,1,6),-3)),"Pr":Decimal((0,(1,4,0,9,0,7,6,5),-5)),"Nd":Decimal((0,(1,4,4,2,4),-2)),"Sm":Decimal((0,(1,5,0,3,6),-2)),"Eu":Decimal((0,(1,5,1,9,6,4),-3)),"Gd":Decimal((0,(1,5,7,2,5),-2)),"Tb":Decimal((0,(1,5,8,9,2,5,3,4),-5)),"Dy":Decimal((0,(1,6,2,5,0,0),-3)),"Ho":Decimal((0,(1,6,4,9,3,0,3,2),-5)),"Er":Decimal((0,(1,6,7,2,5,9),-3)),"Tm":Decimal((0,(1,6,8,9,3,4,2,1),-5)),"Yb":Decimal((0,(1,7,3,0,4),-2)),"Lu":Decimal((0,(1,7,4,9,6,7),-3)),"Hf":Decimal((0,(1,7,8,4,9),-2)),"Ta":Decimal((0,(1,8,0,9,4,7,9),-4)),"W":Decimal((0,(1,8,3,8,4),-2)),"Re":Decimal((0,(1,8,6,2,0,7),-3)),"Os":Decimal((0,(1,9,0,2,3),-2)),"Ir":Decimal((0,(1,9,2,2,1,7),-3)),"Pt":Decimal((0,(1,9,5,0,7,8),-3)),"Au":Decimal((0,(1,9,6,9,6,6,5,5),-5)),"Hg":Decimal((0,(2,0,0,5,9),-2)),"Tl":Decimal((0,(2,0,4,3,8,3,3),-4)),"Pb":Decimal((0,(2,0,7,2),-1)),"Bi":Decimal((0,(2,0,8,9,8,0,3,8),-5)),"Th":Decimal((0,(2,3,2,0,3,8,1),-4)),"Pa":Decimal((0,(2,3,1,0,3,5,8,8),-5)),"U":Decimal((0,(2,3,8,0,2,8,9,1),-5))}
        for row in dbData:
          # sanity check that formula index exist
          if len(row) > formula_index:
            formula = row[formula_index]
            # remove null and empty formulas
            if formula and not formula.isspace():
              # create list of atomic elements and counts
              # e.g., "CH4" becomes ["C","H","4"]
              formulaList = re.findall('[A-Z][a-z]?|[0-9]+', formula)
              # max precission in IUPAC's measurements
              significantExponent = -7
              formulaMass = Decimal((0,(0,0),significantExponent))
              i = 0;
              while i < len(formulaList):
                if formulaList[i] in elementDictionary:
                  elementMass = elementDictionary[formulaList[i]]
                  elementExponent = elementMass.as_tuple().exponent
                  if elementExponent > significantExponent : significantExponent = elementExponent
                  # elementCount == 1
                  if i+1 == len(formulaList) or formulaList[i+1].isalpha():
                    formulaMass += elementMass
                  else:
                    elementCount = int(formulaList[i+1])
                    formulaMass += elementMass * elementCount
                    i+=1
                # fail softly
                else:
                  print("The formula '%s' has an unknown element: '%s'", (formula, formulaList[i]))
                  formulaMass = "error"
                  i = len(formulaList)
                i+=1
              if formulaMass != "error":
                # round non-significant figures
                formulaMass = formulaMass.quantize(Decimal((0,(1,0),significantExponent)))
                # store tuple
                database.append((formulaMass, formula.strip().replace(" ","")))
      # remove duplicates and sort
      database = sorted(set(database))
      try:
        with open(OUTPUT, 'w') as dbFileOut: 
          dbFileOut.write("mass\tformula\n")
          for pair in database:
            dbFileOut.write(str(pair[0])+'\t'+pair[1]+'\n')
      except ValueError:
        print('Error while writing the %s database file.' % dbFileOut)
  except ValueError:
    print('Error while reading the %s tabular database file.' % dbFileIn)
  exit()

def polaritySanitizer(sample_polarity):
  if sample_polarity.lower() in {'positive','pos','+'}:
    sample_polarity = 'positive'
  elif sample_polarity.lower() in {'negative', 'neg', '-'}:
    sample_polarity = 'negative'
  else:
    print('A sample has an unknown polarity type: %s. Polarity in the XCMS sample metadata should be set to "negative" or "positive".' % sample_polarity)
    raise ValueError
  return sample_polarity

# read input modes
POLARITY = getattr(args, "polarity")
if MODE == "tsv":
  tsvFile = getattr(args, "input")
  try:
    with open(tsvFile, 'r') as f:
      next(f) # skip hearder line
      tsvData = csv.reader(f, delimiter='\t')
      for row in tsvData:
        feature = Feature(row[0], polaritySanitizer(row[1]), float(row[2]), float(row[3]), float(row[4]))
        featureList.append(feature)
  except ValueError:
    print('The %s data file could not be read.' % tsvFile)
else: # MODE == "xcms"
  # extract sample polarities
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
  # extract variable mz & rt
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
  # extract intensity and build Feature objects
  xcmsDataMatrixFile = getattr(args, "data_matrix")
  try:
    with open(xcmsDataMatrixFile, 'r') as f:
      xcmsDataMatrix = csv.reader(f, delimiter='\t')
      samples = next(xcmsDataMatrix, None)
      # remove empty columns, bug?
      samples = [x for x in samples if x is not '']
      for row in xcmsDataMatrix:
        row = [x for x in row if x is not '']
        i = 1
        while(i < len(row)):
          intensity = row[i] # keep as string for test
          if intensity not in {"NA", "#DIV/0!", '0'}:
            variable = row[0]
            sample_id = samples[i]
            feature = Feature(sample_id, polarity[sample], mz[variable], rt[variable], float(intensity))
            featureList.append(feature)
          i+=1
  except ValueError:
    print('The %s data file could not be read.' % xcmsDataMatrixFile)

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
      MASS.append(mass)
      FORMULA.append(formula)
except ValueError:
  print('The %s database could not be loaded.' % DATABASE)
MAX_MASS_INDEX = len(MASS)-1

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
  mid = int(((right - left) / 2) + left)
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
  neighbors = [Prediction(float(MASS[prediction]),FORMULA[prediction],(float(MASS[prediction])-mass))]
  while prediction+i+1 <= MAX_MASS_INDEX:
    neighbor = prediction+i+1
    delta = float(MASS[neighbor])-mass
    if uncertainty >= abs(delta):
      neighbors.append(Prediction(float(MASS[neighbor]),FORMULA[neighbor],delta))
      i += 1
    else:
      break
  i = 0
  while prediction + i - 1 >= 0:
    neighbor = prediction + i - 1
    delta = float(MASS[neighbor]) - mass
    if uncertainty >= abs(delta):
      neighbors.append(Prediction(float(MASS[neighbor]),FORMULA[neighbor],delta))
      i -= 1
    else:
      break
  # sort predictions by lowest absolute delta
  neighbors.sort(key=lambda x: abs(x.delta))
  return neighbors

# predict formulas by the mass of a feature
def featurePrediction(feature):
  if NEUTRAL:
    mass = feature.mz
  else:
    mass = adjust(feature.mz, feature.polarity) # mz & polarity
  uncertainty = mass * MASS_ERROR / 1e6
  # initialPrediction is an index
  initialPrediction = predict(mass, uncertainty, 0, MAX_MASS_INDEX) 
  if initialPrediction != -1: # else feature if forgotten
    feature.predictions = predictNeighbors(mass, uncertainty, initialPrediction)
    if not ALTERNATE and len(feature.predictions) > 1:
      return
    # calculate elemental ratios for the lowest abs(delta) prediction
    for p in feature.predictions:
      formulaList = re.findall('[A-Z][a-z]?|[0-9]+', p.formula)
      i = 0;
      while i < len(formulaList):
        # create key
        if formulaList[i] not in p.elementCount:
          p.elementCount[formulaList[i]] = 0
        # if there is only one of this element
        if i+1 == len(formulaList) or formulaList[i+1].isalpha():
          p.elementCount[formulaList[i]] += 1
        else:
          p.elementCount[formulaList[i]] += int(formulaList[i+1])
          i+=1
        i+=1
      if 'C' in p.elementCount:
        if 'H' in p.elementCount : p.hc = float(p.elementCount['H'])/float(p.elementCount['C'])
        if 'O' in p.elementCount : p.oc = float(p.elementCount['O'])/float(p.elementCount['C'])
        if 'N' in p.elementCount : p.nc = float(p.elementCount['N'])/float(p.elementCount['C'])
    return(feature)
 
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
        if ALTERNATE:
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
          jsonElement = jsonElement[:-4]+",\n  \"alternate_predictions\": [\n"+jsonAppend[:-2]+"\n  ]\n},\n"
        tsvFile.writelines(tsvRow)
        json += jsonElement
    json = json[:-2] # remove final comma # [:-1] ??
  except IOError as error:
    print('IOError while writing tabular output: %s' % error.strerror)
  try: 
    with open(DIRECTORY+'d3.html','r',encoding='utf-8') as htmlTemplate,\
         open(OUTPUT+'.html','w',encoding='utf-8') as htmlFile:
      for line in htmlTemplate:
        line = re.sub('^var data.*$', 'var data = ['+json+']', line, flags=re.M)
        htmlFile.write(line)
  except IOError as error:
    print('IOError while writing HTML output or reading HTML template: %s' % error.strerror)
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
