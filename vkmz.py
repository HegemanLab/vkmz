'''
based on the BMRB compound database which can be found at:
http://www.bmrb.wisc.edu/ftp/pub/bmrb/relational_tables/metabolomics/Chem_comp.csv
'''

import re
import argparse
import multiprocessing
from multiprocessing import Pool
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--input',     '-i', nargs='?', type=str, required=True,  help='Data file must be a csv with the columns mz, polarity, intensity, & rt.')
parser.add_argument('--output',    '-o', nargs='?', type=str, required=True,  help='Specify name of output file.')
parser.add_argument('--database',  '-d', nargs='?', default='databases/bmrb-light.tsv', help='Select database.')
parser.add_argument('--error',     '-e', nargs='?', type=int, default=5,      help='Error in PPM for identification.')
parser.add_argument('--multiprocessing', '-m', action='store_true', help='Use flag to turn on multiprocessing.')
args = parser.parse_args()

# read arguments and define globals
vkInputFile = getattr(args, "input")
vkInput = []
try:
  with open(vkInputFile, 'r') as tsv:
    next(tsv) # skip first line
    tsvData = csv.reader(tsv, delimiter='\t')
    for row in tsvData:
      vkInput.append([float(row[0]),row[1],float(row[2]),float(row[3]),[]])
except ValueError:
  print('The %s data file could not be loaded.' % vkInput)

vkError = getattr(args, "error")

vkMultiprocessing = getattr(args, "multiprocessing")

vkDatabaseFile = getattr(args, "database")
vkMass = []
vkFormula = []
try:
  with open(vkDatabaseFile, 'r') as tsv:
    next(tsv) # skip first row
    for row in tsv:
      (mass, formula) = row.split()
      vkMass.append(mass)
      vkFormula.append(formula)
except ValueError:
  print('The %s database could not be loaded.' % vkDatabase)
vkMaxIndex = len(vkMass)-1

vkOutput = getattr(args, "output")

# control predictions
def forecaster(vkInput):
  if vkMultiprocessing:
    try:
      pool = Pool()
      vkOutputList = pool.map(featurePrediction, vkInput)
    except Exception as e:
      print("Error during multirpocessing: "+str(e))
    finally:
      pool.close()
      pool.join()
  else:
    vkOutputList = map(featurePrediction, vkInput)
  vkOutputList = [x for x in vkOutputList if x is not None]
  return(vkOutputList)

# predict feature formulas and creates output list
def featurePrediction(feature):
  mass = adjust(feature[0], feature[1])
  prediction = predict(mass, 0, vkMaxIndex)
  if prediction != -1:
    predictions = predictNeighbors(mass, prediction)
    feature[4] = predictions
    predictionClosest = predictions[0]
    formula = predictionClosest[1]
    formulaList = re.findall('[A-Z][a-z]?|[0-9]+', formula)
    formulaDictionary = {'C':0,'H':0,'O':0,'N':0}
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
    predictionClosest.append(formulaDictionary)
    return(feature)

# adjust observed mass for polarity
def adjust(mass, polarity):
  # value to adjust by
  hm = 1.007276
  if polarity == 'pos':
    mass += hm
  elif polarity == 'neg':
    mass -= hm
  return mass

# BST to match observed mass to known masses within error
def predict(mass, left, right):
  mid = left + (right - left) / 2
  if left <= mid <= right and mid <= vkMaxIndex:
    delta = float(vkMass[mid]) - mass
    if vkError >= abs(delta):
      return mid
    elif vkError < delta:
      return predict(mass, left, mid-1)
    else:
      return predict(mass, mid+1, right)
  return -1
  
# find and sort known masses within error limit of observed mass
def predictNeighbors(mass, prediction):
  i = 0
  neighbors = [[vkMass[prediction],vkFormula[prediction],(float(vkMass[prediction])-mass)],]
  while prediction+i+1 <= vkMaxIndex:
    neighbor = prediction+i+1
    delta = float(vkMass[neighbor])-mass
    if vkError >= abs(delta):
      neighbors.append([vkMass[neighbor],vkFormula[neighbor],delta])
      i += 1
    else:
      break
  i = 0
  while prediction+i-1 >= 0:
    neighbor = prediction+i-1
    delta = float(vkMass[neighbor])-mass
    if vkError >= abs(delta):
      neighbors.append([vkMass[neighbor],vkFormula[neighbor],(float(vkMass[neighbor])-mass)])
      i -= 1
    else:
      break
  neighbors = sorted(neighbors, key = (lambda delta: abs(delta[2])))
  return neighbors

# write output file
def saveForcast(vkOutputList):
  try:
    with open(vkOutput, 'w') as f: 
      f.writelines(str("mz\tpolarity\tintensity\tretention time\tpredicted formulas\tH:C\tO:c\tN:C\tdelta") + '\n')
      for feature in vkOutputList:
        predictedFormula = feature[4][0][3]
        predictedDelta = feature[4][0][2]
        f.writelines(str(feature[0])+'\t'+str(feature[1])+'\t'+str(feature[2])+'\t'+str(feature[3])+'\t'+str(feature[4])+'\t'+str(float(predictedFormula['H'])/float(predictedFormula['C']))+'\t'+str(float(predictedFormula['O'])/float(predictedFormula['C']))+'\t'+str(float(predictedFormula['N'])/float(predictedFormula['C']))+'\t'+str(predictedDelta)+'\n')
  except ValueError:
    print('"%s" could not be saved.' % filename)

saveForcast(forecaster(vkInput))
