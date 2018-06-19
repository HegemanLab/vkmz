import re
import argparse
import multiprocessing
from multiprocessing import Pool
import csv
import math
import pandas as pd
from plotly import __version__
import plotly.offline as py
import plotly.graph_objs as go
 
parser = argparse.ArgumentParser()
inputSubparser = parser.add_subparsers(help='Select input type:', dest='input-type')
parse_tsv = inputSubparser.add_parser('tsv', help='Use tabular data as input.')
parse_tsv.add_argument('--input', '-i', required=True, help='Path to tabular file. Must include columns: sample ID, mz, polarity, intensity, & retention time.')
parse_tsv.add_argument('--no-plot', '-np', action='store_true', help='Disable plot generation.')
parse_xcms = inputSubparser.add_parser('xcms', help='Use XCMS data as input.')
parse_xcms.add_argument('--data-matrix', '-xd', required=True, nargs='?', type=str, help='Path to XCMS dataMatrix file.')
parse_xcms.add_argument('--sample-metadata', '-xs', required=True, nargs='?', type=str, help='Path to XCMS sampleMetadata file.')
parse_xcms.add_argument('--variable-metadata', '-xv', required=True, nargs='?', type=str, help='Path to XCMS variableMetadata file.')
parse_xcms.add_argument('--no-plot', '-np', action='store_true', help='Disable plot generation.')
# Dprecated by future D3 plot
#parse_plot = inputSubparser.add_parser('plot', help='Only plot data.')
#parse_plot.add_argument('--input', '-i', required=True, nargs='?', type=str, help='Path to VKMZ generated tabular file.')
for inputSubparser in [parse_tsv, parse_xcms]:
  inputSubparser.add_argument('--output',   '-o', nargs='?', type=str, required=True, help='Specify output file path.')
  inputSubparser.add_argument('--error',    '-e', nargs='?', type=float, required=True, help='Mass error of mass spectrometer in parts-per-million.')
  inputSubparser.add_argument('--database', '-d', nargs='?', default='databases/bmrb-light.tsv', help='Select database of known formula masses.')
  inputSubparser.add_argument('--directory', nargs='?', default='', type=str, help='Define path of tool directory. Assumes relative path if unset.')
  inputSubparser.add_argument('--polarity', '-p', choices=['positive','negative'], help='Force polarity mode to positive or negative. Overrides variables in input file.')
  inputSubparser.add_argument('--no-adjustment', '-na', action='store_true', help='Set flag to diasble charged mass adjustment.')
  inputSubparser.add_argument('--unique', '-u', action='store_true', help='Set flag to remove features with multiple predictions.')
  inputSubparser.add_argument('--multiprocessing', '-m', action='store_true', help='Set flag to enable multiprocessing.')
# currently non-functional
#  inputSubparser.add_argument('--plottype', '-t', nargs='?', default='scatter-2d', choices=['scatter-2d', 'scatter-3d'], help='Select plot type.')
  inputSubparser.add_argument('--size',     '-s', nargs='?', default=5, type=int, help='Set maxium size of plot symbols.')
  inputSubparser.add_argument('--size-algorithm', '-sa', choices=['uniform', 'relative-intensity', 'relative-log-intensity'], default='uniform', help='Set size algorithm selector.')
args = parser.parse_args()

vkInputType = getattr(args, "input-type")
vkPolarity = getattr(args, "polarity")
vkError = getattr(args, "error")


vkUnique = getattr(args, "unique")

vkMultiprocessing = getattr(args, "multiprocessing")

vkNoAdjustment = getattr(args, "no_adjustment")

vkDatabaseFile = getattr(args, "database")
vkDirectory = getattr(args, "directory")

vkMass = []
vkFormula = []
try:
  with open(vkDirectory+vkDatabaseFile, 'r') as tsv:
    next(tsv) # skip first row
    for row in tsv:
      mass, formula = row.split()
      vkMass.append(mass)
      vkFormula.append(formula)
except ValueError:
  print('The %s database could not be loaded.' % vkDatabaseFile)
vkMaxIndex = len(vkMass)-1

vkOutput = getattr(args, "output")

# Dprecated by future D3 plot
#vkPlotType = getattr(args, 'plottype')

vkSize = getattr(args, 'size')

vkSizeAlgorithm = getattr(args, 'size_algorithm')

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
  if vkNoAdjustment:
    mass = feature[2]
  else:
    mass = adjust(feature[2], feature[1]) # mz & polarity
  uncertainty = mass * vkError / 1e6
  prediction = predict(mass, uncertainty, 0, vkMaxIndex)
  if prediction != -1:
    predictions = predictNeighbors(mass, uncertainty, prediction)
    if vkUnique and len(predictions) > 1:
      return
    feature.append(predictions) # feature[5]
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
    hc = float(formulaDictionary['H'])/float(formulaDictionary['C'])
    oc = float(formulaDictionary['O'])/float(formulaDictionary['C'])
    nc = float(formulaDictionary['N'])/float(formulaDictionary['C'])
    feature += [hc, oc, nc]
    return(feature)

# adjust charged mass to a neutral mass
def adjust(mass, polarity):
  # value to adjust by
  proton = 1.007276
  if polarity == 'positive':
    mass -= proton
  else: # sanitized to negative
    mass += proton
  return mass

# Binary search to match observed mass to known mass within error
# https://en.wikipedia.org/wiki/Binary_search_tree
def predict(mass, uncertainty, left, right):
  mid = ((right - left) / 2) + left
  if left <= mid <= right and mid <= vkMaxIndex:
    delta = float(vkMass[mid]) - mass
    if uncertainty >= abs(delta):
      return mid
    elif uncertainty < delta:
      return predict(mass, uncertainty, left, mid-1)
    else:
      return predict(mass, uncertainty, mid+1, right)
  return -1
  
def polaritySanitizer(sample_polarity):
  if sample_polarity.lower() in {'positive','pos','+'}:
    sample_polarity = 'positive'
  elif sample_polarity.lower() in {'negative', 'neg', '-'}:
    sample_polarity = 'negative'
  else:
    print('A sample has an unknown polarity type: %s. Polarity in the XCMS sample metadata should be set to "negative" or "positive".' % sample_polarity)
    raise ValueError
  return sample_polarity

def plotData(vkData):
  max_rt = 0.0
  max_intensity = 0.0
  for row in vkData:
    intensity = row[4]
    if intensity > max_intensity:
      max_intensity = intensity
    rt = row[3]
    if rt > max_rt:
      max_rt = rt
  if vkSizeAlgorithm == 'uniform':
    for row in vkData:
      row.append(vkSize)
  elif vkSizeAlgorithm == 'relative-intensity':
    alpha = vkSize/max_intensity
    for row in vkData:
      intensity = row[4]
      row.append(alpha*intensity)
  else: # relative-log-intensity
    alpha = vkSize/math.log(max_intensity+1)
    for row in vkData:
      intensity = row[4]
      row.append(alpha*math.log(intensity+1))
  for row in vkData:
    rt = row[3]
    row.append(rt/max_rt)
  return vkData

# find and sort known masses within error limit of observed mass
def predictNeighbors(mass, uncertainty, prediction):
  i = 0
  neighbors = [[vkMass[prediction],vkFormula[prediction],(float(vkMass[prediction])-mass)],]
  while prediction+i+1 <= vkMaxIndex:
    neighbor = prediction+i+1
    delta = float(vkMass[neighbor])-mass
    if uncertainty >= abs(delta):
      neighbors.append([vkMass[neighbor],vkFormula[neighbor],delta])
      i += 1
    else:
      break
  i = 0
  while prediction+i-1 >= 0:
    neighbor = prediction+i-1
    delta = float(vkMass[neighbor])-mass
    if uncertainty >= abs(delta):
      neighbors.append([vkMass[neighbor],vkFormula[neighbor],(float(vkMass[neighbor])-mass)])
      i -= 1
    else:
      break
  neighbors = sorted(neighbors, key = (lambda delta: abs(delta[2])))
  return neighbors

# write output file
def saveForcast(vkOutputList):
  try: 
    with open(vkOutput+'.tsv', 'w') as f: 
      f.writelines(str("sample_id\tpolarity\tmz\trt\tintensity\tpredictions\thc\toc\tnc\tsize\tcolor") + '\n')
      for feature in vkOutputList:
        f.writelines(feature[0]+'\t'+feature[1]+'\t'+str(feature[2])+'\t'+str(feature[3])+'\t'+str(feature[4])+'\t'+str(feature[5])+'\t'+str(feature[6])+'\t'+str(feature[7])+'\t'+str(feature[8])+'\t'+str(feature[9])+'\t'+str(feature[10])+'\n')
  except ValueError:
    print('"%s" could not be saved.' % filename)

def plotRatios(vkData):
  max_rt = 0
  max_hc = 0
  max_oc = 0
  max_nc = 0
  for row in vkData:
    if row[3] > max_rt:
      max_rt = row[3]
    if row[7] > max_hc:
      max_hc = row[7]
    if row[8] > max_oc:
      max_oc = row[8]
    if row[9] > max_nc:
      max_nc = row[9]
  labels = ['sampleID', 'polarity', 'mz', 'rt', 'intensity', 'predictions', 'hc', 'oc', 'nc', 'symbol_size', 'color']
  df = pd.DataFrame.from_records(vkData, columns=labels)
  sampleIDs = df.sampleID.unique()
  data = []
  menus = []
  i = 0
  for sampleID in sampleIDs:
    dfSample = df.loc[df['sampleID'] == sampleID]
    size = dfSample.symbol_size
    trace = go.Scatter(
      x = dfSample.oc,
      y = dfSample.hc,
      text = dfSample.predictions.apply(lambda x: "Prediction: "+str(x[0][1])+"<br>mz: " +str(x[0][0])+"<br>Delta: "+str(x[0][2])),
      line = dict(width = 0.5),
      mode = 'markers',
      marker = dict(
        size = size,
        sizemode = "area",
        color = dfSample.rt,
        colorscale = 'Viridis',
        cmin = 0,
        cmax = max_rt,
        colorbar=dict(title='Retention Time (s)'),
        line = dict(width = 0.5),
        showscale = True
      ),
      opacity = 0.8
    )
    data.append(trace)
    vision = []
    j = 0
    while j < len(sampleIDs):
      if j != i:
        vision.append(False)
      else:
        vision.append(True)
      j += 1
    menu = dict(
      method = 'update',
      label = sampleID,
      args = [{'visible': vision}, {'title': sampleID}]
    )
    menus.append(menu)
    i += 1
  updatemenus = list([
    dict(
      active = -1,
      buttons = menus
    )
  ])
  layout = go.Layout(
    title = "Van Krevelen Diagram",
    showlegend = False,
    xaxis = dict(
      title = 'Oxygen to Carbon Ratio',
      zeroline = False,
      gridcolor = 'rgb(183,183,183)',
      showline = True,
      range = [0, max_oc]
    ),
    yaxis = dict(
      title = 'Hydrogen to Carbon Ratio',
      zeroline = False,
      gridcolor = 'rgb(183,183,183)',
      showline = True,
      range = [0, max_hc]
    ),
    margin = dict(r=0, b=100, l=100, t=100),
    updatemenus = updatemenus
  )     
  fig = go.Figure(data=data, layout=layout)
  py.plot(fig, auto_open=False, show_link=False, filename=vkOutput+'.html')
 
# main
if vkInputType == "tsv":
  vkInput = []
  tsvFile = getattr(args, "input")
  try:
    with open(tsvFile, 'r') as f:
      next(f) # skip hearder line
      tsvData = csv.reader(f, delimiter='\t')
      for row in tsvData:
        vkInput.append([row[0],polaritySanitizer(row[1]),float(row[2]),float(row[3]),float(row[4])])
  except ValueError:
    print('The %s data file could not be read.' % tsvFile)
  vkData = forecaster(vkInput)
  vkData = plotData(vkData)
  saveForcast(vkData)
  plotRatios(vkData)
# Dprecated by future D3 plot
#elif vkInputType == "xcms":
else: # vkInputType == "xcms"
  vkInput = []
  xcmsSampleMetadataFile = getattr(args, "sample_metadata")
  try:
    polarity = {}
    with open(xcmsSampleMetadataFile, 'r') as f:
      xcmsSampleMetadata = csv.reader(f, delimiter='\t')
      next(xcmsSampleMetadata, None) # skip header
      for row in xcmsSampleMetadata:
        sample = row[0]
        if vkPolarity:
          polarity[sample] = vkPolarity
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
  vkData = forecaster(vkInput)
  vkData = plotData(vkData)
  saveForcast(vkData)
  plotRatios(vkData)
# Dprecated by future D3 plot
#else:
#  vkData = []
#  tsvPlotvFile = getattr(args, "input")
#  try:
#    with open(tsvPlotFile, 'r') as f:
#      next(f) # skip header line
#      plotData = csv.reader(f, delimiter='\t')
#      for row in plotData:
#        vkData.append([row[0],row[1],float(row[2]),float(row[3]),float(row[4]),list(row[4]),float(row[5]),float(row[6]),float(row[7]),float(row[8])])
#  except ValueError:
#    print('The %s data file could not be read.' % tsvFile)
#  plotRatios(vkData)

