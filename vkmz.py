'''
based on the BMRB compound database which can be found at:
http://www.bmrb.wisc.edu/ftp/pub/bmrb/relational_tables/metabolomics/Chem_comp.csv
'''

import re
import argparse
import multiprocessing
from multiprocessing import Pool
import csv

import numpy as np
import math
import pandas as pd
from plotly import __version__
import plotly.offline as py
import plotly.graph_objs as go
 
parser = argparse.ArgumentParser()
inputSubparser = parser.add_subparsers(help='Select input type:', dest='input-type')
parse_tsv = inputSubparser.add_parser('tsv', help='Use tabular data as input.')
parse_tsv.add_argument('--input', '-i', required=True, help='Path to tabular file. Must include columns: sample ID, mz, polarity, intensity, & retention time.')
parse_tsv.add_argument('--no-plot', '-n', action='store_true', help='Disable plot generation.')
parse_xcms = inputSubparser.add_parser('xcms', help='Use xcms data as input.')
parse_xcms.add_argument('--data-matrix', '-xd', required=True, nargs='?', type=str, help='Path to xcms dataMatrix file.')
parse_xcms.add_argument('--sample-metadata', '-xs', required=True, nargs='?', type=str, help='Path to xcms sampleMetadata file.')
parse_xcms.add_argument('--variable-metadata', '-xv', required=True, nargs='?', type=str, help='Path to xcms variableMetadata file.')
parse_xcms.add_argument('--no-plot', '-n', action='store_true', help='Disable plot generation.')
parse_plot = inputSubparser.add_parser('plot', help='Only plot data.')
parse_plot.add_argument('--input', '-i', required=True, nargs='?', type=str, help='Path to vkmz generated tabular file.')
for inputSubparser in [parse_tsv, parse_xcms]:
  inputSubparser.add_argument('--output',   '-o', nargs='?', type=str, required=True, help='Specify output file path.')
  inputSubparser.add_argument('--error',    '-e', nargs='?', type=int, default=5,     help='PPM error for identification.')
  inputSubparser.add_argument('--database', '-d', nargs='?', default='databases/bmrb-light.tsv', help='Select database.')
  inputSubparser.add_argument('--multiprocessing', '-m', action='store_true', help='Use flag to turn on multiprocessing.')
  inputSubparser.add_argument('--plottype', '-p', nargs='?', default='scatter-2d', choices=['scatter-2d', '2d', 'scatter-3d', '3d', 'heatmap'], help='Set to "2d" or "3d". Default is "2d".')
  inputSubparser.add_argument('--size',     '-s', nargs='?', default=5, type=int, help='Set size of of dots. size+2*log(size*peak/(highest_peak/lowest_peak')
  inputSubparser.add_argument('--sizealgo', '-a', nargs='?', default=0, type=int, choices=[0,1,2],help='Size algorithm selector. Algo 0: size, Algo 1: size+2*log(size*peak/(highest_peak/lowest_peak, Algo 2: size+2*size*peak/(highest_peak-lowest_peak)')
args = parser.parse_args()

vkInputType = getattr(args, "input-type")

# read inputs, arguments and define globals
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

vkPlotType = getattr(args, 'plottype')
if vkPlotType == '2d': vkPlotType = 'scatter-2d'
if vkPlotType == '3d': vkPlotType = 'scatter-3d'

vkSize = getattr(args, 'size')

vkSizeAlgo = getattr(args, 'sizealgo')

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
  mass = adjust(feature[2], feature[1]) # mz & polarity
  prediction = predict(mass, 0, vkMaxIndex)
  if prediction != -1:
    predictions = predictNeighbors(mass, prediction)
    feature[5] = predictions
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

# Binary search to match observed mass to known mass within error
# https://en.wikipedia.org/wiki/Binary_search_tree
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
      f.writelines(str("sample id\tpolarity\tmz\tretention time\tintensity\tpredictions\tdelta\tH:C\tO:C\tN:C") + '\n')
      for feature in vkOutputList:
        predictedDelta = feature[5][0][2]
        predictedFormula = feature[5][0][3]
        hc = float(predictedFormula['H'])/float(predictedFormula['C'])
        oc = float(predictedFormula['H'])/float(predictedFormula['C'])
        nc = float(predictedFormula['H'])/float(predictedFormula['C'])
        f.writelines(feature[0]+'\t'+feature[1]+'\t'+str(feature[2])+'\t'+str(feature[3])+'\t'+str(feature[4])+'\t'+str(feature[5])+'\t'+str(predictedDelta)+'\t'+str(hc)+'\t'+str(oc)+'\t'+str(nc)+'\n')
  except ValueError:
    print('"%s" could not be saved.' % filename)

def plotRatios(identified, type):
  traces = []
  trace_count = 0
  lowest_peak = 10.0**10
  highest_peak = 0.0
  feature_rts =[]
  highest_rt = 0
  feature_formulas = []
  feature_size = []
  x=[]
  y=[]
  # 3d data is given to all plots for added plot.ly functionality
  z=[]
  for feature in identified:
    feature_peak = feature[2]
    if feature_peak > highest_peak:
      highest_peak = feature_peak
    elif feature_peak < lowest_peak:
      lowest_peak = feature_peak
    if highest_rt < feature[3]:
      higher_rt = feature[3]
  # assign metadata to each feature
  for feature in identified:
    feature_peak = feature[2]
    # changes retention time from seconds to minutes
    feature_rts.append(feature[3]/60)
    # feature_size algorithm is not complete
    if vkSizeAlgo == 0:
      feature_size.append(vkSize)
    elif vkSizeAlgo == 1:
      feature_size.append(vkSize+4*vkSize*feature_peak/(highest_peak-lowest_peak))
    else: 
      feature_size.append(vkSize+2*math.log(vkSize*feature_peak/(highest_peak-lowest_peak)))
    feature_formulas.append(feature[4])
    x.append(feature[6]) # Oxygen / Carbon
    y.append(feature[5]) # Hydrogen / Carbon
    z.append(feature[7]) # Nitrogen / Carbon
  if type == 'scatter-3d':
    feature_trace = go.Scatter3d(
      x = x,
      y = y,
      z = z,
      mode='markers',
      text=feature_formulas,
      marker=dict(
        size=feature_size,
        color=feature_rts,
        colorscale='Viridis',
        colorbar=dict(title='Retention Time (m)'),
        line=dict(width=0.5),
        opacity=0.8
      )
    )
    traces.append(feature_trace)
    layout = go.Layout(
      title="Van Krevelen Diagram", 
      scene = dict(
        xaxis= dict(
          title= 'Carbon to Oxygen Ratio',
          zeroline= False,
          gridcolor='rgb(183,183,183)',
          showline=True
        ),
        yaxis= dict(
          title= 'Hydrogen to Carbon Ratio',
          zeroline= False,
          gridcolor='rgb(183,183,183)',
          showline=True
        ),
        zaxis=dict(
          title= 'Carbon to Nitrogen Ratio',
          zeroline= False,
          gridcolor='rgb(183,183,183)',
          showline=True
        ),
      ), 
      margin=dict(r=0, b=0, l=0, t=100)
    )
    fig = go.Figure(data=traces, layout=layout)
    py.plot(fig, filename=vkOutput, auto_open=False)
  if type == 'scatter-2d':
    feature_trace = go.Scatter(
      x = x,
      y = y,
      mode='markers',
      text=feature_formulas,
      marker=dict(
        size=feature_size,
        color=feature_rts,
        colorscale='Viridis',
        colorbar=dict(title='Retention Time (m)'),
        line=dict(width=0.5),
        opacity=0.8
      )
    )
    traces.append(feature_trace)
    layout = go.Layout(
      title="Van Krevelen Diagram", 
      xaxis= dict(
        title= 'Carbon to Oxygen Ratio',
        zeroline= False,
        gridcolor='rgb(183,183,183)',
        showline=True
      ),
      yaxis= dict(
        title= 'Hydrogen to Carbon Ratio',
        zeroline= False,
        gridcolor='rgb(183,183,183)',
        showline=True
      ),
      margin=dict(r=0, b=100, l=100, t=100)
    )
    fig = go.Figure(data=traces, layout=layout)
    py.plot(fig, filename=vkOutput, auto_open=False)


# main
if vkInputType == "tsv":
  vkInput = []
  tsvFile = getattr(args, "input")
  try:
    with open(tsvFile, 'r') as f:
      next(f) # skip hearder line
      tsvData = csv.reader(f, delimiter='\t')
      for row in tsvData:
        vkInput.append([row[0],row[1],float(row[2]),float(row[3]),float(row[4]),[]])
  except ValueError:
    print('The %s data file could not be read.' % tsvFile)
  vkmzData = forecaster(vkInput)
  saveForcast(vkmzData)
  plot = getattr(args, "no_plot")
  if plot:
    plotRatios(vkmzData, vkPlotType)
elif vkInputType == "xcms":
  vkInput = []
  xcmsSampleMetadataFile = getattr(args, "sample_metadata")
  try:
    polarity = {}
    with open(xcmsSampleMetadataFile, 'r') as f:
      xcmsSampleMetadata = csv.reader(f, delimiter='\t')
      for row in xcmsSampleMetadata:
        polarity[row[0]] = row[2]
  except ValueError:
    print('The %s data file could not be read.' % xcmsSampleMetadataFile)
  xcmsVariableMetadataFile = getattr(args, "variable_metadata")
  try:
    mz = {}
    rt = {}
    with open(xcmsVariableMetadataFile, 'r') as f:
      xcmsVariableMetadata = csv.reader(f, delimiter='\t')
      for row in xcmsVariableMetadata:
        mz[row[0]] = row[2]
        rt[row[0]] = row[3]
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
              variable = row[0]
              if variable != "NA":
                intensity = row[i]
                sample = sample_id[i]
                vkInput.append([sample, polarity[sample], float(mz[variable]), float(rt[variable]), intensity, []])
            i+=1
  except ValueError:
    print('The %s data file could not be read.' % xcmsDataMatrixFile)
  vkmzData = forecaster(vkInput)
  saveForcast(vkmzData)
  plot = getattr(args, "no_plot")
  if plot:
    plotRatios(vkmzData, vkPlotType)
else:
  vkmzData = []
  tsvPlotvFile = getattr(args, "input")
  try:
    with open(tsvPlotFile, 'r') as f:
      next(f) # skip header line
      plotData = csv.reader(f, delimiter='\t')
      for row in plotData:
        vkmzData.append([row[0],row[1],float(row[2]),float(row[3]),float(row[4]),list(row[4]),float(row[5]),float(row[6]),float(row[7]),float(row[8])])
  except ValueError:
    print('The %s data file could not be read.' % tsvFile)
  plotRatios(vkmzData, vkPlotType)
