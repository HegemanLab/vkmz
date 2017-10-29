import numpy as np
import re
import os
import time
import extractNeededElementalData
import processElementalData
import bmrbLookup as bmrb
import argparse
import multiprocessing
from multiprocessing import Pool
from flexPlot import plotVanK
from functools import partial
import csv

parser = argparse.ArgumentParser()
#parser.add_argument('--load',      '-l', nargs='?', default='',     help='Load a previously generated ratio table. Set file path. Disabled by default.')
parser.add_argument('--database',  '-d', nargs='*', default='bmrb.csv',   help='Select database(s).')
parser.add_argument('--input',     '-i', nargs='?', type=str, required=True, help='Enter data source. Data file must be a comma seperated list containing four elements: mz,polarity,intensity,rt.')
parser.add_argument('--output',    '-o', action='store_true',       help='Call variable to ouput a ratio file.')
parser.add_argument('--polarity',  '-p', nargs='?', default='both', type=str, choices=['both', 'pos', 'neg'], help='Set to "pos", "neg", or "both". Default is "both". Only one plot type can be set.')
parser.add_argument('--plottype', '-pt', nargs='*', default=['scatter'], choices=['scatter', 'heatmap', '3d'], help='Set to "scatter", "heatmap", or "3d". Default is "scatter".')
parser.add_argument('--multiprocessing', '-m', action='store_true', help='Call variable to use multiprocessing. One process per core.')
args = parser.parse_args()

# read load argument
#vkLoad = getattr(args, "load")

# read load argument
vkDatabase = getattr(args, "database")
#lt = []
#print(type(lt))
#for database in vkDatabase:
#  lt = bmrb.getLookupTable('databases/' + database)
lt = bmrb.getLookupTable('databases/' + vkDatabase)

# read input argument(s)
vkInput = getattr(args, "input")
vkInputMzs = []
try:
  with open(vkInput, 'rb') as csvfile:
    foo = csv.reader(csvfile)
    for row in foo:
      # each row of csvfile becomes a tuple element inside a list
      # placeholder elements are added to each tuple
      vkInputMzs.append([float(row[0]),row[1],float(row[2]),float(row[3]),[],float])
except ValueError:
  print('The %s data file could not be loaded.' % vkInput)

# read output argument
vkOutput = getattr(args, "output")

# read polarity argument
vkPolarity = getattr(args, 'polarity')

# read plottype argument
vkPlotTypes = getattr(args, 'plottype')

# read multiprocessing argument
vkMultiprocessing = getattr(args, "multiprocessing")

# vk ratio dataset builder for each polarity
# dataset is list with four elements.
# input mass to charge ratios checked against database of masses with known chemical structure
# each identified mass is saved across the four elements, at a given index, of the dataset list
# the first  element represents ratio of hydrogens in identified masses structure
# the second element represents ratio of carbons                    #FLAG VERIRFY
# the third  element represents if a nitrogen is in the strcutre    #FLAG remove element, have plotter check last element != 0
# the fourth element represents ratio of nitrogens
def buildRatios(vkInputMzs):
  # get lookup table.
  # set up. elements could be changed but would need to do some editing elsewhere.
  elements = ['C', 'H', 'O', 'N']
  error = 5
  if vkMultiprocessing:
    try:
      pool = Pool()
      multiprocessMzsArgs = partial(multiprocessMzs, polarity, error)
      identified = pool.map(multiprocessMzsArgs, vkInputMzs[index])
    except Exception as e:
      print(str(e))
    finally:
      pool.close()
      pool.join()
  else:
    identified = []
    for centroid in vkInputMzs:
      identity = bmrb.getFormulaFromMass(bmrb.adjust(centroid[0], centroid[1]), lt, tolerance=error)
      if identity != 'No Match':
        identity_dict = {'C':0,'H':0,'O':0,'N':0} # necessary keys for 3d
        for element in identity:
          regex = re.compile(element+'([0-9]*)')
          value = regex.findall(identity)
          if element.isalpha():
            if value == ['']:
                value = ['1']
            result = int(''.join(value))
            identity_dict[element] = result
        centroid[4] = [identity_dict]
        identified.append(centroid)
        # this would be a good place to add unsaturation (2+2(carbons)+2(nitrogens)-hydrogens)/2
  print(identified)
  if vkOutput: 
    saveRatios(identifiedRatios, polarity)
  for type in vkPlotTypes:
    plotRatios(identified, type)

def multiprocessMzs(polarity, error, inputMz): # recieves a single Mz
  return bmrb.getFormulaFromMass(bmrb.adjust(inputMz, str(polarity)), lt, tolerance=error)

# write vk ratios as csv file
def saveRatios(ratios, polarity):
  try:
    filename = 'ratios-' + time.strftime("%Y%m%d%H%M%S-") + str(polarity) + '.csv'
    with open(filename, 'w') as f: 
      for ratio in ratios:
        f.writelines(str(ratio).strip('[]') + '\n')
  except ValueError:
    print('"%s" could not be saved.' % filename)

## load VK ratio csv file
## FLAG: add multiple file support
#def loadRatios(vkLoad):
#  try:
#    # read in ratios
#    with open(vkLoad, 'r') as f:
#      ratios = f.readlines()
#    # split ratios into proper value
#    for i in range(0, len(ratios)):
#      ratios[i] = ratios[i].split(', ')
#    # Cast to correct type
#    ratios[0] = map(lambda x: float(x), ratios[0])
#    ratios[1] = map(lambda x: float(x), ratios[1])
#    ratios[3] = map(lambda x: float(x), ratios[3])
#    for type in vkPlotTypes:
#      plotRatios(ratios, type)
#  except ValueError:
#    print('The %s data file could not be loaded.' % vkLoad)

def plotRatios(identified, type):
  import pandas as pd
  import plotly as py
  import plotly.graph_objs as go
  #if type == 'scatter':
  #  x_axis = []
  #  for known_id in identified:
  #      x_axis.append(known_id[5][0]
  #  trace1 = Scatter(x=ratios[1], y=ratios[0], mode = 'markers')
  #  layout = Layout(title="<b>Van Krevelin Diagram</b>", 
  #       xaxis= dict(
  #         title= 'Oxygen to Carbon Ratio',
  #         zeroline= False,
  #         gridcolor='rgb(183,183,183)',
  #         showline=True
  #       ),
  #       yaxis=dict(
  #         title= 'Hydrogen to Carbon Ratio',
  #         gridcolor='rgb(183,183,183)',
  #         zeroline=False,
  #         showline=True
  #       ))
  #  plotly.offline.plot({"data": [trace1], "layout": layout}, filename='vk-scatter.html', image='jpeg') 
  if type == '3d':
    from plotly import __version__
    #from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
    import plotly.offline as py
    import plotly.graph_objs as go
    import numpy as np
    traces = []
    trace_count = 0
    for feature in identified:
      feature_formula = feature[4][0]
      if feature_formula['O'] == 0:
        x = 0
      else:
        x = float(feature_formula["C"]/feature_formula["O"])
      if feature_formula['N'] == 0:
        y = 0
      else:
        y = float(feature_formula["C"]/feature_formula["N"])
      if feature_formula['N'] == 0:
        z = 0
      else:
        z = float(feature_formula["C"]/feature_formula["H"])
      x = [x]
      y = [y]
      z = [z]
      feature_trace = go.Scatter3d(
        x = x,
        y = y,
        z = z,
        mode='markers',
        marker=dict(
          size=12,
          line=dict(
          color='rgba(217, 217, 217, 0.14)',
            width=0.5
          ),
          opacity=0.8
        )
      )
      traces.append(feature_trace)
    print(traces)
    layout = go.Layout(
      margin=dict(
        l=0,
        r=0,
        b=0,
        t=0
      )
    )
    fig = go.Figure(data=traces, layout=layout)
    py.plot(fig, filename='simple-3d-scatter.html')

buildRatios(vkInputMzs)

