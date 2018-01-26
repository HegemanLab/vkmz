import argparse
import csv
import numpy as np
import math
import pandas as pd
from plotly import __version__
import plotly.offline as py
import plotly.graph_objs as go
 
parser = argparse.ArgumentParser()
parser.add_argument('--input',    '-i', nargs='?', required=True, help='Load a tabular file with ratio information.')
parser.add_argument('--output',   '-o', nargs='?', required=True, help='Load a tabular file with ratio information.')
parser.add_argument('--plottype', '-p', nargs='?', default='scatter-2d', choices=['scatter-2d', '2d', 'scatter-3d', '3d', 'heatmap'], help='Set to "2d" or "3d". Default is "2d".')
parser.add_argument('--size',     '-s', nargs='?', default=5, type=int, help='Set size of of dots. size+2*log(size*peak/(highest_peak/lowest_peak')
parser.add_argument('--sizealgo', '-a', nargs='?', default=0, type=int, choices=[0,1,2],help='Size algorithm selector. Algo 0: size, Algo 1: size+2*log(size*peak/(highest_peak/lowest_peak, Algo 2: size+2*size*peak/(highest_peak-lowest_peak)')
args = parser.parse_args()

# read input argument(s)
vkInput = getattr(args, "input")
identified = []
try:
  with open(vkInput, 'r') as tsv:
    next(tsv) # skip first row
    tsvData = csv.reader(tsv, delimiter='\t')
    for row in tsvData:
      identified.append([float(row[0]),row[1],float(row[2]),float(row[3]),row[4],float(row[5]),float(row[6]),float(row[7])])
except ValueError:
  print('The %s data file could not be loaded.' % vkInput)

vkOutput = getattr(args, "output")

vkPlotType = getattr(args, 'plottype')
if vkPlotType == '2d': vkPlotType = 'scatter-2d'
if vkPlotType == '3d': vkPlotType = 'scatter-3d'
print vkPlotType

vkSize = getattr(args, 'size')

vkSizeAlgo = getattr(args, 'sizealgo')

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
    py.plot(fig, filename=vkOutput+'.html')
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
    py.plot(fig, filename=vkOutput+'.html')

plotRatios(identified, vkPlotType)
