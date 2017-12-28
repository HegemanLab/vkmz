import numpy as np
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--input',    '-i', nargs='?', default='', required=True, help='Load a previously generated ratio table. Set file path. Disabled by default.')
parser.add_argument('--plottype', '-p', nargs='?', default=['scatter'], choices=['scatter', 'heatmap', '3d'], help='Set to "scatter", "heatmap", or "3d". Default is "scatter".')
args = parser.parse_args()

# read input argument(s)
vkInput = getattr(args, "input")
identified = []
try:
  with open(vkInput, 'r') as tsv:
    next(tsv) # skip first row
    tsvData = csv.reader(tsv, delimiter='\t')
    for row in tsvData:
      identified.append([float(row[0]),row[1],float(row[2]),float(row[3]),float(row[5]),float(row[6]),float(row[7])])
except ValueError:
  print('The %s data file could not be loaded.' % vkInput)

# read plottype argument
vkPlotType = getattr(args, 'plottype')

def plotRatios(identified, type):
  import pandas as pd
  import numpy as np
  import plotly as py
  import plotly.graph_objs as go
  from plotly import __version__
  import plotly.offline as py
  import plotly.graph_objs as go
  traces = []
  trace_count = 0
  lowest_peak = 10.0**10
  highest_peak = 0.0
  highest_rt = identified[-1][3]
  feature_rts =[]
  feature_peaks = []
  #feature_names = []
  x=[]
  y=[]
  z=[]
  feature_names = []
  for feature in identified:
    feature_peak = feature[2]
    if feature_peak > highest_peak:
      highest_peak = feature_peak
    elif feature_peak < lowest_peak:
      lowest_peak = feature_peak
  for feature in identified:
    # changes retention time from seconds to minutes
    feature_rts.append(feature[3]/60)
    feature_peak = feature[2]
    feature_peaks.append(10+20*(feature_peak/(highest_peak-lowest_peak)))
    #feature_formula = feature[4][0]
    #feature_name = ''
    #for i in feature_formula:
    #   feature_name+=i+str(feature_formula[i])
    #feature_names.append(feature_name)
    x.append(feature[5]) # Oxygen / Carbon
    y.append(feature[4]) # Hydrogen / Carbon
    # 3d data is given to all plots for added plot.ly functionality
    z.append(feature[6]) # Nitrogen / Carbon
  if type == '3d':
    feature_trace = go.Scatter3d(
      x = x,
      y = y,
      z = z,
      mode='markers',
      text=feature_names,
      marker=dict(
        size=feature_peaks,
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
    py.plot(fig, filename='simple-3d-scatter.html')
  if type == 'scatter':
    feature_trace = go.Scatter(
      x = x,
      y = y,
      mode='markers',
      text=feature_names,
      marker=dict(
        size=feature_peaks,
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
    py.plot(fig, filename='simple-3d-scatter.html')

plotRatios(identified, vkPlotType)
