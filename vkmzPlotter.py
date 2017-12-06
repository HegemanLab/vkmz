import numpy as np

def plotRatios(identified, type):
  import pandas as pd
  import plotly as py
  import plotly.graph_objs as go
  if type == '3d':
    from plotly import __version__
    #from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
    import plotly.offline as py
    import plotly.graph_objs as go
    import numpy as np
    traces = []
    trace_count = 0
    lowest_peak = 10.0**10  # ???
    highest_peak = 0.0      # ???
    highest_rt = identified[-1][3]  # rt is ordered
    print(highest_rt)
    feature_rts =[]
    feature_peaks =[]
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
      feature_rts.append(feature[3]/60) # turn into minutes
      feature_peak = feature[2]
      feature_peaks.append(10+20*(feature_peak/(highest_peak-lowest_peak)))
      feature_formula = feature[4][0]
      feature_name = ''
      for i in feature_formula:
         feature_name+=i+str(feature_formula[i])
      feature_names.append(feature_name)
      # Ratio Builder
      if feature_formula['O'] == 0:
        x.append(0)
      else:
        x.append(float(feature_formula["C"])/feature_formula["O"])
      if feature_formula['N'] == 0:
        y.append(0)
      else:
        y.append(float(feature_formula["C"])/feature_formula["N"])
      if feature_formula['H'] == 0:
        z.append(0) #what?
      else:
        z.append(float(feature_formula["H"])/feature_formula["C"])
      # intensity builder
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
    print(traces)
    layout = go.Layout(
      title="Van Krevelen Diagram", 
      scene = dict(
        xaxis= dict(
          title= 'Carbon to Oxygen Ratio',
          zeroline= False,
          gridcolor='rgb(183,183,183)',
          showline=True
        ),
        zaxis=dict(
          title= 'Hydrogen to Carbon Ratio',
          zeroline= False,
          gridcolor='rgb(183,183,183)',
          showline=True
        ),
        yaxis= dict(
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
    from plotly import __version__
    #from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
    import plotly.offline as py
    import plotly.graph_objs as go
    import numpy as np
    traces = []
    trace_count = 0
    lowest_peak = 10.0**10
    highest_peak = 0.0
    highest_rt = identified[-1][3]
    feature_rts =[]
    feature_peaks =[]
    x=[]
    y=[]
    feature_names = []
    for feature in identified:
      feature_peak = feature[2]
      if feature_peak > highest_peak:
        highest_peak = feature_peak
      elif feature_peak < lowest_peak:
        lowest_peak = feature_peak
    for feature in identified:
      print(feature)
      feature_rts.append(feature[3]/60) # turn into minutes
      feature_peak = feature[2]
      feature_peaks.append(10+20*(feature_peak/(highest_peak-lowest_peak)))
      feature_formula = feature[4][0]
      feature_name = ''
      for i in feature_formula:
         feature_name+=i+str(feature_formula[i])
      feature_names.append(feature_name)
      # Ratio Builder
      if feature_formula['O'] == 0:
        x.append(0)
      else:
        x.append(float(feature_formula["C"])/feature_formula["O"])
      if feature_formula['H'] == 0:
        y.append(0)
      else:
        y.append(float(feature_formula["H"])/feature_formula["C"])
      # intensity builder
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
    print(traces)
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
