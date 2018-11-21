#!/usr/bin/env python

from arguments import args, JSON, METADATA, MODE, SQL
import read
from predict import predict
import write

def main():
  '''main of vkmz
  
     Read input data into feature objects. Results in dictionary for samples
     and features.

     Make predictions for features; pruning features without predictions.
  '''
  if MODE == "tabular":
    # read arguments here incase "input" is undeclared
    tabular_f = getattr(args, "input")
    samples, features = read.tabular(tabular_f)
  else: # MODE == "xcms"
    sample_f = getattr(args, "sample_metadata")
    variable_f = getattr(args, "variable_metadata")
    matrix_f = getattr(args, "data_matrix")
    samples, features = read.xcmsTabular(sample_f, variable_f, matrix_f)
  # make predictions for all features
  features = {k: predict(v) for k, v in features.items()}
  # remove features without a prediction
  features = {k: v for k, v in features.items() if v is not None}
  # remove sample feature intensities without a feature
  for s in samples.values():
    s.sfis = [x for x in s.sfis if len(x.feature.predictions) > 0] 
  # remove samples without a sample feature intensity
  samples = {k: v for k, v in samples.items() if len(v.sfis) > 0}
  # write outputs
  write.tabular(samples)
  json = write.generateJson(samples)
  if JSON:
    write.json(json)
  write.html(json)
  if SQL:
    write.sql(samples, features)
  if METADATA:
    write.metadata()
  # NOTE: DEBUG FOO
  foo_s = 0
  foo_sfis = 0
  for s in samples.values():
    foo_s += 1
    for sfi in s.sfis:
      foo_sfis += 1
  print("s:\t\t", foo_s)
  print("sfis:\t\t", foo_sfis)
  foo_f = 0
  for f in features.values():
    foo_f += 1
  print("f:\t\t", foo_f)
  # NOTE: END DEBUG FOO

if __name__ == '__main__':
  main()
