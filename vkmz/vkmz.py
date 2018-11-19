#!/usr/bin/env python

import argparse
import csv
import math
import re
import sqlite3

parser = argparse.ArgumentParser()
sub_parser = parser.add_subparsers(help='Select mode:', dest='mode')

# Tabular mode arguments
parse_tabular = sub_parser.add_parser('tabular',
                                      help='Use tabular data as input.'
                                     )
parse_tabular.add_argument('--input', '-i',
                           required=True,
                           help='Path to tabular file.'
                          )

# XCMS-tabular mode arguments
parse_xcms = sub_parser.add_parser('xcms',
                                   help='Use XCMS data as input.'
                                  )
parse_xcms.add_argument('--data-matrix', '-xd',
                        required=True, nargs='?', type=str,
                        help='Path to XCMS data matrix file.'
                       )
parse_xcms.add_argument('--sample-metadata', '-xs',
                        required=True, nargs='?', type=str,
                        help='Path to XCMS sample metadata file.'
                       )
parse_xcms.add_argument('--variable-metadata', '-xv',
                        required=True, nargs='?', type=str,
                        help='Path to XCMS variable metadata file.'
                       )

# All modes arguments
for sub_parser in [parse_tabular, parse_xcms]:
  sub_parser.add_argument('--output', '-o',
                          required=True, nargs='?', type=str,
                          help='Specify output file path.'
                         )
  sub_parser.add_argument('--error', '-e',
                          required=True, nargs='?', type=float,
                          help='Mass error of MS data in parts-per-million.'
                         )
  sub_parser.add_argument('--json', '-j',
                          action='store_true',
                          help='Set JSON flag to save JSON output.'
                         )
  sub_parser.add_argument('--sql', '-s',
                          action='store_true',
                          help='Set SQL flag to save SQL output.'
                         )
  sub_parser.add_argument('--metadata', '-m',
                          action='store_true',
                          help='Set metadata flag to save the tools metadata.'
                         )
  sub_parser.add_argument('--database', '-db',
                          nargs='?', default='databases/bmrb-light.tsv',
                          help='Define database of known formula-mass pairs.'
                         )
  sub_parser.add_argument('--directory', '-dir',
                          nargs='?', default='', type=str,
                          help='Define tool directory path.'
                         ) # default to '.' ?
  sub_parser.add_argument('--polarity', '-p',
                          choices=['positive', 'negative'],
                          help='Set all polarities to positive or negative.'
                         )
  sub_parser.add_argument('--neutral', '-n',
                          action='store_true',
                          help='Set masses in input data are neutral.'
                         )
  sub_parser.add_argument('--alternate', '-a',
                          action='store_true',
                          help='Set to keep features with multiple predictions.'
                         )
  sub_parser.add_argument('--charge', '-c',
                          action='store_true',
                          help='Set if input data contains charge.'
                         )
args = parser.parse_args()

class Sample(object):
  '''Sample Class

  Represents a single sample injection on a LC-MS. Stores the sample's name and
  a list of SampleFeatureIntensity objects belonging to the sample.
  '''
  def __init__(self, name):
    self.name = name
    self.sample_feature_intensities = []

class SampleFeatureIntensity(object):
  '''SampleFeatureIntensity Class

  Represents a feature's intensity from a single Sample.
  '''
  def __init__(self, intensity, feature):
    self.intensity = intensity
    self.feature = feature

class Feature(object):
  '''Feature Class

  Represents a feature from LC-MS data.

  WARNING: If a feature's charge is not specified in the input data a charge of
  1 is assumed to make predictions.
  '''
  def __init__(self, name, samples, polarity, mz, rt, charge = None):
    self.name = name
    self.samples = [samples]
    self.polarity = polarity
    self.mz = mz
    self.rt = rt
    self.predictions = []
    self.charge = charge

class Prediction(object):
  '''Prediction Class

  Represents a prediction based on a given feature's mass, polarity and charge.

  element_count is a dictionary of the formula with elemental symbols as keys.
  '''
  def __init__(self, mass, formula, delta, element_count, hc, oc, nc):
    self.formula = formula
    self.mass = mass
    self.delta = delta
    self.element_count = element_count
    self.hc = hc
    self.oc = oc
    self.nc = nc

# parameter constants
MODE = getattr(args, 'mode')
MASS_ERROR = getattr(args, 'error')
OUTPUT = getattr(args, 'output')
JSON = getattr(args, 'json')
SQL = getattr(args, 'sql')
POLARITY = getattr(args, 'polarity')
ALTERNATE = getattr(args, 'alternate')
NEUTRAL = getattr(args, 'neutral')
DATABASE = getattr(args, 'database')
DIRECTORY = getattr(args, 'directory')

# DEV NOTE
# parameter constants not fully implemented
CHARGE = getattr(args, 'charge')
METADATA = getattr(args, 'metadata')

# generated constants
# MASS and FORMULA are used as indexable key-value pairs
MASS = []
FORMULA = []
try:
  with open(DIRECTORY+DATABASE, 'r') as tabular:
    next(tabular)
    for row in tabular:
      mass, formula = row.split()
      MASS.append(float(mass))
      FORMULA.append(formula)
except:
  print(f'An error occured while reading the {DATABASE} database.')
  raise
MAX_MASS_INDEX = len(MASS)-1

def polaritySanitizer(polarity: str):
  '''Sanitize input polarity values.

  Renames the, case-insensitive, values 'positive', 'pos', or '+' to 'positive'
  and 'negative', 'neg', or '-' to 'negative'.

  Errors on unrecognized polarity value.
  '''
  if polarity.lower() in ['positive', 'pos', '+']:
    polarity = 'positive'
  elif polarity.lower() in ['negative', 'neg', '-']:
    polarity = 'negative'
  else:
    raise ValueError(f'{polarity} is not recognized as a polarity type.')
  return polarity

def readTabular(tabular_file):
  '''Read a tabular file and create objects.

  Reads columns named "sample_name", "polarity", "mz", "rt", "intensity", and
  optionally "charge" to create Sample, SampleFeatureIntensity, and Feature
  objects.

  Results in two dictionaries which are name-object pairs for samples and
  objects.

  Should check for feature name in tabular.
  '''
  samples = {}
  features = {}
  try:
    with open(tabular_file, 'r') as f:
      tabular_data = csv.reader(f, delimiter='\t')
      header = next(tabular_data)
      try:
        sample_name_index = header.index('sample_name')
        polarity_index = header.index('polarity')
        mz_index = header.index('mz')
        rt_index = header.index('rt')
        intensity_index = header.index('intensity')
        if CHARGE:
          charge_index = header.index('charge')
      except ValueError:
        print('An expected column was not found in the tabular file.\n'
              'The tabular file must contain columns named: "sample_name",  '
              '"polarity", "mz", and "rt".\n'
              'The tabular file must include a "charge" column if the --charge'
              ' argument was used.'
             )
        raise
      for row in tabular_data:
        sample_name = row[sample_name_index]
        polarity = polaritySanitizer(row[polarity_index])
        mz = float(row[mz_index])
        rt = float(row[rt_index])
        intensity = float(row[intensity_index])
        feature_name = polarity + '-' + str(rt) + '-' + str(mz)
        if sample_name not in samples:
          samples[sample_name] = Sample(sample_name)
        if feature_name not in features:
          feature = Feature(feature_name, sample_name, polarity, mz, rt)
          features[feature_name] = feature
        else:
          feature = features[feature_name]
          feature.samples.append(sample_name)
        # CHARGE untested
        if CHARGE:
          feature.charge = row[charge_index]
        samples[sample_name].sample_feature_intensities.append(
                                              SampleFeatureIntensity(intensity,
                                                                     feature)
                                                             )
  except IOError:
    print(f'Error while reading {tabular_file}.')
    raise
  return samples, features

def readXcmsTabular(sample_file, variable_file, matrix_file):
  '''Read W4M's XCMS tabular files and return a list of features.

  Reads sample metadata to create a dictionary of sample ids keys with,
  sanitized, polarity values.

  Reads variable metadata to create two dictionaries with variable ids as keys
  and either mass to charge or retention time as values.

  Finally, read data matrix and create all Feature objects and append to list.
  '''
  samples = {}
  features = {}
  # extract sample polarities
  try:
    polarity = {}
    with open(sample_file, 'r') as f:
      sample_data = csv.reader(f, delimiter='\t')
      next(sample_data) # skip header
      for row in sample_data:
        sample = row[0]
        if POLARITY:
          polarity[sample] = POLARITY
        else:
          polarity[sample] = polaritySanitizer(row[2])
  except IOError:
    print(f'Error while reading the XCMS tabular file {sample_file}')
    raise
  # extract variable mz & rt
  try:
    mz_rt = {}
    mz_index = int()
    rt_index = int()
    with open(variable_file, 'r') as f:
      variable_data = csv.reader(f, delimiter = '\t')
      header = next(variable_data)
      mz_index = header.index('mz')
      rt_index = header.index('rt')
      for row in variable_data:
        mz_rt[row[0]] = (float(row[mz_index]), float(row[rt_index]))
  except IOError:
    print(f'Error while reading the XCMS tabular file {variable_file}.')
    raise
  # extract intensity and build Feature objects
  try:
    with open(matrix_file, 'r') as f:
      matrix_data = csv.reader(f, delimiter = '\t')
      header = next(matrix_data) # list of samples
      # remove empty columns
      # check w/ fresh input
      header = [x for x in header if x is not '']
      for row in matrix_data:
        # remove empty columns
        row = [x for x in row if x is not '']
        feature_name = row[0]
        i = 1
        while(i < len(row)):
          intensity = row[i] # keep type string for test
          if intensity not in {'NA', '#DIV/0!', '0'}:
            sample_name = header[i]
            if sample_name not in samples:
              samples[sample_name] = Sample(sample_name)
            if feature_name not in features:
              feature_polarity = polarity[sample_name]
              feature_mz = mz_rt[feature_name][0]
              feature_rt = mz_rt[feature_name][1]
              intensity = float(intensity)
              feature = Feature(feature_name, sample_name, feature_polarity,
                                feature_mz, feature_rt
                               )
              features[feature_name] = feature
            else:
              feature = features[feature_name]
              feature.samples.append(sample_name)
            samples[sample_name].sample_feature_intensities.append(
                                     SampleFeatureIntensity(intensity, feature)
                                                                 )
          i+=1
  except IOError:
    print(f'Error while reading the XCMS tabular file {matrix_file}.')
    raise
  return samples, features

def adjust(mass, polarity, charge):
  '''
  Adjust a charged mass to a neutral mass.

  Mass of an electrons (1.007276) is multiplied by the charge and subtracted
  from positively charged ions and added to negatively charged ions.

  WARNING: If a feature's charge is not specified in the input data a charge of
  1 is assumed.
  '''
  # value to adjust by
  proton = 1.007276
  # if charge is not given, assume 1
  if charge == None:
    charge = 1
  if polarity == 'positive':
    mass -= (proton * charge)
  else: # sanitized to negative
    mass += (proton * charge)
  return mass

def predict(mass, uncertainty, left, right):
  '''Search for a matching mass within the known-mass list.

  Uses binary search to match a given mass to a known mass within a given
  uncertainty. Upon match returns known mass index.

  If no match is made returns -1.
  '''
  mid = int(((right - left) / 2) + left)
  if left <= mid <= right and mid <= MAX_MASS_INDEX:
    delta = mass - MASS[mid]
    if uncertainty >= abs(delta):
      return mid
    elif uncertainty > delta:
      return predict(mass, uncertainty, left, mid-1)
    else:
      return predict(mass, uncertainty, mid+1, right)
  return -1

def predictAll(mass, uncertainty, init_index):
  '''Search for all matching masses within the known-mass list.

  Checks adjacent indexes from a given index of known-masses which are within a
  given uncertinty of a given mass.
  '''
  matches = [init_index]
  i = 0
  while init_index + i + 1 <= MAX_MASS_INDEX:
    m = init_index + i + 1
    delta = mass - MASS[m]
    if uncertainty >= abs(delta):
      matches.append(m)
      i += 1
    else:
      break
  i = 0
  while init_index + i - 1 >= 0:
    m = init_index + i - 1
    delta = float(MASS[m]) - mass
    if uncertainty >= abs(delta):
      matches.append(m)
      i -= 1
    else:
      break
  return matches

def parseFormula(formula):
  '''Parse molecular formula by it's constituent elements.

  Parses formula into a list of element symbols and integers of element count.
  From this list a dictionary, element_count, is created and returned with
  element symbol keys and element count values.

  Ratios for H:C, O:C, N:C are calculated and also returned.
  '''
  formula_list = re.findall('[A-Z][a-z]?|[0-9]+', formula)
  element_count = {}
  hc = float()
  oc = float()
  nc = float()
  i = 0;
  while i < len(formula_list):
    # create keys
    if formula_list[i] not in element_count:
      element_count[formula_list[i]] = 0
    # if there is only one of this element
    if i+1 == len(formula_list) or formula_list[i+1].isalpha():
      element_count[formula_list[i]] += 1
    else:
      element_count[formula_list[i]] += int(formula_list[i + 1])
      i += 1
    i += 1
  if 'C' in element_count:
    if 'H' in element_count : hc = element_count['H']/float(element_count['C'])
    if 'O' in element_count : oc = element_count['O']/float(element_count['C'])
    if 'N' in element_count : nc = element_count['N']/float(element_count['C'])
  return element_count, hc, oc, nc

def featurePrediction(feature):
  '''Make predictions for a feature.

  Reads a Feature as input and, if possible, returns it with a list of
  Prediction objects.

  A feature is assumed to be charged by default. The observed charged mass of
  the feature is converted to a neutral mass through adjust(). The --neutral
  flag disables adjustment.

  predict() returns an index for the MASS/FORMULA lists  if a known mass is
  within the mass error uncertainty of the observed, neutral, mass.  Features
  without a prediction are thrown out.

  On match, predictAll() searches for matches adjacent to the initial match.

  By default, features with multiple predictions are thrown out unless the
  --alternate flag is set. Alternate matches are sorted by absolute delta.

  For each match an element_count dictionary is parsed and elemental ratios are
  calculated.

  Prediction objects are made for each match and added to the features
  predictions list before returning the feature object.
  '''
  if NEUTRAL:
    mass = feature.mz
  else:
    mass = adjust(feature.mz, feature.polarity, feature.charge)
  # uncertainty is the mass error in parts per million
  uncertainty = mass * MASS_ERROR / 1e6
  init_index = predict(mass, uncertainty, 0, MAX_MASS_INDEX) 
  if init_index != -1:
    matches = predictAll(mass, uncertainty, init_index)
    # remove feature if multiple predictions are made and --alternate not set
    if not ALTERNATE and len(matches) > 1:
      return
    for m in matches:
      delta = mass - MASS[m] # check with Stephen
      formula = FORMULA[m]
      element_count, hc, oc, nc = parseFormula(formula)
      feature.predictions.append(Prediction(MASS[m], formula, delta,\
                                 element_count, hc, oc, nc)
                                )
    # sort alternate matches by lowest absolute delta
    if not ALTERNATE and len(matches) > 1:
      feature.predictions.sort(key=lambda m: abs(m.delta))
    return(feature)
  # no prediction was made
  return
 
# write output file
def write(samples):
  json = ''
  try: 
    with open(OUTPUT + '.tabular', 'w') as tabular_file:
      tabularHeader = ('sample_name\t'
                       'feature_name\t'
                       'polarity\t'
                       'mz\t'
                       'rt\t'
                       'intensity\t'
                       'predicted_mass\t'
                       'predicted_delta\t'
                       'predicted_formula\t'
                       'predicted_element_count\t'
                       'predicted_hc\t'
                       'predicted_oc\t'
                       'predicted_nc\n'
                      )
      if ALTERNATE:
        tabularHeader = tabularHeader[:-1] + '\talternate_predictions\n'
      tabular_file.writelines(tabularHeader)
      for s in samples.values():
        for sfi in s.sample_feature_intensities:
          f = sfi.feature
          p = f.predictions[0]
          tabular_row = (f'{f.name}\t{f.polarity}\t{f.mz}\t{f.rt}\t'
                         f'{sfi.intensity}\t{p.mass}\t{p.delta}\t'
                         f'{p.formula}\t{p.element_count}\t{p.hc}\t'
                         f'{p.oc}\t{p.nc}\n'
                        )
          json_element_count = ''
          for e in p.element_count:
            json_element_count += f'      "{e}": {p.element_count[e]},\n'
          json_element_count = json_element_count[:-2]
          json_element = (f'{{\n'
                          f'  "sample_name": "{s.name}",\n'
                          f'  "feature_name": "{f.name}",\n'
                          f'  "polarity": "{f.polarity}",\n'
                          f'  "mz": {f.mz},\n'
                          f'  "rt": {f.rt},\n'
                          f'  "intensity": {sfi.intensity},\n'
                          f'  "prediction": {{\n'
                          f'    "mass": {p.mass},\n'
                          f'    "delta": {p.delta},\n'
                          f'    "formula": "{p.formula}",\n'
                          f'    "element_count": {{\n{json_element_count}\n'
                          f'    }},\n'
                          f'    "hc": {p.hc},\n'
                          f'    "oc": {p.oc},\n'
                          f'    "nc": {p.nc}\n'
                          f'  }}\n'
                          f'}},\n'
                         )
          if ALTERNATE and len(f.predictions) > 1:
            tabular_append = []
            json_append = str()
            for a in f.predictions[1:]:
              tabular_append.append((a.mass, a.formula, a.delta))
              json_append += (f'    {{\n'
                              f'      "mass": {a.mass},\n'
                              f'      "delta": {a.delta},\n'
                              f'      "formula": "{a.formula}",\n'
                              f'      "element_count": "{a.element_count}",\n'
                              f'      "hc": {a.hc},\n'
                              f'      "oc": {a.oc},\n'
                              f'      "nc": {a.nc}\n'
                              f'    }},\n'
                             )
            tabular_row = tabular_row[:-1] + '\t' + str(tabular_append) + '\n'
            json_element = json_element[:-4] + ',\n  "alternate_predictions"'+ \
                           ': [\n' + json_append[:-2] + '\n  ]\n},\n'
          tabular_file.writelines(tabular_row)
          json += json_element
    json = json[:-2] # remove final comma # [:-1] ??
  except IOError as error:
    print('IOError while writing tabular output')
    raise
  try: 
    with open(DIRECTORY+'d3.html', 'r', encoding='utf-8') as htmlTemplate,\
         open(OUTPUT+'.html', 'w', encoding='utf-8') as htmlFile:
      for line in htmlTemplate:
        line = re.sub('^var data.*$','var data = ['+json+']', line, flags=re.M)
        htmlFile.write(line)
  except IOError as error:
    print('IOError while writing HTML output or reading HTML template')
    raise
  if METADATA:
    try: 
      with open(OUTPUT + '_metadata.tabular', 'w') as metadataFile:
        metadata = (f'Mode\tMass\tOutput\tJSON\tSQL\tPolarity\t'
                    f'Neutral\tDatabase\tDirectory\tCharge\n'
                    f'{MODE}\t{MASS_ERROR}\t{OUTPUT}\t{SQL}\t{POLARITY}\t'
                    f'{NEUTRAL}\t{DATABASE}\t{DIRECTORY}\t{CHARGE}\n'
                   )
        metadataFile.write(metadata)
    except IOError as error:
      print('IOError while writing metadata output: %s' % error.strerror)
  if JSON:
    try: 
      with open(OUTPUT + '.json', 'w') as jsonFile:
        jsonFile.write(json)
    except IOError as error:
      print('IOError while writing JSON output: %s' % error.strerror)
  if SQL:
    con = sqlite3.connect(OUTPUT + '.db')
    c = con.cursor()
    # create tables
    c.execute('''
      CREATE TABLE Sample (
        Id INTEGER PRIMARY KEY,
        Name TEXT
      )
      ''')
    c.execute('''
      CREATE TABLE Feature (
        Id INTEGER PRIMARY KEY,
        Name TEXT,
        Polarity TEXT,
        Mz REAL,
        Rt REAL,
        Charge INTEGER
      )
      ''')
    c.execute('''
      CREATE TABLE Prediction (
        Id INTEGER PRIMARY KEY AUTOINCREMENT,
        Formula TEXT,
        Mass TEXT,
        Delta REAL,
        ElementCount TEXT,
        Hc REAL,
        Oc REAL,
        Nc REAL,
        FeatureId INTEGER,
        FOREIGN KEY(FeatureId) REFERENCES Feature(Id)
      )
      ''')
    c.execute('''
      CREATE TABLE SampleFeatureIntensity (
        Id INTEGER PRIMARY KEY AUTOINCREMENT,
        Intensity REAL,
        SampleId INTEGER,
        FeatureId INTEGER,
        FOREIGN KEY(SampleId) REFERENCES Sample(Id),
        FOREIGN KEY(FeatureId) REFERENCES Feature(Id)
      )
      ''')
    # add Sample values
    sample_sql = []
    i = 1 # unique Id
    for sample_name in samples.keys():
      sample_sql.append((i, sample_name))
      i += 1
    c.executemany('INSERT INTO Sample VALUES (?, ?)', (sample_sql))
    # add Feature and Prediction values
    feature_sql = []
    prediction_sql = []
    i = 1
    for f in features.values():
      feature_sql.append((i, f.name, f.polarity, f.mz, f.rt, f.charge))
      for p in f.predictions:
        prediction_sql.append((p.formula, p.mass, p.delta, str(p.element_count),
                               p.hc, p.oc, p.nc, i)
                             )
      i += 1
    c.executemany('INSERT INTO Feature \
                   VALUES (?, ?, ?, ?, ?, ?)',
                  (feature_sql)
                 )
    c.executemany('INSERT INTO Prediction(Formula, Mass, Delta, ElementCount, \
                   Hc, Oc, Nc, FeatureId) \
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                  (prediction_sql)
                 )
    # add SampleFeatureIntensity values
    sample_feature_intensity_sql = []
    sample_key = 1
    for sample in samples.values():
      for sfi in sample.sample_feature_intensities:
        feature_key = list(features).index(sfi.feature.name) + 1
        sample_feature_intensity_sql.append(
                                            (sfi.intensity, sample_key,
                                             feature_key)
                                           )
      sample_key += 1
    c.executemany('INSERT INTO SampleFeatureIntensity(Intensity, SampleId, \
                   FeatureId) \
                   VALUES (?,  ?, ?)',
                  (sample_feature_intensity_sql)
                 )
    if METADATA:
      # add Metadata table and values
      c.execute('''
        CREATE TABLE Metadata (
          Mode,
          MassError,
          Output,
          Json,
          Sql,
          Polarity,
          Neutral,
          Database,
          Directory,
          Charge
        )
        ''')
      c.execute('INSERT INTO Metadata(Mode, MassError, Output, Json, Sql, \
                 Polarity, Neutral, Database, Directory, Charge) \
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                (MODE, MASS_ERROR, OUTPUT, JSON, SQL, POLARITY, NEUTRAL,
                 DATABASE, DIRECTORY, CHARGE)
               )
    con.commit()
    con.close()

# read input
if MODE == "tabular":
  tabular_file = getattr(args, "input")
  samples, features = readTabular(tabular_file)
else: # MODE == "xcms"
  sample_file = getattr(args, "sample_metadata")
  variable_file = getattr(args, "variable_metadata")
  matrix_file = getattr(args, "data_matrix")
  samples, features = readXcmsTabular(sample_file, variable_file, matrix_file)

# make predictions
features = {k: featurePrediction(v) for k, v in features.items()}
features = {k: v for k, v in features.items() if v is not None}

# remove sample features without predictions
for sample in samples.values():
  sample.sample_feature_intensities = [x for x in \
                                       sample.sample_feature_intensities if
                                       len(x.feature.predictions) > 0
                                      ] 
 
# sort by intensity so that D3 draws largest symbols first
# removed after object implementation
#predicted_list.sort(key=lambda x: x.intensity, reverse=True)

write(samples)
