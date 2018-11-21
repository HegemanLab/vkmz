#!/usr/bin/env python

import csv
import re
import sqlite3

from arguments import ALTERNATE, CHARGE, DATABASE, DIRECTORY, JSON, \
                      MASS_ERROR, METADATA, MODE, NEUTRAL, OUTPUT, \
                      POLARITY, SQL

def tabular(samples):
  '''Write results to tabular'''
  try: 
    # NOTE: DEBUG FOO
    i = 0
    with open(OUTPUT + '.tabular', 'w') as t_file:
      t_header = ('sample_name\tfeature_name\tpolarity\tmz\trt\tintensity\t'
                  'predicted_mass\tpredicted_delta\tpredicted_formula\t'
                  'predicted_element_count\tpredicted_hc\tpredicted_oc\t'
                  'predicted_nc\n')
      if ALTERNATE:
        t_header = t_header[:-1] + '\talternate_predictions\n'
      t_file.writelines(t_header)
      for s in samples.values():
        for sfi in s.sfis:
          # NOTE: DEBUG FOO
          i += 1
          f = sfi.feature
          p = f.predictions[0]
          t_row = (f'{f.name}\t{f.polarity}\t{f.mz}\t{f.rt}\t'
                   f'{sfi.intensity}\t{p.mass}\t{p.delta}\t{p.formula}\t'
                   f'{p.element_count}\t{p.hc}\t{p.oc}\t{p.nc}\n')
          if ALTERNATE and len(f.predictions) > 1:
            t_append = []
            for a in f.predictions[1:]:
              t_append.append((a.mass, a.formula, a.delta))
            t_row = t_row[:-1] + '\t' + str(t_append) + '\n'
          t_file.writelines(t_row)
      # NOTE: DEBUG FOO
      print("tabular i:", i)
  except IOError as error:
    print('IOError while writing tabular output')
    raise

def generateJson(samples):
  '''Convert results to JSON'''
  json = ''
  for s in samples.values():
    for sfi in s.sfis:
      f = sfi.feature
      p = f.predictions[0]
      j_element_count = ''
      for e in p.element_count:
        j_element_count += f'      "{e}": {p.element_count[e]},\n'
      # remove ',\n'
      j_element_count = j_element_count[:-2]
      j_element = (f'{{\n'
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
                   f'    "element_count": {{\n{j_element_count}\n'
                   f'    }},\n'
                   f'    "hc": {p.hc},\n'
                   f'    "oc": {p.oc},\n'
                   f'    "nc": {p.nc}\n'
                   f'  }}\n'
                   f'}},\n'
                  )
      if ALTERNATE and len(f.predictions) > 1:
        j_append = str()
        for a in f.predictions[1:]:
          j_append += (f'    {{\n'
                       f'      "mass": {a.mass},\n'
                       f'      "delta": {a.delta},\n'
                       f'      "formula": "{a.formula}",\n'
                       f'      "element_count": "{a.element_count}",\n'
                       f'      "hc": {a.hc},\n'
                       f'      "oc": {a.oc},\n'
                       f'      "nc": {a.nc}\n'
                       f'    }},\n'
                      )
        j_element = j_element[:-4] + ',\n  "alternate_predictions"'+ \
                    ': [\n' + j_append[:-2] + '\n  ]\n},\n'
      json += j_element
    json = json[:-2]
    return json

def json(json):
  '''Write results to JSON'''
  try: 
    with open(OUTPUT + '.json', 'w') as jsonFile:
      jsonFile.write(json)
  except IOError as error:
    print('IOError while writing JSON output: %s' % error.strerror)

def html(json):
  '''Write results to html webpage'''
  try: 
    with open(DIRECTORY+'d3.html', 'r', encoding='utf-8') as h_template,\
         open(OUTPUT+'.html', 'w', encoding='utf-8') as h_file:
      for line in h_template:
        line = re.sub('^var data.*$','var data = ['+json+']', line, flags=re.M)
        h_file.write(line)
  except IOError as error:
    print('IOError while writing HTML output or reading HTML template')
    raise

def metadata():
  '''Write VKMZ parameters to tabular file'''
  if METADATA:
    try: 
      with open(OUTPUT + '_metadata.tabular', 'w') as m_file:
        metadata = (f'Mode\tMass\tOutput\tJSON\tSQL\tPolarity\t'
                    f'Neutral\tDatabase\tDirectory\tCharge\n'
                    f'{MODE}\t{MASS_ERROR}\t{OUTPUT}\t{SQL}\t{POLARITY}\t'
                    f'{NEUTRAL}\t{DATABASE}\t{DIRECTORY}\t{CHARGE}\n'
                   )
        m_file.write(metadata)
    except IOError as error:
      print('IOError while writing metadata output: %s' % error.strerror)

def sql(samples, features):
  '''Write results to sqlit3 database'''
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
  s_sql = []
  i = 1 # unique Id
  for sample_name in samples.keys():
    s_sql.append((i, sample_name))
    i += 1
  c.executemany('INSERT INTO Sample VALUES (?, ?)', (s_sql))
  # add Feature and Prediction values
  f_sql = []
  p_sql = []
  i = 1
  for f in features.values():
    f_sql.append((i, f.name, f.polarity, f.mz, f.rt, f.charge))
    for p in f.predictions:
      p_sql.append((p.formula, p.mass, p.delta, str(p.element_count),
                    p.hc, p.oc, p.nc, i))
    i += 1
  c.executemany('INSERT INTO Feature \
                 VALUES (?, ?, ?, ?, ?, ?)',
                (f_sql))
  c.executemany('INSERT INTO Prediction(Formula, Mass, Delta, ElementCount, \
                 Hc, Oc, Nc, FeatureId) \
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                (p_sql))
  # add SampleFeatureIntensity values
  sfi_sql = []
  s_id = 1
  for sample in samples.values():
    for sfi in sample.sfis:
      f_id = list(features).index(sfi.feature.name) + 1
      sfi_sql.append((sfi.intensity, s_id,
                     f_id))
    s_id += 1
  c.executemany('INSERT INTO SampleFeatureIntensity(Intensity, SampleId, \
                 FeatureId) \
                 VALUES (?,  ?, ?)',
                (sfi_sql))
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
              )''')
    c.execute('INSERT INTO Metadata(Mode, MassError, Output, Json, Sql, \
               Polarity, Neutral, Database, Directory, Charge) \
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
              (MODE, MASS_ERROR, OUTPUT, JSON, SQL, POLARITY, NEUTRAL,
               DATABASE, DIRECTORY, CHARGE))
  con.commit()
  con.close()
