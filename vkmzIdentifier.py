import numpy as np
import pandas as pd
import re
import os
import time
import extractNeededElementalData
import processElementalData
import bmrbLookup as bmrb
import argparse
import multiprocessing
from multiprocessing import Pool
from functools import partial
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--database',  '-d', nargs='?', default='bmrb-light.csv', help='Select database(s).')
parser.add_argument('--input',     '-i', nargs='?', type=str, required=True,  help='Data file must be a csv with the columns mz, polarity, intensity, & rt.')
parser.add_argument('--output',    '-o', nargs='?', type=str, required=True,  help='Specify name of output file.')
parser.add_argument('--error',     '-e', nargs='?', type=int, default=5,      help='Error in PPM for identification.')
parser.add_argument('--multiprocessing', '-m', action='store_true',           help='Call to use multiprocessing. One process per core.')
args = parser.parse_args()

vkDatabase = getattr(args, "database")
lt = bmrb.getLookupTable('databases/' + vkDatabase)

vkInput = getattr(args, "input")
vkInputMzs = []
try:
  with open(vkInput, 'r') as tsv:
    next(tsv) # skip first line
    tsvData = csv.reader(tsv, delimiter='\t')
    for row in tsvData:
      vkInputMzs.append([float(row[0]),row[1],float(row[2]),float(row[3]),[],str])
except ValueError:
  print('The %s data file could not be loaded.' % vkInput)

vkOutput = getattr(args, "output")

vkError = getattr(args, "error")

vkMultiprocessing = getattr(args, "multiprocessing")

def buildRatios(vkInputMzs):
  # elements to plot
  elements = ['C', 'H', 'O', 'N']
  identified = []
  if vkMultiprocessing:
    try:
      pool = Pool()
      multiprocessMzsArgs = partial(multiprocessMzs, vkError)
      identified = pool.map(multiprocessMzsArgs, vkInputMzs)
      identified = [x for x in identified if x is not None]
    except Exception as e:
      print(str(e))
    finally:
      pool.close()
      pool.join()
  else:
    for centroid in vkInputMzs:
      identity = bmrb.getFormulaFromMass(bmrb.adjust(centroid[0], centroid[1]), lt, vkError)
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
        centroid[5] = identity
        identified.append(centroid)
        # this would be a good place to add unsaturation (2+2(carbons)+2(nitrogens)-hydrogens)/2
  return identified

def multiprocessMzs(vkError, centroid): # recieves a single Mz
  identity = bmrb.getFormulaFromMass(bmrb.adjust(centroid[0], centroid[1]),lt,vkError)
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
    centroid[5] = identity
    return centroid

# write vk identified as tsv file
# this file can be read by vkmzPlotter.py
def saveRatios(identified):
  filename = vkOutput+'.tsv'
  try:
    with open(filename, 'w') as f: 
      f.writelines(str("mz\tpolarity\tintensity\tretention time\tidentified structure\tH:C\tC:O\tC:N") + '\n')
      for ratio in identified:
        elements = ratio[4][0]
        f.writelines(str(ratio[0])+'\t'+str(ratio[1])+'\t'+str(ratio[2])+'\t'+str(ratio[3])+'\t'+str(ratio[5])+'\t'+str(float(elements['H'])/float(elements['C']))+'\t'+str(float(elements['O'])/float(elements['C']))+'\t'+str(float(elements['N'])/float(elements['C']))+'\n')
  except ValueError:
    print('"%s" could not be saved.' % filename)

saveRatios(buildRatios(vkInputMzs))
