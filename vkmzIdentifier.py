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
from functools import partial
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--database',  '-d', nargs='*', default='bmrb-light.csv',   help='Select database(s).')
parser.add_argument('--input',     '-i', nargs='?', type=str, required=True, help='Enter data source. Data file must be a comma seperated list containing four elements: mz,polarity,intensity,rt.')
parser.add_argument('--output',    '-o', action='store_true',       help='Call variable to ouput a ratio file.')
parser.add_argument('--error',     '-e', nargs='?', type=int, default=5, help='Error in PPM for identification.')
parser.add_argument('--multiprocessing', '-m', action='store_true', help='Call variable to use multiprocessing. One process per core.')
args = parser.parse_args()

vkDatabase = getattr(args, "database")
#lt = []
#print(type(lt))
#for database in vkDatabase:
#  lt = bmrb.getLookupTable('databases/' + database)
lt = bmrb.getLookupTable('databases/' + vkDatabase)

vkInput = getattr(args, "input")
vkInputMzs = []
try:
  with open(vkInput, 'rb') as csvfile:
    foo = csv.reader(csvfile)
    for row in foo:
      # each row of csvfile becomes a tuple element inside a list
      # placeholder elements are added to each tuple
      vkInputMzs.append([float(row[0]),row[1],float(row[2]),float(row[3]),[],str])
except ValueError:
  print('The %s data file could not be loaded.' % vkInput)

# read output argument
vkOutput = getattr(args, "output")

# read PPM error argument
vkError = getattr(args, "error")

# read multiprocessing argument
vkMultiprocessing = getattr(args, "multiprocessing")

def buildRatios(vkInputMzs):
  # elements to plot
  elements = ['C', 'H', 'O', 'N']
  if vkMultiprocessing:
    try:
      pool = Pool()
      multiprocessMzsArgs = partial(multiprocessMzs, vkError)
      identified = pool.map(multiprocessMzsArgs, vkInputMzs[index])
    except Exception as e:
      print(str(e))
    finally:
      pool.close()
      pool.join()
  else:
    identified = []
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
  if vkOutput: 
    saveRatios(identified)

def multiprocessMzs(vkError, inputMz): # recieves a single Mz
  return bmrb.getFormulaFromMass( bmrb.adjust( inputMz, lt, vkError ) )

# write vk ratios as csv file
def saveRatios(ratios):
  try:
    filename = 'ratios-' + time.strftime("%Y%m%d%H%M%S") + '.tsv'
    with open(filename, 'w') as f: 
      f.writelines(str("mz\tpolarity\tintensity\tretention time\tidentified structure\tH:C\tC:O\tC:N") + '\n')
      for ratio in ratios:
        elements = ratio[4][0]
        f.writelines(str(ratio[0])+'\t'+str(ratio[1])+'\t'+str(ratio[2])+'\t'+str(ratio[3])+'\t'+str(ratio[5])+'\t'+str(float(elements['H'])/float(elements['C']))+'\t'+str(float(elements['O'])/float(elements['C']))+'\t'+str(float(elements['N'])/float(elements['C']))+'\n')
  except ValueError:
    print('"%s" could not be saved.' % filename)

buildRatios(vkInputMzs)
