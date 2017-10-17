from operator import itemgetter
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--load',      '-l', nargs='?', default='',     
               help='Load a previously generated ratio table. Set file path. Disabled by default.')
args = parser.parse_args()

# read load argument
vkdbLoad = getattr(args, "load")

# read output argument
#vkDBOutput = getattr(args, "output")

def loadDB(vkdbLoad):
  try:
    with open(vkdbLoad, 'r') as f:
      db = f.readlines()
      print(db)
      sorted(db,key=itemgetter(3))
      print(db)
  except ValueError:
    print('The %s data file could not be loaded.' % vkLoad)

#loadDB(vkdbLoad)

with open(vkdbLoad) as f:
  print sorted(list(csv.reader(f)), key=itemgetter(3))
