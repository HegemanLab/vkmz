import sys
import base64  # Imports a binary converter package
import struct
import gzip  # For file handling
import xml.parsers.expat
from MSScan import MS1Scan, MS2Scan
import time
import extractNeededElementalData
import processElementalData
import bmrbLookup as bmrb
import argparse
import pymzml

# this code is mostly unmaintained. It is better to feed vkmzDriver from a seperate mzml/mzxml parser such as xcms
# the data strucutre will likely change in future vkmz releases

parser = argparse.ArgumentParser()
parser.add_argument('--input',     '-i', nargs='*', required=True,          help='Enter full mzXML/mzML file paths. For multiple files seperate files with a space.')
parser.add_argument('--output',    '-o', action='store_true',               help='Call variable to ouput a ratio file.')
parser.add_argument('--threshold', '-t', nargs='?', default='10', type=int, help='Set threshold as a percent integer from 0 to 100. eg. 15 for 15%%.')
args = parser.parse_args()

# this scripts can concatenate multiple files into one
# if you want an output file for each input file run this script for each input individually
vkMZMLInput = getattr(args, "input")
for file in vkMZMLInput:
  if file.lower().endswith(('.mzml', '.mzxml')) == False:
    raise ValueError('Input file "%s" is not a mzML or mzXML file.' % (file))
    exit()

vkMZMLOutput = getattr(args, "output")

vkMZMLThreshold = getattr(args, "threshold")
if 0 <= vkMZMLThreshold <= 100:
  vkMZMLThreshold = vkMZMLThreshold * 0.01
else:
  raise ValueError("The given threshold, %i, is out of bounds." % (vkMZMLThreshold))
  exit()

'''
Class for dealing with mzXML files and handling the data. Particularly useful for decoding the
peak lists.  Found from: https://code.google.com/p/massspec-toolbox/source/browse/#svn/trunk/mzxml
and minimal updates made.
'''
# Defines a class of objects of MzXML type. Will initialize with the following values.
class MzXML():
    def __init__(self):
        self.msLevel = 0
        self.current_tag = ''
        self.tag_level = 0
        self.MS1_list = []
        self.MS2_list = []

    '''
    Function that decodes a single line. ***NOTE*** This is where files being run through a different
    mzXML converter broke down. MM_File_Conversion3 was used in testing to convert .raw files (not folders)
    to mzXML files. If this converter is used, this code should work but currently looking into ways to
    make this work with other, more broadly used converters. mzML is another file format that should work
    more consistently if this gives you issues.
    '''
    def decode_spectrum(self, line):
        decoded = base64.decodestring(line)
        # Determines unpack format which is specific to the type of data being examined
        tmp_size = len(decoded) / 4
        unpack_format1 = ">%dL" % tmp_size
        idx = 0
        # Declares list for mzs and intensities
        mz_list = []
        intensity_list = []
        # Loops through the decoded, unpacked line and breaks them apart into mz and intensity lists
        for tmp in struct.unpack(unpack_format1, decoded):
            tmp_i = struct.pack("I", tmp)
            tmp_f = struct.unpack("f", tmp_i)[0]
            if (idx % 2 == 0):
                mz_list.append(float(tmp_f))
            else:
                intensity_list.append(float(tmp_f))
            idx += 1
        # Returns the lists of intensities and mzs
        return mz_list, intensity_list

    def _start_element(self, name, attrs):
        # Increments the tag_level for the MzXML object and updates the current_tag
        self.tag_level += 1
        self.current_tag = name
        # If it's a precursorMz it adjusts accordingly
        if name == 'precursorMz':
            self.MS2_list[-1].precursor_intensity = float(attrs['precursorIntensity'])
            self.MS2_list[-1].precursor_charge = 0
            if attrs.has_key('precursorCharge'):
                self.MS2_list[-1].precursor_charge = int(attrs['precursorCharge'])
        # If the element being read in is a scan, checks what the scan level is and initializes a list of that type
        if name == 'scan':
            self.msLevel = int(attrs['msLevel'])
            if self.msLevel == 1:
                tmp_ms = MS1Scan()
                '''
                Note, the below and above code is critical and will need to be adjusted depending on how your machine
                does scans. This setting will impact numerous other factors so it is critical to get this right.
                This same change needs to be made on line 96 as called out in the other comment.
                The main thing to be aware of is that msLevel == 1 is used for negative scans and msLevel == 0
                is used for positive mode scans which may vary by machine. Positive and negative modes are also
                handled by double checking polarity but if your machine uses something other than msLevel 0 or 1
                this will still need to be changed.
                '''
            elif (self.msLevel == 0):
                tmp_ms = MS2Scan()
            else:
                print("What is it?", attrs)
                sys.exit(1)
            # Assigns attributes to their logical properties
            tmp_ms.id = int(attrs['num'])
            tmp_ms.peak_count = int(attrs['peaksCount'])
            tmp_ms.retention_time = float(attrs['retentionTime'].strip('PTS'))
            tmp_ms.low_mz = float(attrs['lowMz'])
            tmp_ms.high_mz = float(attrs['highMz'])
            tmp_ms.base_peak_mz = float(attrs['basePeakMz'])
            tmp_ms.base_peak_intensity = float(attrs['basePeakIntensity'])
            tmp_ms.total_ion_current = float(attrs['totIonCurrent'])
            tmp_ms.list_size = 0
            tmp_ms.encoded_mz = ''
            tmp_ms.encoded_intensity = ''
            tmp_ms.mz_list = []
            tmp_ms.intensity_list = []
            tmp_ms.polarity = attrs['polarity']
            # Adds the scan to the correct list of scans
            if self.msLevel == 1:
                self.MS1_list.append(tmp_ms)
            elif self.msLevel == 0:  # *************** changed ms level from == 2 to ==0**************
                self.MS2_list.append(tmp_ms)

    # Reduces the tag level, sets current_tag to '' and msLevel at 0
    def _end_element(self, name):
        self.tag_level -= 1
        self.current_tag = ''
        self.msLevel == 0

    def _char_data(self, data):
        # if self.current_tag == 'precursorMz':
        #     self.MS2_list[-1].precursor_mz = float(data)
        if self.current_tag == 'peaks':
            mz_list, intensity_list = self.decode_spectrum(data)
            mz_string = ''.join([struct.pack('>f', i) for i in mz_list])
            intensity_string = ''.join([struct.pack('>f', i) for i in intensity_list])
            if self.msLevel == 1:
                self.MS1_list[-1].list_size += len(mz_list)
                self.MS1_list[-1].encoded_mz += base64.encodestring(mz_string)
                self.MS1_list[-1].encoded_intensity += base64.encodestring(intensity_string)
                self.MS1_list[-1].mz_list += mz_list
                self.MS1_list[-1].intensity_list += intensity_list
            elif self.msLevel == 0:
                self.MS2_list[-1].list_size += len(mz_list)
                self.MS2_list[-1].encoded_mz += base64.encodestring(mz_string)
                self.MS2_list[-1].encoded_intensity += base64.encodestring(intensity_string)
                self.MS2_list[-1].mz_list = mz_list
                self.MS2_list[-1].intensity_list = intensity_list

    def parse_file(self, filename_xml):
        sys.stderr.write("Reading %s ... " % filename_xml)
        f_xml = open(filename_xml, 'r')
        if filename_xml.endswith('.gz'):
            f_xml = gzip.open(filename_xml, 'rb')
        content_list = []
        for line in f_xml:
            content_list.append(line)
        f_xml.close()
        expat = xml.parsers.expat.ParserCreate()
        expat.StartElementHandler = self._start_element
        expat.EndElementHandler = self._end_element
        expat.CharacterDataHandler = self._char_data
        expat.Parse("".join(content_list))
        sys.stderr.write("Done\n")

'''
Takes an mzML object which contains a list of intensities, a
list of mzs, and a value for the threshold intensity for filtering.

pymzml citation
Bald, T., Barth, J., Niehues, A., Specht, M., Hippler, M., and Fufezan, C. (2012) pymzML - Python module for high throughput bioinformatics on mass spectrometry data, Bioinformatics, doi: 10.1093/bioinformatics/bts066
'''


def process_mzs(file_name, threshold=.1):  # What fraction of the max intensity will be used for a threshold
    # Sets up mzML reader
    ms_run = pymzml.run.Reader(
            file_name,
            obo_version='1.2.0',
            extraAccessions=[
                ('MS:1000129', ['value']),
                ('MS:1000130', ['value']),
            ]
    )
    # For elements that clear the filter
    keepers_neg_mz = []
    keepers_pos_mz = []

    '''
    Loops through positive list and negative list and adds mz values to the keeper list when the intensity is
    above a threshold and catches any errors that are thrown (test data still had some issues with encoded data
    not decoding properly but the try and except handles that).
    '''
    # Variable to make sure loop continues even if an error is thrown
    # for each spectrum in a "run" (iteration)
    for spectrum in ms_run:
        try:
            # Calculate the max value and then generate a threshold based on that.
            maxS = max(spectrum.i)
            thresh = maxS * threshold
            # Looks at each peak intensity, and if it is past the threshold, adds it to the neg or pos list
            for peak in spectrum.peaks:
                if peak[1] > thresh:
                    # Check for polarity once intensity is determined to be high enough
                    negative_polarity = spectrum.get('MS:1000129', False)
                    if negative_polarity == '':
                        keepers_neg_mz.append(peak[0])
                    positive_polarity = spectrum.get('MS:1000130', False)
                    if positive_polarity == '':
                        keepers_pos_mz.append(peak[0])
        except:
            pass
    # Removes duplicates
    filtered_neg_mz = list(set(keepers_neg_mz))
    filtered_pos_mz = list(set(keepers_pos_mz))
    # Combines list where negatives are in the 0 position and positives in the 1
    combo_set = [filtered_neg_mz, filtered_pos_mz]
    return combo_set

'''
Churns through an mzXML object to generate lists of positive and negative
mz values.
'''
# Takes an mzXML object which contains a list of intensities, a
# list of mzs, and a value for the threshold intensity for filtering.
def process_mzs(mzXML_obj, threshold=.1):  # What fraction of the max intensity you want to use for a threshold
    # For elements that clear the filter
    keepers_neg_mz = []
    keepers_pos_mz = []
    # Loops through positive list and negative list and adds mz values to the keeper list when the intensity is
    # above a threshold
    for scan in mzXML_obj.MS1_list:
        max_peak = max(scan.intensity_list)
        thresh = max_peak * threshold
        i = 0
        # Look through each peak in each scan
        for peak in scan.intensity_list:
            # if the intensity is great enough
            if peak > thresh:
                # Determine the polarity
                peak_percent = (peak - thresh) / (max_peak - thresh)
                if scan.polarity == '-':
                    keepers_neg_mz.append((scan.mz_list[i], peak_percent))
                elif scan.polarity == '+':
                    keepers_pos_mz.append((scan.mz_list[i], peak_percent))
            i += 1
    # Second list of scans found in some mzXML objects
    for scan in mzXML_obj.MS2_list:
        thresh = max(scan.intensity_list) * threshold
        i = 0
        # Look through each peak in each scan
        for peak in scan.intensity_list:
            # if the intensity is great enough
            if peak > thresh:
                # Determine the polarity
                if scan.polarity == '-':
                    keepers_neg_mz.append((scan.mz_list[i], peak))
                elif scan.polarity == '+':
                    keepers_pos_mz.append((scan.mz_list[i], peak))
            i += 1
    # Removes duplicates
    filtered_neg_mz = list(set(keepers_neg_mz))
    filtered_pos_mz = list(set(keepers_pos_mz))
    # Combines list where negatives are in the 0 position and positives in the 1
    combo_list = [filtered_neg_mz, filtered_pos_mz]
    return combo_list

# parse input files
# foobar testing code
def dataParser(vkMZMLInput, vkMZMLThreshold, vkMZMLOutput):
  vkInputMzs = [[],[]]
  for file in vkMZMLInput:
    if file.lower().endswith('.mzxml'):
      mzXML = MzXML()
      mzXML.parse_file(file)
      vkInputMzsTemp = process_mzs(mzXML, threshold=vkMZMLThreshold)
      vkInputMzs[0] += vkInputMzsTemp[0] # this code looks redundant
      vkInputMzs[1] += vkInputMzsTemp[1]
    elif file.lower().endswith('.mzml'):   #FLAG, test this filetype
      vkInputMzs = process_mzs(file, threshold=vkMZMLThreshold)
    # Removes all duplicates from both neg and pos lists
    #   this may have been done in process_mzs
    vkInputMzs[0] = list(set(vkInputMzs[0]))
    vkInputMzs[1] = list(set(vkInputMzs[1]))
    # create tuple of mz, polarity, intensity and retention time
    #   for each element
    # list compresion to convert list elements to tuples
    #posValues = [(x,'pos',None,None) for x in vkInputMzs[0]]
    posValues = []
    for element in vkInputMzs[0]:
       posValues.append((element[0],'pos',element[1],None))
    #negValues = [(x,'neg',None,None) for x in vkInputMzs[1]]
    negValues = []
    for element in vkInputMzs[1]:
       posValues.append((element[0],'neg',element[1],None))
    # sort tuples by mass value
    vkInputMzs = sorted(posValues+negValues, key=(lambda x: x[0]))
    try:
      import csv
      #filename = file + time.strftime("-%Y%m%d%H%M%S") + '.csv'
      with open(file+'.csv', 'w') as out: 
        csv_out=csv.writer(out)
        csv_out.writerow(['mass','polarity','retention time','intensity'])
        for row in vkInputMzs:
          csv_out.writerow(row)
    except ValueError:
      return

dataParser(vkMZMLInput, vkMZMLThreshold, vkMZMLOutput)

