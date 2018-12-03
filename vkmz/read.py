#!/usr/bin/env python
"""Convert input data into objects

Can either parse a tabular file or W4M's XCMS tabular as input.
Input MS data can be given in two "modes", (1) tabular or (2) Workflow4Metabolomics'
XCMS for Galaxy (W4M-XCMS) files.

Tabular mode requires a single tabular file as input and  must include the columns
"sample_name", "polarity", "mz", "rt", and "intensity". Each row represents a 
feature. Optionally a "charge" column can exist.

W4M-XCMS mode requires the sample metadata, variable metadata, and data matrix
files generated with W4M-XCMS. Feature charge infomration can be read from the
variable metadata file if it has been annotated with CAMERA.

If feature charge information is present, features without charge information
will be removed. If CAMERA annotation is present, only monoisotopic features
will be kept.
"""


import csv
import re
from vkmz.arguments import IMPUTE, POLARITY
from vkmz.objects import Sample, SampleFeatureIntensity, Feature


def polaritySanitizer(polarity):
    """Sanitize input polarity values.

    Renames the, case-insensitive, values 'positive', 'pos', or '+' to 'positive'
    and 'negative', 'neg', or '-' to 'negative'.

    Errors on unrecognized polarity value.

    Arguments:
        polarity (str): unsanitized polarity type
    """
    if polarity.lower() in ["positive", "pos", "+"]:
        polarity = "positive"
    elif polarity.lower() in ["negative", "neg", "-"]:
        polarity = "negative"
    else:
        raise ValueError(f"{polarity} is not recognized as a polarity type.")
    return polarity


def tabular(tabular_file):
    """Read a tabular file and create objects.

    Reads columns named "sample_name", "polarity", "mz", "rt", and "intensity"
    to create Sample, SampleFeatureIntensity, and Feature objects. If these
    required columns do not exist function will error.

    Optionally reads a charge column.

    Results in two name-object dictionaries for samples and features.

    Should check for feature name in tabular.

    Arguments:
        tabular_file (str): path to input tabular file
    """
    samples = {}
    features = {}
    try:
        with open(tabular_file, "r") as f:
            tabular_data = csv.reader(f, delimiter="\t")
            header = next(tabular_data)
            try:
                sample_name_index = header.index("sample_name")
                if not POLARITY:  # --polarity argument is not used
                    polarity_index = header.index("polarity")
                mz_index = header.index("mz")
                rt_index = header.index("rt")
                intensity_index = header.index("intensity")
                #
                charge_index = bool("charge" in header)
                if charge_index:
                    charge_index = header.index("charge")
            except ValueError:
                print(
                    """An expected column was not found in the tabular file.
                    The tabular file must contain columns named: "sample_name",
                    "polarity", "mz", "rt", "intensity", and "charge".'
                    """
                )
                raise
            for row in tabular_data:
                sample_name = row[sample_name_index]
                if POLARITY:  # --polarity argument is used
                    polarity = POLARITY
                else:
                    polarity = polaritySanitizer(row[polarity_index])
                mz = float(row[mz_index])
                rt = float(row[rt_index])
                feature_name = f"{polarity}-{rt}-{mz}"
                intensity = float(row[intensity_index])
                charge = None
                if charge_index:
                    charge = row[charge_index]
                    # charge data is null and impute flag is set
                    if not charge and IMPUTE:
                        charge = 1
                    else:
                        break
                if sample_name not in samples:
                    samples[sample_name] = Sample(sample_name)
                if feature_name not in features:
                    feature = Feature(
                        feature_name, sample_name, polarity, mz, rt, charge
                    )
                    features[feature_name] = feature
                else:
                    feature = features[feature_name]
                    feature.samples.append(sample_name)
                sfi = SampleFeatureIntensity(intensity, feature)
                samples[sample_name].sfis.append(sfi)
    except IOError:
        print(f"Error while reading {tabular_file}.")
        raise
    return samples, features


# TODO: break up function
def xcmsTabular(sample_file, variable_file, matrix_file):
    """Read W4M's XCMS tabular files and return a list of features.

    Reads sample metadata to create a dictionary of sample ids keys with,
    sanitized, polarity values.

    Reads variable metadata to create dictionaries, with feature names as keys,
    for mass-to-charge ratio, retention time, and charge. Mass-to-charge and
    retention time are stored together as a tuple.

    Finally, read data matrix and create all Feature objects and append to list.

    Arguments:
        sample_file (str): path to input sample metadata file
        variable_file (str): path to input variable metadata file
        matrix_file (str): path to input data matrix file
    """
    samples = {}
    features = {}
    # extract sample polarities
    try:
        polarity = {}
        with open(sample_file, "r") as f:
            sample_data = csv.reader(f, delimiter="\t")
            next(sample_data)  # skip header
            for row in sample_data:
                sample = row[0]
                if POLARITY:
                    polarity[sample] = POLARITY
                else:
                    polarity[sample] = polaritySanitizer(row[2])
    except IOError:
        print(f"Error while reading the XCMS tabular file {sample_file}")
        raise
    # extract variable mz & rt
    try:
        mz_rt = {}
        mz_index = int()
        rt_index = int()
        charges = {}
        with open(variable_file, "r") as f:
            variable_data = csv.reader(f, delimiter="\t")
            header = next(variable_data)
            mz_index = header.index("mz")
            rt_index = header.index("rt")
            isotopes_index = False
            if "isotopes" in header:
                isotopes_index = header.index("isotopes")
            for row in variable_data:
                feature_name = row[0]
                mz = float(row[mz_index])
                rt = float(row[rt_index])
                mz_rt[feature_name] = (mz, rt)
                if isotopes_index:
                    charges[feature_name] = row[isotopes_index]
                else:  # CAMERA / charge data does not exist
                    charges[feature_name] = None
    except IOError:
        print(f"Error while reading the XCMS tabular file {variable_file}.")
        raise
    # if CAMEARA data exists
    if isotopes_index:
        camera_pattern = re.compile(r"^\[\d+\]\[M\+1\]")
        for c in charges:
            charge = charges[c]
            monoisotopic = bool(camera_pattern.search(charge))
            if monoisotopic:
                # TODO: extract correct charge
                #       currently impute 1
                #       need multi-charge, + &, - test data
                # parse monoisotopic charge
                charges[c] = 1
            elif IMPUTE:
                charges[c] = 1
            else:
                charges[c] = "remove"
    # extract intensity and build Feature objects
    try:
        with open(matrix_file, "r") as f:
            matrix_data = csv.reader(f, delimiter="\t")
            header = next(matrix_data)  # list of samples
            # remove empty columns
            # required for W4M-XCMS 1.7
            # TODO: check W4M-XCMS 3.0
            header = [x for x in header if x is not ""]
            for row in matrix_data:
                # remove empty columns
                row = [x for x in row if x is not ""]
                feature_name = row[0]
                i = 1
                while i < len(row):
                    feature_charge = charges[feature_name]
                    # if CAMERA data exists, remove non-monoisotopic features
                    if feature_charge is "remove":
                        break
                    intensity = row[i]  # keep as string type for test
                    if intensity not in {"NA", "#DIV/0!", "0"}:
                        sample_name = header[i]
                        if sample_name not in samples:
                            samples[sample_name] = Sample(sample_name)
                        if feature_name not in features:
                            feature_polarity = polarity[sample_name]
                            feature_mz = mz_rt[feature_name][0]
                            feature_rt = mz_rt[feature_name][1]
                            intensity = float(intensity)
                            feature = Feature(
                                feature_name,
                                sample_name,
                                feature_polarity,
                                feature_mz,
                                feature_rt,
                                feature_charge,
                            )
                            features[feature_name] = feature
                        else:
                            feature = features[feature_name]
                            feature.samples.append(sample_name)
                        samples[sample_name].sfis.append(
                            SampleFeatureIntensity(intensity, feature)
                        )
                    i += 1
    except IOError:
        print(f"Error while reading the XCMS tabular file {matrix_file}.")
        raise
    return samples, features
