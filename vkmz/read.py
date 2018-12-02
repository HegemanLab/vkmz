#!/usr/bin/env python
"""Convert input data into objects

Can either parse a tabular file or W4M's XCMS tabular as input.
"""


import csv
from vkmz.arguments import POLARITY
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
                charge = row[charge_index]
                if sample_name not in samples:
                    samples[sample_name] = Sample(sample_name)
                if feature_name not in features:
                    feature = Feature(feature_name, sample_name, polarity, mz, rt)
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


def xcmsTabular(sample_file, variable_file, matrix_file):
    """Read W4M's XCMS tabular files and return a list of features.

    Reads sample metadata to create a dictionary of sample ids keys with,
    sanitized, polarity values.

    Reads variable metadata to create two dictionaries with variable ids as keys
    and either mass to charge or retention time as values.

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
        with open(variable_file, "r") as f:
            variable_data = csv.reader(f, delimiter="\t")
            header = next(variable_data)
            mz_index = header.index("mz")
            rt_index = header.index("rt")
            for row in variable_data:
                mz_rt[row[0]] = (float(row[mz_index]), float(row[rt_index]))
    except IOError:
        print(f"Error while reading the XCMS tabular file {variable_file}.")
        raise
    # extract intensity and build Feature objects
    try:
        with open(matrix_file, "r") as f:
            matrix_data = csv.reader(f, delimiter="\t")
            header = next(matrix_data)  # list of samples
            # remove empty columns
            # required for W4M-XCMS 1.7, 3.0 not checked
            header = [x for x in header if x is not ""]
            for row in matrix_data:
                # remove empty columns
                row = [x for x in row if x is not ""]
                feature_name = row[0]
                i = 1
                while i < len(row):
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
