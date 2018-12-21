#!/usr/bin/env python
"""vkmz.predict module

Functions to predict a MS feature's molecular structure.
"""


import re
from vkmz.arguments import ALTERNATE, FORMULA, MASS, MASS_ERROR, MAX_MASS_INDEX, NEUTRAL
from vkmz.objects import Prediction

PROTON = 1.00727646677


def adjust(mass, polarity, charge):
    """Adjust a charged mass to a neutral mass.

    Observed charge mass is adjusted by adding or removing the mass of protons
    based on polarity and charge of feature.

    WARNING:
        If a feature's charge is not specified a charge of 1 is imputed.

    Arguments:
        mass (float): a charged feature's mass
        polarity (str): ionization mode
        charge (float): electric charge
    """
    # if charge is not given, impute 1
    if charge is None:
        charge = 1
    if polarity == "positive":
        mass -= PROTON * charge
    else:  # polarity == "negative"
        mass += PROTON * charge
    return mass


def predictInit(mass, uncertainty, left, right):
    """Search for a matching mass within the known-mass list.

    Uses binary search to match a given mass to a known mass within a given
    uncertainty. Upon match returns known mass index.

    If no match is made returns -1.

    Arguments:
        mass (float): observed neutral molecular mass in daltons
        uncertainty (float): mass error range in daltons
        left (int): left index of MASS list
        right (int): right index of MASS list
    """
    mid = int(((right - left) / 2) + left)
    if left <= mid <= right and mid <= MAX_MASS_INDEX:
        delta = mass - MASS[mid]
        if uncertainty >= abs(delta):
            return mid
        elif uncertainty > delta:
            return predictInit(mass, uncertainty, left, mid - 1)
        else:
            return predictInit(mass, uncertainty, mid + 1, right)
    return -1


def predictAll(mass, uncertainty, init_index):
    """Search for all matching masses within the known-mass list.

    Checks adjacent indexes from a given index of known-masses which are within a
    given uncertainty of a given mass.

    Arguments:
        mass (float): observed neutral mass
        uncertainty (float): mass error range
        init_index (int): initial index in MASS list to begin search
    """
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
    """Parse molecular formula by it's constituent elements.

    Parses formula into a list of element symbols and integers of element count.
    From this list a dictionary, element_count, is created and returned with
    element symbol keys and element count values.

    Ratios for H:C, O:C, N:C are calculated and also returned.

    Formula's must use proper, case-sensitive, chemical symbols.
        (e.g., Copper is 'Cu', not 'CU')

    Arguments:
        formula (string): molecular formula
    """
    formula_list = re.findall("[A-Z][a-z]?|[0-9]+", formula)
    element_count = {}
    hc = float()
    oc = float()
    nc = float()
    i = 0
    while i < len(formula_list):
        # create keys
        if formula_list[i] not in element_count:
            element_count[formula_list[i]] = 0
        # if there is only one of this element
        if i + 1 == len(formula_list) or formula_list[i + 1].isalpha():
            element_count[formula_list[i]] += 1
        else:
            element_count[formula_list[i]] += int(formula_list[i + 1])
            i += 1
        i += 1
    if "C" in element_count:
        if "H" in element_count:
            hc = element_count["H"] / float(element_count["C"])
        if "O" in element_count:
            oc = element_count["O"] / float(element_count["C"])
        if "N" in element_count:
            nc = element_count["N"] / float(element_count["C"])
    return element_count, hc, oc, nc


def predict(feature):
    """Make predictions for a feature.

    Reads a Feature as input and, if possible, returns it with a list of
    Prediction objects.

    A feature is assumed to be charged by default. The observed charged mass of
    the feature is converted to a neutral mass through adjust(). The --neutral
    flag disables adjustment.

    predictInit() returns an index for the MASS/FORMULA lists  if a known mass is
    within the mass error uncertainty of the observed, neutral, mass.  Features
    without a prediction are thrown out.

    On match, predictAll() searches for matches adjacent to the initial match.

    By default, features with multiple predictions are thrown out unless the
    --alternate flag is set. Alternate matches are sorted by absolute delta.

    For each match an element_count dictionary is parsed and elemental ratios are
    calculated.

    Prediction objects are made for each match and added to the features
    predictions list before returning the feature object.

    Arguments:
        feature (Feature): feature to make a prediction for
    """
    if NEUTRAL:
        mass = feature.mz
    else:
        mass = adjust(feature.mz, feature.polarity, feature.charge)
    # uncertainty is the mass error in parts per million
    uncertainty = mass * MASS_ERROR / 1e6
    init_index = predictInit(mass, uncertainty, 0, MAX_MASS_INDEX)
    if init_index != -1:
        matches = predictAll(mass, uncertainty, init_index)
        # remove feature if multiple predictions are made and --alternate not set
        if not ALTERNATE and len(matches) > 1:
            return
        for m in matches:
            delta = mass - MASS[m]  # check with Stephen
            formula = FORMULA[m]
            element_count, hc, oc, nc = parseFormula(formula)
            feature.predictions.append(
                Prediction(MASS[m], formula, delta, element_count, hc, oc, nc)
            )
        # sort by lowest absolute delta
        if len(matches) > 1:
            feature.predictions.sort(key=lambda m: abs(m.delta))
        return feature
    # no prediction was made
    return
