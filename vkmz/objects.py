#!/usr/bin/env python


class Sample(object):
    """A MS sample.

    Attributes:
        name (str): samples name
        sfis (list): SampleFeatureIntensity objects
    """

    def __init__(self, name):
        self.name = name
        self.sfis = []


class SampleFeatureIntensity(object):
    """A feature's intensity observed by a sample.

    Attributes:
        intensity (float): intensity of a feature for a sample
        feature (Feature): feature with the observed intensity
    """

    def __init__(self, intensity, feature):
        self.intensity = intensity
        self.feature = feature


class Feature(object):
    """A feature from a MS dataset.

    Attributes:
        name (str): feature's name
        samples (list): Samples feature was observed in
        polarity (str): ionization mode
        mz (float): observed mass-to-charge ratio
        rt (float): retention time
        charge (float): electric charge
        predictions (list): Prediction objects
    """

    def __init__(self, name, samples, polarity, mz, rt, charge=None):
        self.name = name
        self.samples = [samples]
        self.polarity = polarity
        self.mz = mz
        self.rt = rt
        self.predictions = []
        self.charge = charge


class Prediction(object):
    """Prediction of a feature.

    Attributes:
        mass (float): neutral mass of predicted formula
        formula (str): molecular formula
        delta (float): absolute differnece between neutral mass of observed and
                       predicted mass
        element_count (dict): molecular formula as element symbol keys with
                              element count values
        hc (float): hydrogen to carbon ratio
        oc (float): oxygen to carbon ratio
        nc (float): nitrogen to carbon ratio
    """

    def __init__(self, mass, formula, delta, element_count, hc, oc, nc):
        self.formula = formula
        self.mass = mass
        self.delta = delta
        self.element_count = element_count
        self.hc = hc
        self.oc = oc
        self.nc = nc
