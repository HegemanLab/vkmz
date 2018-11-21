#!/usr/bin/env python


class Sample(object):
    """Sample Class

  Represents a single LC-MS sample injection. Object stores the sample's name
  and a list of SampleFeatureIntensity objects from the sample.
  """

    def __init__(self, name):
        self.name = name
        self.sfis = []


class SampleFeatureIntensity(object):
    """SampleFeatureIntensity Class

  Represents a feature's intensity from a single Sample.
  """

    def __init__(self, intensity, feature):
        self.intensity = intensity
        self.feature = feature


class Feature(object):
    """Feature Class

  Represents a feature from LC-MS data.
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
    """Prediction Class

  Represents a prediction based on a given feature's mass, polarity and charge.

  element_count attribute is a dictionary of the formula with elemental symbols
  as keys.
  """

    def __init__(self, mass, formula, delta, element_count, hc, oc, nc):
        self.formula = formula
        self.mass = mass
        self.delta = delta
        self.element_count = element_count
        self.hc = hc
        self.oc = oc
        self.nc = nc
