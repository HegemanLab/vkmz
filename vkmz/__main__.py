#!/usr/bin/env python



def main():
    """Main flow control of vkmz
  
    Read input data into feature objects. Results in dictionary for samples
    and features.

    Make predictions for features; pruning features without predictions.
    """
    import os
    from vkmz.arguments import parser
    from vkmz.read import tabular as readTabular, xcmsTabular as readXcmsTabular
    from vkmz.predict import predict
    import vkmz.write as write

    # create constants
    args = parser.parse_args()
    ALTERNATE = getattr(args, "alternate")
    DATABASE = getattr(args, "database")
    CHARGE = getattr(args, "charge")
    JSON = getattr(args, "json")
    MASS_ERROR = getattr(args, "error")
    METADATA = getattr(args, "metadata")
    MODE = getattr(args, "mode")
    NEUTRAL = getattr(args, "neutral")
    OUTPUT = getattr(args, "output")
    POLARITY = getattr(args, "polarity")
    PREFIX = getattr(args, "prefix")
    # TODO: check if not against PEP8
    if not PREFIX:
        PREFIX = os.path.abspath(os.path.dirname(__file__))
    SQL = getattr(args, "sql")
    # MASS and FORMULA are used as indexable dictionaries
    MASS = []
    FORMULA = []
    try:
        with open(os.path.join(PREFIX, DATABASE), "r") as tabular:
            next(tabular)  # skip header
            for row in tabular:
                mass, formula = row.split()
                MASS.append(float(mass))
                FORMULA.append(formula)
    except:
        print(f"An error occured while reading the {DATABASE} database.")
        raise
    MAX_MASS_INDEX = len(MASS) - 1

    # read input
    if MODE == "tabular":
        # read arguments here incase "input" is undeclared
        tabular_f = getattr(args, "input")
        samples, features = readTabular(tabular_f)
    else:  # MODE == "xcms"
        sample_f = getattr(args, "sample_metadata")
        variable_f = getattr(args, "variable_metadata")
        matrix_f = getattr(args, "data_matrix")
        samples, features = readXcmsTabular(sample_f, variable_f, matrix_f)

    # make predictions for all features
    features = {k: predict(v) for k, v in features.items()}

    # remove features without a prediction
    features = {k: v for k, v in features.items() if v is not None}
    # remove sample feature intensities without a feature
    for s in samples.values():
        s.sfis = [x for x in s.sfis if len(x.feature.predictions) > 0]
    # remove samples without a sample feature intensity
    samples = {k: v for k, v in samples.items() if len(v.sfis) > 0}

    # write results
    write.tabular(samples)
    json = write.generateJson(samples)
    if JSON:
        write.json(json)
    write.html(json)
    if SQL:
        write.sql(samples, features)
    if METADATA:
        write.metadata()

    # NOTE: DEBUG FOO
    foo_s = 0
    foo_sfis = 0
    for s in samples.values():
        foo_s += 1
        j = 0
        y = []
        for sfi in s.sfis:
            foo_sfis += 1
            j += 1
            y.append(sfi.feature)
            print(sfi.intensity)
        print(s.name, j, len(y), len(set(y)))
    print("s:\t", foo_s)
    print("sfis:\t", foo_sfis)
    foo_f = 0
    for f in features.values():
        foo_f += 1
    print("f:\t", foo_f)
    # NOTE: END DEBUG FOO


if __name__ == "__main__":
    main()
