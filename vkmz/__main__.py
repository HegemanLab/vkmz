#!/usr/bin/env python


def main():
    """Main flow control of vkmz
  
    Read input data into feature objects. Results in dictionaries for samples
    and features.

    Then, make predictions for features. Features without predictions are removed
    by default.

    Finally, write results.
    """
    from vkmz.arguments import args, JSON, METADATA, MODE, SQL
    from vkmz.read import tabular as readTabular, xcmsTabular as readXcmsTabular, formulas as readFormulas
    from vkmz.predict import predict
    import vkmz.write as write

    # read input
    if MODE == "tabular":
        # read arguments here in case "input" is undeclared
        tabular_f = getattr(args, "input")
        samples, features = readTabular(tabular_f)
        print(features)
    elif MODE == "xcms":
        sample_f = getattr(args, "sample_metadata")
        variable_f = getattr(args, "variable_metadata")
        matrix_f = getattr(args, "data_matrix")
        samples, features = readXcmsTabular(sample_f, variable_f, matrix_f)
    else:  # MODE == "formula"
        formula_f = getattr(args, "input")
        samples, features = readFormulas(formula_f)

    if MODE == "tabular" or MODE == "xmcs":
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
    j_objs = write.generateJson(samples)
    if JSON:
        write.json_write(j_objs)
    write.html(j_objs)
    if SQL:
        write.sql(samples, features)
    if METADATA:
        write.metadata()


if __name__ == "__main__":
    main()
