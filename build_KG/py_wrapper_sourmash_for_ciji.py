import os
import sourmash as sm
import argparse
import numpy as np
import pandas as pd


def local_variable():
    outfile = 'test'
    query_db = os.path.abspath("./test_pyapi_2db/gp1.db.zip")
    ref_db = os.path.abspath("./test_pyapi_2db/gp2.db.zip")
    k_list = [21,31,51]
    similarity_method = "CI"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PyAPI for sourmash matrix when database is too large. No abundance, CI and JI only.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-q', '--query_db', type=str, help="Path to file containing ALL query signatures")
    parser.add_argument('-r', '--ref_db', type=str, help="Path to file containing ALL ref signatures")
    parser.add_argument('-k', '--kvalue', type=str, help="K values to use, comma sep, e.g. 21,31,51")
    parser.add_argument('-o', '--outname', type=str, help="Output CSV file name", default="test")
    parser.add_argument('-m', '--similarity_method', type=str, help="CI or JI")

    ### read parameters
    args = parser.parse_args()
    query_db = args.query_db
    ref_db = args.ref_db
    k_list = [int(x) for x in args.kvalue.split(",")]
    outfile = args.outname
    similarity_method = args.similarity_method
    if similarity_method == "CI":
        print("Generating CI matrix")
    elif similarity_method == "JI":
        print("Generating JI matrix")
    else:
        raise Exception("similarity_method must be either 'CI' or 'JI'")

    for k_value in k_list:
        print("Processing k value %s" % k_value)
        out_name = outfile + "_k_" + str(k_value) + ".csv"

        # read ref: might be better to use list for ref as we will retrive values for len(query) times
        iter_ref = sm.load_file_as_signatures(ref_db, ksize=k_value)
        list_ref = [x for x in iter_ref]

        # colnames
        out_colname = [x.filename for x in list_ref]
        out_colname.insert(0, similarity_method + "_of_row_in_col")

        # initial outfile header, will overwrite if file exists
        with open(out_name, 'w') as fp:
            fp.write(",".join(out_colname) + '\n')

        # loop every row in query_db and generate 1 row for CI to append to file (so we get 1 row at a time)
        iter_query = sm.load_file_as_signatures(query_db, ksize=k_value)
        for sig in iter_query:
            temp_row_name = sig.filename
            # time-limiting step
            if similarity_method == "CI":
                temp_containment_record = [str(round(sig.contained_by(x), 4)) for x in list_ref]
            else:  #similarity_method == "JI"
                temp_containment_record = [str(round(sig.jaccard(x), 4)) for x in list_ref]
            temp_containment_record.insert(0, temp_row_name)
            # append new line to outfile
            with open(out_name, 'a') as fp:
                fp.write(",".join(temp_containment_record) + '\n')

        # end of loop


