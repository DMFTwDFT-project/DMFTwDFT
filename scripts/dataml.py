#!/usr/bin/env python

import argparse
import os

import numpy as np
import pandas as pd


def store_data(args):
    """
	This stores the DMFT energies in an excel sheet.
	"""

    # creating dataframe
    df = pd.DataFrame(
        columns=["Configuration", "Etot (Migdal-Galisky)", "Etot (ctqmc sampling)"]
    )

    # iterating over folders
    pathlist = sorted([int(d) for d in os.listdir(args.path) if os.path.isdir(d)])
    print(pathlist)
    dirname = os.getcwd().split("/")[-1]

    for path in pathlist:

        pathstr_infotime = str(path) + os.sep + "DMFT" + os.sep + "INFO_TIME"
        pathstr_infoiter = str(path) + os.sep + "DMFT" + os.sep + "INFO_ITER"

        # first check if calculation is complete
        if os.path.exists(pathstr_infotime):
            fi = open(pathstr_infotime, "r")
            done_word = fi.readlines()[-1]
            fi.close()

            if done_word.split()[0] == "Calculation":

                # opening INFO_ITER if calculation is done
                fi = open(pathstr_infoiter, "r")
                lastline = fi.readlines()[-1]
                fi.close()

                etot1 = lastline.split()[6]
                etot2 = lastline.split()[7]

            else:
                print("Calculation incomplete.")
                etot1 = ""
                etot2 = ""

        else:
            print("Calculation incomplete.")
            etot1 = ""
            etot2 = ""

        # appending data to dataframe
        df = df.append(
            {
                "Configuration": path,
                "Etot (Migdal-Galisky)": etot1,
                "Etot (ctqmc sampling)": etot2,
            },
            ignore_index=True,
        )

    # store in spreadsheet
    filestr = "mldata_" + str(dirname) + ".xlsx"
    if os.path.exists(filestr):
        os.remove(filestr)
    df.to_excel(filestr, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script stores the DMFT energies in a spreadsheet."
    )
    parser.add_argument("path", type=str, default="./", help="Path to DMFT directory")
    args = parser.parse_args()
    store_data(args)
