#!/usr/bin/env python

import argparse
import os

import numpy as np
import pandas as pd
import statistics
from scipy.stats import sem


def store_data(args):
    """
    This stores the DMFT energies in an excel sheet.
    """

    # creating dataframe
    df = pd.DataFrame(
        columns=[
            "Configuration",
            "Avg Etot (Migdal-Galisky)",
            "Avg Etot (ctqmc sampling)",
            "E1_std",
            "E2_std",
            "E1_sem",
            "E2_sem",
        ]
    )

    # directory
    dirpath = os.getcwd()

    # iterating over folders
    pathlist = sorted(
        [int(d) for d in os.listdir(dirpath) if os.path.isdir(d) and d.isnumeric()]
    )
    print(pathlist)
    dirname = os.getcwd().split("/")[-1]

    for path in pathlist:

        pathstr_infotime = str(path) + os.sep + "DMFT" + os.sep + "INFO_TIME"
        pathstr_infoiter = str(path) + os.sep + "DMFT" + os.sep + "INFO_ITER"

        # Initializing variables
        etot1_list = []
        etot2_list = []

        # Check if calculation is complete
        if os.path.exists(pathstr_infotime):
            fi = open(pathstr_infotime, "r")
            done_word = fi.readlines()[-1]
            fi.close()

            # Check if calculation is completed or not
            if done_word.split()[0] == "Calculation":
                print("Calculation complete at %s." % path)
            else:
                print("Calculation incomplete at %s." % path)
                etot1 = ""
                etot2 = ""

        else:
            print("INFO_TIME does not exist at %s." % path)
            etot1 = ""
            etot2 = ""

        # Getting total energy
        if args.navg == 1:
            # opening INFO_ITER if exists
            if os.path.exists(pathstr_infoiter):
                fi = open(pathstr_infoiter, "r")
                lastline = fi.readlines()[-1]
                fi.close()

                etot1 = float(lastline.split()[6])
                etot2 = float(lastline.split()[7])
            else:
                print("INFO_ITER does not exist at %s" % path)

        else:
            # opening INFO_ITER if exists
            if os.path.exists(pathstr_infoiter):

                fi = open(pathstr_infoiter, "r")
                lastlines = fi.readlines()[-args.navg :]
                fi.close()

                for i in range(args.navg):
                    etot1_list.append(float(lastlines[i].split()[6]))
                    etot2_list.append(float(lastlines[i].split()[7]))

                # Averaging
                etot1 = sum(etot1_list) / len(etot1_list)
                etot2 = sum(etot2_list) / len(etot2_list)

            else:
                print("INFO_ITER does not exist at %s" % path)

        # appending data to dataframe
        df = df.append(
            {
                "Configuration": path,
                "Avg Etot (Migdal-Galisky)": etot1,
                "Avg Etot (ctqmc sampling)": etot2,
                "E1_std": statistics.stdev(etot1_list),
                "E2_std": statistics.stdev(etot2_list),
                "E1_sem": sem(etot1_list),
                "E2_sem": sem(etot2_list),
            },
            ignore_index=True,
        )

    # store in spreadsheet
    filestr = "DMFT-total-energy_" + str(dirname) + ".xlsx"
    if os.path.exists(filestr):
        os.remove(filestr)
    df.to_excel(filestr, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script stores the DMFT energies in a spreadsheet."
    )
    # parser.add_argument("path", type=str, default=".", help="Path to root directory.")
    parser.add_argument(
        "-navg", type=int, default="1", help="Number of last iterations to average."
    )
    args = parser.parse_args()
    store_data(args)
