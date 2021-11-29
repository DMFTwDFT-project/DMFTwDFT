#!/usr/bin/env python

import argparse
import os

import numpy as np
import pandas as pd
import statistics
from scipy.stats import sem
from statistics import StatisticsError


def store_data(args):
    """
    This stores the DMFT energies in an excel sheet.
    """

    if args.batch:
        # Batch job assuming folders are named as 1, 2, .... n for each DMFT calculation.

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
                    print("Calculation in-progress at %s." % path)
                    etot1 = ""
                    etot2 = ""

            else:
                # print("INFO_TIME does not exist at %s." % path)
                print("Calculation waiting at %s." % path)
                etot1 = ""
                etot2 = ""

            # Getting total energy
            if args.navg == 1:
                # opening INFO_ITER if exists
                if os.path.exists(pathstr_infoiter):
                    try:
                        fi = open(pathstr_infoiter, "r")
                        lastline = fi.readlines()[-1]
                        fi.close()

                        etot1 = float(lastline.split()[6])
                        etot2 = float(lastline.split()[7])
                    except ValueError:
                        etot1 = ""
                        etot2 = ""

                    # appending data to dataframe
                    df = df.append(
                        {
                            "Configuration": path,
                            "Avg Etot (Migdal-Galisky)": etot1,
                            "Avg Etot (ctqmc sampling)": etot2,
                        },
                        ignore_index=True,
                    )

                else:
                    pass
                    # print("INFO_ITER does not exist at %s" % path)

            else:
                # opening INFO_ITER if exists
                if os.path.exists(pathstr_infoiter):

                    try:
                        fi = open(pathstr_infoiter, "r")
                        lastlines = fi.readlines()[-args.navg :]
                        fi.close()

                        for i in range(args.navg):
                            etot1_list.append(float(lastlines[i].split()[6]))
                            etot2_list.append(float(lastlines[i].split()[7]))

                        # Averaging
                        etot1 = sum(etot1_list) / len(etot1_list)
                        etot2 = sum(etot2_list) / len(etot2_list)

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

                    except ValueError:

                        # appending data to dataframe
                        df = df.append(
                            {
                                "Configuration": path,
                                "Avg Etot (Migdal-Galisky)": "",
                                "Avg Etot (ctqmc sampling)": "",
                                "E1_std": "",
                                "E2_std": "",
                                "E1_sem": "",
                                "E2_sem": "",
                            },
                            ignore_index=True,
                        )

                else:
                    pass
                    # print("INFO_ITER does not exist at %s" % path)


        # store in spreadsheet
        filestr = "DMFT-total-energy_" + str(dirname) + ".xlsx"
        if os.path.exists(filestr):
            os.remove(filestr)
        df.to_excel(filestr, index=False)

    else:
        # Single DMFT calculation

        dirname = os.getcwd() #.split("/")[-1]

        infotime = str(dirname) + os.sep + "DMFT" + os.sep + "INFO_TIME"
        infoiter = str(dirname) + os.sep + "DMFT" + os.sep + "INFO_ITER"

        # Initializing variables
        etot1_list = []
        etot2_list = []

        # Check if calculation is complete
        if os.path.exists(infotime):
            fi = open(infotime, "r")
            done_word = fi.readlines()[-1]
            fi.close()

            # Check if calculation is completed or not
            if done_word.split()[0] != "Calculation":
                print("Calculation incomplete!")
                etot1 = ""
                etot2 = ""

        else:
            print("INFO_TIME does not exist!")
            etot1 = ""
            etot2 = ""

        # Getting total energy
        if args.navg == 1:
            # opening INFO_ITER if exists
            if os.path.exists(infoiter):
                fi = open(infoiter, "r")
                lastline = fi.readlines()[-1]
                fi.close()

                etot1 = float(lastline.split()[6])
                etot2 = float(lastline.split()[7])
            else:
                print("INFO_ITER does not exist!")
        else:
            # opening INFO_ITER if exists
            if os.path.exists(infoiter):

                fi = open(infoiter, "r")
                lastlines = fi.readlines()[-args.navg :]
                fi.close()

                for i in range(args.navg):
                    etot1_list.append(float(lastlines[i].split()[6]))
                    etot2_list.append(float(lastlines[i].split()[7]))

                # Averaging
                etot1 = sum(etot1_list) / len(etot1_list)
                etot2 = sum(etot2_list) / len(etot2_list)

            else:
                print("INFO_ITER does not exist!")

        # Printing information

        print("Avg Etot (Migdal-Galisky) : {:6.4f} eV".format(etot1))
        print("Avg Etot (ctqmc sampling) : {:6.4f} eV".format(etot2))
        try:
            print("E1_std : {:6.4f} eV".format(statistics.stdev(etot1_list)))
            print("E2_std : {:6.4f} eV".format(statistics.stdev(etot2_list)))
            print("E1_sem : {:6.4f} eV".format(sem(etot1_list)))
            print("E2_sem : {:6.4f} eV".format(sem(etot2_list)))
        except(StatisticsError):
            print("Statistical analysis requires at least two iterations!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script stores the DMFT energies in a spreadsheet."
    )
    # parser.add_argument("path", type=str, default=".", help="Path to root directory.")
    parser.add_argument(
        "-navg", type=int, default="1", help="Number of last iterations to average."
    )
    parser.add_argument(
        "-batch", action="store_true", help="Batch job for multiple DMFT calculations. Assumes DMFT calculations are in folders with sequential integer names."
    )

    args = parser.parse_args()
    store_data(args)
