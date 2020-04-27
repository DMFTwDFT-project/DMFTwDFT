#!/usr/bin/env python2
import pyfiglet

version = "1.0"
date = "April 23, 2020"


def welcome():
    print(pyfiglet.figlet_format("DMFTwDFT"))
    print(
        "- An open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages.\n"
    )
    print("Version %s created on %s\n" % (version, date))
    print(
        "Please cite: Vijay Singh, Uthpala Herath, Benny Wah, Xingyu Liao, Aldo H. Romero, Hyowon Park,. DMFTwDFT: An open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages,. arXiv:2002.00068 [cond-mat.str-el].\n"
    )

    print("-------------------------------------------------------")
    print("Starting calculation...\n")

    return
