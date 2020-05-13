#!/usr/bin/env python2

import sys

version = "1.1"
date = "May 11th, 2020"


def welcome():

    art = """
 ____  __  __ _____ _____          ____  _____ _____
|  _ \|  \/  |  ___|_   _|_      _|  _ \|  ___|_   _|
| | | | |\/| | |_    | | \ \ /\ / / | | | |_    | |
| |_| | |  | |  _|   | |  \ V  V /| |_| |  _|   | |
|____/|_|  |_|_|     |_|   \_/\_/ |____/|_|     |_|

    """
    pversion = ".".join(map(str, sys.version_info[:3]))
    print(art)
    print("Python 2.x version running on Python %s." % pversion)
    print(
        "\n--- An open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages. ---"
    )
    print("\nVersion %s updated on %s.\n" % (version, date))
    print(
        "Please cite:\nVijay Singh, Uthpala Herath, Benny Wah, Xingyu Liao, Aldo H. Romero, Hyowon Park,. DMFTwDFT: An open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages,. arXiv:2002.00068 [cond-mat.str-el].\n"
    )

    separator_art = """
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        """
    print(separator_art)

    return
