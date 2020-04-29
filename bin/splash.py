#!/usr/bin/env python2

version = "1.0"
date = "April 23, 2020"


def welcome():

    art = """
 ____  __  __ _____ _____          ____  _____ _____
|  _ \|  \/  |  ___|_   _|_      _|  _ \|  ___|_   _|
| | | | |\/| | |_    | | \ \ /\ / / | | | |_    | |
| |_| | |  | |  _|   | |  \ V  V /| |_| |  _|   | |
|____/|_|  |_|_|     |_|   \_/\_/ |____/|_|     |_|

Python 2.x version.
    """

    print(art)
    print(
        "--- An open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages. ---"
    )
    print("\nVersion %s updated on %s\n" % (version, date))
    print(
        "Please cite:\nVijay Singh, Uthpala Herath, Benny Wah, Xingyu Liao, Aldo H. Romero, Hyowon Park,. DMFTwDFT: An open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages,. arXiv:2002.00068 [cond-mat.str-el].\n"
    )

    print(
        "----------------------------------------------------------------------------------------\n"
    )

    return
