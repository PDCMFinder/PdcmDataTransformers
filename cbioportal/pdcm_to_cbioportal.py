#!/usr/bin/env python
from utilities.cbio_class import cBioPortal
from utilities.other import *

def run(args):
    # Path to data # output path
    cBioPortal(args[1], args[2]).main()


if len(sys.argv) > 0:
    run(sys.argv)