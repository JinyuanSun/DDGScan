#!/usr/bin/env python

from time import time
import utlis.io as io
from utlis import grape_phaseI

if __name__ == "__main__":
    args = io.Parser().get_args()
    grape_phaseI.main1(args)
    print("[INFO]: Ended at %s" % (time.ctime()))
