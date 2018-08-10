"""
usage:
$ python normtoxsec.py xsec infile.yoda outfile.yoda
"""

import yoda
from sys import argv

xsec = float(argv[1])
objs = yoda.read(argv[2])

xsec_orig = objs["/_XSEC"].points[0].x
scale = xsec/xsec_orig

for obj in objs.values():
    if hasattr(obj, "scaleW"):
        obj.scaleW(scale)

objs["/_XSEC"].scaleX(scale)

yoda.write(objs, argv[3])
