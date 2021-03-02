"""
Usage:

python slimMCTTBBDecayed.py infile.yoda outfile.yoda
"""

import yoda

from sys import argv

infname = argv[1]
outfname = argv[2]


d = yoda.read(infname)

dnew = {}
for (k, v) in d.iteritems():
  if k.startswith("/RAW"):
    continue

  if "MCTTBBDecayed" not in k:
    dnew[k] = v
    continue

  spl = k.split("/")

  prefix = "/".join(spl[:-3])
  variation = spl[-3]
  region = spl[-2]
  histname = spl[-1]

  # keep all nominal histograms
  if variation == "nominal":
    newpath = "/".join([prefix, region, histname])
    v.setPath(newpath)
    dnew[newpath] = v

  # keep all "nodet" histograms
  elif variation == "nodet":
    newpath = "/".join([prefix, variation, region, histname])
    v.setPath(newpath)
    dnew[newpath] = v

  # remove weight variations for non-nominal histograms
  elif "[" in histname and "]" in histname:
    pass

  else:
    newpath = "/".join([prefix, region, histname]) + "[%s]" % variation
    v.setPath(newpath)
    dnew[newpath] = v

  continue


yoda.write(dnew, outfname)
