#!/usr/bin/env python
# Created on: Mar 30, 2023 at 7:04:08 PM
__author__ = "Michael Cuoco"

from TElocal_Toolkit.TEindex import TEfeatures
import pickle

ind = TEfeatures()
ind.build(snakemake.input[0])

with open(snakemake.output[0], "wb") as f:
    pickle.dump(ind, f)
