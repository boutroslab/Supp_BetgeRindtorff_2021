from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

# initialise the entry point
ent = entry_point()

infile = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/data/processed/PhenotypeSpectrum/mofa_model_all_drugs_unscaled.hdf5"
ent.save(infile)
