from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

## NOT RUN

# initialise the entry point
ent = entry_point()

D = [1000,1000] # Number of features per view
M = len(D)      # Number of views
K = 5           # Number of factors
N = [100,100]   # Number of samples per group
G = len(N)      # Number of groups

data_dt = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz", sep="\t")
data_dt.head()

# defining model options
ent.set_data_options(
    scale_groups = False, 
    scale_views = False
)

ent.set_model_options(
    factors = 10, 
    spikeslab_weights = True, 
    ard_factors = True,
    ard_weights = True
)

ent.set_train_options(
    iter = 1000, 
    convergence_mode = "fast", 
    startELBO = 1, 
    freqELBO = 1, 
    dropR2 = 0.001, 
    gpu_mode = True, 
    verbose = False, 
    seed = 1
)

ent.build()

ent.run()

# Save the output
outfile = "/Users/ricard/data/mofaplus/hdf5/test.hdf5"
ent.save(outfile)

