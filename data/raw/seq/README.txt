I put the data into the B110-Isilon3 partition (the one that is for protected data).
Those are in the directory valentini/amplicon_analysis/data and valentini/amplicon_analysis/data_03_2018
It makes sense not to have data scattered but those are protected data since are genomic data from patients (in theory) so I'm not sure if it makes sense to have them in an "unprotected" directory.
If you don't mind I will raise the question at the next data management meeting (next week) and ask their opinion about it.

It's located at smb://gpcf.ilmn-isilonx.dkfz-heidelberg.de/b110_data/B110_Isilon3/valentini/amplicon_analysis/data and smb://gpcf.ilmn-isilonx.dkfz-heidelberg.de/b110_data/B110_Isilon3/valentini/amplicon_analysis/data_03_2018.
The whole analysis and the Snakemake files I used to obtain it are also in the same amplicon_analysis/ directory.
