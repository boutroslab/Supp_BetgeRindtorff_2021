# Semi-Supervised Segmentation of Organoids (S3O)

This package is for the semi-supervised segmentation of fluorescent microscopy images of organoids. It is designed to be abort-safe, i.e. saves intermediate results. As a result, it requires a significant amount of disk space

## Data Structure

The package is meant to be flexible with regards to formats. The package supports a nested folder structure, i.e.

Base Directory
    * Folder 1
        * Picture 1
        * Picture 2
        * ...
    * Folder 2
        * ...

This is to reflect the usual setup of high-throughput screens of plates and wells. However, all images are treated equally within this package, so that all pictures could essentially be placed into one folder. Note that a subfolder is essential, i.e. images saved in the base directory are ignored.

## Altering code

The code is designed to be easy to alter. Each functional module, e.g. intensity segmentation, contains an API function (named accordingly). This should be the only function called from other modules to permit easy alteration of code. In the spirit of the code being abort-safe, all modules save their output to disk, i.e. there is absolutely no passing of results between modules.

