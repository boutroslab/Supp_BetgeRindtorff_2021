# PROMISE Organoid Inference

The organoid segmentation is provided here as a Docker image for inference (no training).

## Building the Docker image

To build the docker image, simply execute the shell script

```
./build-docker-image.sh
```

Both the `DNN_weights.h5` and `DNN_weights_3_channel.h5` files can be used, but as I don't have any organoid images available, I cannot test their outputs. To switch out the weights file being used, simply change line 47 in the file `organoid-segmentation/Dockerfile` from:

```
COPY weights/DNN_weights.h5 DNN_weights.h5
```

to

```
COPY weights/DNN_weights_3_channel.h5 DNN_weights.h5
```

and run the build script as above.

## Run inference

To obtain predictions, place all images to be analyzed into the folder `io/` and then run the shell script `./run-inference.sh`. The Docker image will write all output files into the same folder. The outputs will have the `"_segmentation.h5"` suffix appended to the file names.

The Docker image expects images to be stored as h5 files with the image data stored as an array with the shape `(number of fields, height, width, 3 channels)` under the key "images". The following code generates a random demo image

```
import numpy as np
import h5py

arr = np.random.randint(0, 256, (4, 2048, 2048, 3), dtype=np.uint8)
with h5py.File("./io/sample.h5", "w-") as h5handle:
    h5handle.create_dataset(name="images", data=arr, compression=3)
```

Note that the Docker image expects to run on a GPU-enabled system. To run it on a CPU system, the following changes must be made.

1. In the file `organoid-segmentation/Dockerfile` change line 39 from `tensorflow-gpu==2.6.2 \` to `tensorflow==2.6.2 \`
2. In the file `run-inference.sh` remove line 5 `--runtime nvidia \` entirely.

Tensorflow still raises some warnings about CUDA missing, so I can make no guarantees for whether the outputs will be correct.