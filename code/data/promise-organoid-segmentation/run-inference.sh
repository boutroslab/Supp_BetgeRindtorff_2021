#!/bin/bash
docker run --rm \
  --name promise-organoid-segmentation \
  --mount type=bind,source="$(pwd)"/io,target=/io \
  --runtime nvidia \
  promise:organoid-segmentation