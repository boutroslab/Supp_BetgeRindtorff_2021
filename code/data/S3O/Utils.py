# This module contains utility functions that don't fit elsewhere
import h5py


def write_log(msg, filename=None, mode="a"):
    """
    Wrapper function to write log entries.
    :param msg: String. The message to write
    :param filename: (Optional) String. The full file path of a log file to
    write to. If this is missing, the string will be printed to the standard
    output.
    :param mode: The 'open()' mode. Can be set to "w" if the log file should
    be overwritten.
    :return:
    """

    if filename is None:
        print(msg)
        return None
    else:
        with open(filename, mode=mode) as f:
            f.write(msg + "\n")


def merge_hdf5_files(in_files, out_file, use_compression=True):
    """
    Merges multiple hdf5 files together into a single file. All in_files must
    have the same datasets

    Merging means that the contents of each dataset are concatenated along
    the first axis. Other axes are currently not supported

    If the hdf5 file is enormous, then compression can be beneficial.
    Keep in mind, however, that this increases reading times during training.

    :param in_files: An iterable of full filenames (Strings)
    :param out_file: String
    :param use_compression: Boolean
    :return:
    """

    # This is here because the code is ALMOST generalizable to any axis
    axis = 0

    # Test datasets
    data_sets = []
    shapes = []
    dtypes = []
    for in_file in in_files:
        with h5py.File(in_file, "r") as h5handle:
            keys = tuple(sorted(h5handle.keys()))
            data_sets.append(keys)
            shapes.append([h5handle[key].shape for key in keys])
            dtypes.append([h5handle[key].dtype for key in keys])

    if len(set(data_sets)) != 1:
        print("'in_files' have differing datasets")
        return None

    keys = data_sets[0]
    num_datasets = set([len(ds) for ds in data_sets]).pop()

    # Check that the dimensions of each dataset are identical (with the
    # exception of 'axis')
    for ii in range(num_datasets):
        sub_shapes = [list(shape[ii]) for shape in shapes]
        for sub_shape in sub_shapes:
            sub_shape.pop(axis)
        sub_shapes = [tuple(sub_shape) for sub_shape in sub_shapes]
        if len(set(sub_shapes)) != 1:
            print("Dataset '{}' not identical in shape across files".format(
                keys[ii]))
            return None

    # Check that dtypes are identical
    for ii in range(num_datasets):
        sub_dtypes = [dtype[ii] for dtype in dtypes]
        if len(set(sub_dtypes)) != 1:
            print("Dataset '{}' not identical in dtype across files".format(
                keys[ii]))
            return None
    dtype = dtypes[0]

    # Count total number of entries for each dataset
    num_entries = [0] * len(keys)
    for ii in range(num_datasets):
        num_entries[ii] += sum([shape[ii][axis] for shape in shapes])

    # Create file
    with h5py.File(out_file, "w-") as h5handle:
        for ii in range(len(keys)):
            ds_shape = list(shapes[0][ii])
            ds_shape[axis] = num_entries[ii]
            if use_compression:
                h5handle.create_dataset(
                    name=keys[ii], shape=ds_shape,
                    dtype=dtype[ii], compression=3)
            else:
                h5handle.create_dataset(
                    name=keys[ii], shape=ds_shape,
                    dtype=dtype[ii])

    # Write to file
    # TODO: Update to generalize for any concat axis
    #  - I don't know how to define the slices to apply to any axis
    counter = [0] * len(keys)
    for in_file in in_files:
        with h5py.File(in_file, "r") as h5in:
            for ii in range(num_datasets):
                ds = h5in[keys[ii]][()]
                indices = slice(counter[ii], counter[ii] + ds.shape[axis])
                with h5py.File(out_file, "r+") as h5out:
                    h5out[keys[ii]][indices, ...] = ds
                counter[ii] += ds.shape[axis]
