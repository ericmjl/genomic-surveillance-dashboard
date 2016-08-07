import numpy as np

def iterable_shape(iterable):
    """
    Returns a numpy array version of the iterable, as well as its shape.
    """
    iterable = np.array(iterable)

    if len(iterable.shape) == 1:
        iterable = iterable.reshape(iterable.shape[0], -1)

    return iterable, iterable.shape


def yerrorbars(data, errors, xdata_label, ydata_label, sort=False):
    """
    Generates y-error bars properly formatted for input into Bokeh.


    """
    data_len = len(data[xdata_label])
    xdata = np.array(data[xdata_label])
    ydata = np.array(data[ydata_label])
    xs = [[i+1]*2 for i in range(data_len)]

    # There are a few cases that we are defining as acceptable.

    # Case 1: Uniform error bars across all x-axis samples.
    if isinstance(errors, int) or isinstance(errors, float):
        ys = [[ydata[i] - errors,
               ydata[i] + errors]
              for i in range(data_len)]

    # Case 2: Error bars are different per sample and are provided as an
    # iterable.
    if hasattr(errors, "__iter__"):
        # Firstly, convert the errors to a numpy array. This allows easy
        # checking of shape.
        if isinstance(errors, list) or isinstance(errors, tuple):
            errors, err_shape = iterable_shape(errors)
            print(errors, err_shape)

        # Next, check to make sure that the errors provided make sense.
        assert err_shape[0] == data_len,\
            "errors provided are not of the same length as data."
        assert err_shape[1] in [1, 2],\
            "errors provided must be a single row or two rows."

        # Case 2a: the error bars are symmetric around the data point.
        if err_shape[1] == 1:
            ys = [[ydata[i] - errors[i],
                   ydata[i] + errors[i]]
                  for i in range(data_len)]
        # Case 2b: the error bars are non-symmetric around the data point.
        # In this case, the errors are provided in the shape of 2xN.
        elif err_shape[1] == 2:
            ys = [[ydata[i] - errors[i, 0],
                   ydata[i] + errors[i, 1]]
                  for i in range(data_len)]

    if sort:
        xs = [x for (label, x) in sorted(zip(data[xdata_label], xs))]
        ys = [y for (label, y) in sorted(zip(data[xdata_label], ys))]

    return xs, ys
