from scipy.io import readsav
from astropy.io import fits
from astropy.time import Time, TimeDelta
import astropy.units as u
import numpy
from math import ceil

class LazyLoader:
    def __init__(self, data_file, loader_function):
        self.data_file = data_file
        self._data = None  # Placeholder for the loaded data
        self.loader_function = loader_function

    @property
    def data(self):
        if self._data is None:
            # Load the data using the provided loader function when it's first requested
            self._data = self.loader_function(self.data_file)
        return self._data
    
def load_rfistat_data(data_file):
    # Load data from a .sav file
    data = readsav(data_file)
    return data


def load_data_from_fits(data_file, stokes = 'V/I'):
    with fits.open(data_file, memmap=True) as hdus:
        len_fits_file_ = len(hdus)
        if len_fits_file_ > 6:
            variable = hdus[2].data
            nt_ = variable.NT[0]
            nf_ = variable.NF[0]
            ns_ = variable.NS[0]


        
            if ns_ == 4:
                datasize = 4 * nt_ * nf_ * ns_
            if ns_ == 3:
                datasize = 3 * nt_ * nf_ * ns_
            datasizemax = 2**31 - 1
            
            
            if stokes.lower() != 'rm':
                if datasize <= datasizemax:
                    data = numpy.array(hdus[3].data.T)
                    k = 0
                else:
                    data = numpy.zeros((ns_, nf_, nt_))
                    for k in range(ns_):
                        data[k, :, :] = numpy.array(hdus[3 + k].data)
                    data = data.T
                
                frequency = numpy.array((hdus[6 + k].data * u.MHz).value)   

                    
            else:
                    if datasize <= datasizemax:
                        k = 0
                    else:
                        k = 3
                    data = numpy.array(hdus[-2].data.T)
                    frequency = numpy.array(hdus[-1].data)
                    
            rfilevel0 = numpy.array(hdus[4 + k].data.T)
            #time = numpy.array((Time(hdus[2].data['timestamp'][0], format='unix') + TimeDelta(hdus[5 + k].data, format='sec')).value)         
            time_ = da.from_array(Time(hdus[2].data["timestamp"][0] + hdus[8].data, format="unix").unix, chunks = chunk_size_time)
        else:
            print,'Something went wrong with the fits file. Length does not correspond. Are you sure of you fits file?'

    return time, frequency, data, rfilevel0

def load_RFI_data_from_fits(data_file):
    rfi_mask = []
    with fits.open(data_file, memmap=True) as hdus:
        ndatasize = hdus[2].data
        ndatabit = numpy.array(hdus[3].data)

        # Convert bits array to bytes array
        ndatabyte = bitarray_to_bytearray(ndatabit, ndatasize)
        
        rfimask_level1 = ndatabyte[:,:, 0]
        rfimask_level2 = ndatabyte[:,:, 1]
        rfimask_level3 = ndatabyte[:,:, 2]


    return rfimask_level1, rfimask_level2, rfimask_level3


def multiply_data(data1, data2):
    """
    Multiply data1 by data2 element-wise.

    Parameters:
    - data1 (Dask array): First array to multiply.
    - data2 (Dask array): Second array to multiply.

    Returns:
    Dask array resulting from element-wise multiplication.
    """
    
    return data1 * data2

def safe_divide(numerator, denominator):
    result = numpy.empty_like(numerator, dtype=float)
    zero_mask = denominator == 0
    result[~zero_mask] = numerator[~zero_mask] / denominator[~zero_mask]
    return result


def bitarray_to_bytearray(xbit, xsize):
    """Converts a bits array to a bytes array.

    Parameters:
        xbit: int
            Input bitarray to be converted to a bytearray.
        xsize: list of int
            Size of the original bytes array dimensions.

    Returns:
        numpy array
            Bytes array obtained by converting the input bits array.

    Notes:
        The function performs a bitwise operation to convert the input bits array to a bytes array. The resulting bytes array
        is reshaped based on the size of the input original bytes array dimensions.
    """


def bitarray_to_bytearray(xbit, xsize):
    ns = len(xsize)
    nx = xsize[-1]
    nx8 = numpy.ceil(nx / 8).astype(numpy.int64) * 8
    xbyte = numpy.reshape(
                numpy.transpose(
                        numpy.reshape(
                                     [(xbit & 0b10000000) // 128, (xbit & 0b01000000) // 64,
                                     (xbit & 0b00100000) // 32, (xbit & 0b00010000) // 16,
                                     (xbit & 0b00001000) // 8, (xbit & 0b00000100) // 4,
                                     (xbit & 0b00000010) // 2, (xbit & 0b00000001) // 1],
                                     (nx8 // 8, 8), order='F'
                                     )
                                ),
                            nx8)
    ndim = xsize[0]
    if ndim == 1:
        xbyte = numpy.reshape(xbyte[:nx], xsize[1], order='F')
    elif ndim == 2:
        xbyte = numpy.reshape(xbyte[:nx], (xsize[1], xsize[2]), order='F')
    elif ndim == 3:
        xbyte = numpy.reshape(xbyte, (xsize[1], xsize[2], xsize[3]), order='F')
    elif ndim == 4:
        xbyte = numpy.reshape(xbyte[:nx], (xsize[1], xsize[2], xsize[3], xsize[4]), order='F')


    return xbyte
