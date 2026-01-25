import numpy as np
from astropy.table import Table
import warnings


class ISOCMD:
    """
    A class representing an isochrone
    """
    
    def __init__(self, filename, verbose=True):
        """
        Initialize ISOCMD from a file
        
        Parameters
        ----------
        filename : str
            The name of the file
        verbose : bool, optional
            Whether to print status messages (default: True)
        """
        if verbose:
            print(f"Reading in: {filename}")
        
        (self.version, self.photo_sys, self.abun, self.Av_extinction, 
         self.rot, self.log_ages, self.num_ages, self.hdr_list, 
         self.isocmds) = read_isocmd_file(filename)
        
        self.filename = filename
    
    def __getitem__(self, log_age):
        """
        Get the isochrone closest to the requested log age
        
        Parameters
        ----------
        log_age : float
            The log age in years
        
        Returns
        -------
        astropy.table.Table
            The isochrone table for the closest age
        """
        log_age = float(log_age)
        idx = age_index(self, log_age)
        return self.isocmds[idx]


def read_isocmd_header(filename):
    """
    Reads the header attributes of the isochrone file
    
    Parameters
    ----------
    filename : str
        Path to the isochrone file
    
    Returns
    -------
    tuple
        version, photo_sys, abun, Av_extinction, rot, num_ages
    """
    with open(filename, 'r') as f:
        lines = [line.split() for line in f.readlines()[:10]]
    
    version = {"MIST": lines[0][-1], "MESA": lines[1][-1]}
    photo_sys = " ".join(lines[2][4:])
    abun = {lines[4][i]: float(lines[5][i]) for i in range(1, 5)}
    rot = float(lines[5][-1])
    num_ages = int(lines[7][-1])
    Av_extinction = float(lines[8][-1])
    
    return version, photo_sys, abun, Av_extinction, rot, num_ages


def read_isocmd_file(filename):
    """
    Read the entire isochrone file
    
    Parameters
    ----------
    filename : str
        Path to the isochrone file
    
    Returns
    -------
    tuple
        version, photo_sys, abun, Av_extinction, rot, log_ages, 
        num_ages, hdr_list, isocmd_set
    """
    version, photo_sys, abun, Av_extinction, rot, num_ages = read_isocmd_header(filename)
    
    with open(filename, 'r') as f:
        data = [line.split() for line in f.readlines()[10:]]
    
    isocmd_set = []
    log_ages = []
    counter = 0
    hdr_list = data[counter + 2][1:]
    
    for i_age in range(num_ages):
        num_eeps = int(data[counter][-2])
        num_cols = int(data[counter][-1])
        hdr_list = data[counter + 2][1:]
        
        isocmd = np.zeros((num_eeps, num_cols))
        for eep in range(num_eeps):
            iso_cmd_chunk = data[counter + 3 + eep]
            iso_cmd_chunk = [float(x) for x in iso_cmd_chunk]
            isocmd[eep, :] = iso_cmd_chunk

        
        # Create Astropy Table
        table = Table(isocmd, names=hdr_list)
        
        isocmd_set.append(table)
        log_ages.append(table['log10_isochrone_age_yr'][0])
        counter += 3 + num_eeps + 2
    
    log_ages = np.array(log_ages)
    
    return version, photo_sys, abun, Av_extinction, rot, log_ages, num_ages, hdr_list, isocmd_set


def age_index(isocmd, log_age):
    """
    Returns the index of the isochrone closest to the requested age
    
    Parameters
    ----------
    isocmd : ISOCMD
        The ISOCMD object
    log_age : float
        The requested log age
    
    Returns
    -------
    int
        The index of the closest age
    
    Raises
    ------
    ValueError
        If the requested age is outside the available range
    """
    diff_arr = np.abs(isocmd.log_ages - log_age)
    idx = np.argmin(diff_arr)
    
    if log_age > np.max(isocmd.log_ages) or log_age < np.min(isocmd.log_ages):
        raise ValueError(
            f"The requested age {log_age} is outside the range. "
            f"Try log age {np.min(isocmd.log_ages)} and {np.max(isocmd.log_ages)}"
        )
    
    log_age_actual = isocmd.log_ages[idx]
    
    if diff_arr[idx] > 0.001:
        warnings.warn(
            f"The requested age is not in the isochrone set. "
            f"Rounding {log_age} to log_age = {log_age_actual}"
        )
    
    return idx

