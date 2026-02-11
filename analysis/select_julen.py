from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy

import os
import sep
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.interpolate import interp1d

from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse


def get_obs_dates(obsname):
    obs_i = int(obsname)
    if obs_i == 1:
        date_g = date_r = date_i = "20230708"
    elif obs_i == 2:
        date_g = date_r = "20230708"
        date_i = "20230709"
    elif obs_i == 3:
        date_g = date_r = date_i = "20230709"
    elif obs_i == 4:
        date_g = date_r = date_i = "20230709"
    elif obs_i == 5:
        date_g = date_r = "20230709"
        date_i = "20230710"
    elif obs_i == 6:
        date_g = date_r = date_i = "20230710"
    elif obs_i == 7:
        date_g = date_r = date_i = "20230724"
    elif obs_i in [8, 9]:
        date_g = date_r = date_i = "20230816"
    else:
        raise Exception(f"observation not known {obsname}")

    return date_g, date_r, date_i


def select_julen_sample(objname, obsname="0001", sigma=6, aperture_radius=1.5 *
                        1.7, aperture_radius_g=None, aperture_radius_r=None,
                        aperture_radius_i=None, fixed_aperture=False,
        error_scale=1, max_xmatch_radius=1
    ):
    ###################################
    # Substitute with your own directory
    ###################################
    # daniel: updated to be more flexible
    if obsname == "stack":
        fits_file_g = fits.open(f'./julen_stacks/{objname}_g_stack.fits')
        fits_file_r = fits.open(f'./julen_stacks/{objname}_r_stack.fits')
        fits_file_i = fits.open(f'./julen_stacks/{objname}_i_stack.fits')
    else:
        date_g, date_r, date_i = get_obs_dates(obsname)

        fits_file_g = fits.open(f'../imaging/julen_stacked//GTC73-23A_{obsname}_Sloan_g_{date_g}_NOSKY_BBI.fits')
        fits_file_r = fits.open(f'../imaging/julen_stacked//GTC73-23A_{obsname}_Sloan_r_{date_r}_NOSKY_BBI.fits')
        fits_file_i = fits.open(f'../imaging/julen_stacked//GTC73-23A_{obsname}_Sloan_i_{date_i}_NOSKY_BBI.fits')

    # daniel: flexible aperture radius
    if aperture_radius_g is None:
        aperture_radius_g = aperture_radius
    if aperture_radius_r is None:
        aperture_radius_r = aperture_radius
    if aperture_radius_i is None:
        aperture_radius_i = aperture_radius


    image_data_g = fits_file_g[0].data
    header_g = fits_file_g[0].header

    image_data_r = fits_file_r[0].data
    header_r = fits_file_r[0].header

    image_data_i = fits_file_i[0].data
    header_i = fits_file_i[0].header

    ###################################
    # I got an error sometimes when 
    # trying to display the data.
    # This should solve it
    ###################################

    # Daniel: updated from Julen to account for changes to numpy
    if image_data_g.dtype.byteorder not in ('=', '|'):
        image_data_g = image_data_g.astype(image_data_g.dtype.newbyteorder("="))
        
    if image_data_r.dtype.byteorder not in ('=', '|'):
        image_data_r = image_data_r.astype(image_data_r.dtype.newbyteorder("="))
        
    if image_data_i.dtype.byteorder not in ('=', '|'):
        image_data_i = image_data_i.astype(image_data_i.dtype.newbyteorder("="))


    #######################################
    # Here I subtract the background noise
    #######################################

    image_data_contiguous_g = np.ascontiguousarray(image_data_g)
    image_data_contiguous_r = np.ascontiguousarray(image_data_r)
    image_data_contiguous_i = np.ascontiguousarray(image_data_i)

    bkg_g = sep.Background(image_data_contiguous_g)
    bkg_r = sep.Background(image_data_contiguous_r)
    bkg_i = sep.Background(image_data_contiguous_i)

    image_data_bkg_subtracted_g = image_data_contiguous_g - bkg_g
    image_data_bkg_subtracted_r = image_data_contiguous_r - bkg_r
    image_data_bkg_subtracted_i = image_data_contiguous_i - bkg_i

    # daniel: updated to use `sigma  instead of hard-coded value (6)
    objects_g = sep.extract(image_data_bkg_subtracted_g, sigma, err = bkg_g.globalrms)
    objects_r = sep.extract(image_data_bkg_subtracted_r, sigma, err = bkg_r.globalrms)
    objects_i = sep.extract(image_data_bkg_subtracted_i, sigma, err = bkg_i.globalrms)

    df_g = pd.DataFrame(objects_g)
    df_r = pd.DataFrame(objects_r)

    # daniel: new ID column
    df_r["id"] = np.arange(len(df_r))
    df_i = pd.DataFrame(objects_i)


    df_g['fwhm'] = 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt((df_g['a']**2 + df_g['b']**2) / 2.0)
    df_r['fwhm'] = 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt((df_r['a']**2 + df_r['b']**2) / 2.0)
    df_i['fwhm'] = 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt((df_i['a']**2 + df_i['b']**2) / 2.0)

    df_g['ellipticity'] = 1 - (df_g['b'] / df_g['a'])
    df_r['ellipticity'] = 1 - (df_r['b'] / df_r['a'])
    df_i['ellipticity'] = 1 - (df_i['b'] / df_i['a'])

    # daniel: instead of 1.7*df_g['fwhm'], now 1.7 is variable and can use
    # fixed aperture instead.

    if fixed_aperture:
        r = aperture_radius_g
    else:
        r = aperture_radius_g * df_g['fwhm']

    flux_g, fluxerr_g, flag_g = sep.sum_circle(image_data_bkg_subtracted_g, objects_g['x'], objects_g['y'], 
                                               r, err = bkg_g.globalrms, gain = 1.0)

    if fixed_aperture:
        r = aperture_radius_r
    else:
        r = aperture_radius_r * df_r['fwhm']
    flux_r, fluxerr_r, flag_r = sep.sum_circle(image_data_bkg_subtracted_r, objects_r['x'], objects_r['y'], 
                                               r, err = bkg_r.globalrms, gain = 1.0)

    if fixed_aperture:
        r = aperture_radius_i
    else:
        r = aperture_radius_i * df_i['fwhm']
    flux_i, fluxerr_i, flag_i = sep.sum_circle(image_data_bkg_subtracted_i, objects_i['x'], objects_i['y'], 
                                               r, err = bkg_i.globalrms, gain = 1.0)


    # daniel: correct stack magnitudes
    if obsname == "stack":
        flux_g /= 3
        flux_r /= 3
        flux_i /= 3
        fluxerr_g /= 3
        fluxerr_r /= 3
        fluxerr_i /= 3


    extinction_coefficient_g = 0.15 
    extinction_coefficient_r = 0.07 
    extinction_coefficient_i = 0.04 

    extinction_coefficient_error_g = 0.02
    extinction_coefficient_error_r = 0.01
    extinction_coefficient_error_i = 0.01

    # daniel: updating to just use the fits file value
    #airmass_array_g = [1.095, 1.090, 1.085, 1.080, 1.075]
    #airmass_array_r = [1.071, 1.067, 1.063, 1.060, 1.057]
    #airmass_array_i = [1.054, 1.051, 1.049, 1.046, 1.045]

    #airmass_g = np.mean(airmass_array_g)
    #airmass_r = np.mean(airmass_array_r)
    #airmass_i = np.mean(airmass_array_i)
    airmass_g = header_g["AIRMASS"]
    airmass_r = header_r["AIRMASS"]
    airmass_i = header_i["AIRMASS"]

    #airmass_error_g = np.std(airmass_array_g)
    #airmass_error_r = np.std(airmass_array_r)
    #airmass_error_i = np.std(airmass_array_i)
    # TODO: fix this arbitrarily set airmass uncert
    airmass_error_g = 0.05
    airmass_error_r = 0.05
    airmass_error_i = 0.05

    # daniel: the below was commented out, now includes 
    # conditionals and variable ra0, dec0...

    if objname == "yasone1":
        zero_point_estimation_g = 28.15227
        zero_point_estimation_error_g = 0.01692 
        
        zero_point_estimation_r = 28.24766 
        zero_point_estimation_error_r = 0.01008
        
        zero_point_estimation_i = 27.81407
        zero_point_estimation_error_i = 0.01047
        ra0, dec0 = 265.52019, 13.17146


    elif objname == "yasone2":
        zero_point_estimation_g = 28.12761 
        zero_point_estimation_error_g = 0.01931
        
        zero_point_estimation_r = 28.13342
        zero_point_estimation_error_r = 0.00882
        
        zero_point_estimation_i = 27.61163
        zero_point_estimation_error_i = 0.01528
        ra0, dec0 = 262.34921, 6.42302

    elif objname == "yasone3":
        
        zero_point_estimation_g = 27.93537
        zero_point_estimation_error_g = 0.01880
        
        zero_point_estimation_r = 28.06096 
        zero_point_estimation_error_r = 0.01128
        
        zero_point_estimation_i = 27.49029
        zero_point_estimation_error_i = 0.01095
        ra0, dec0 = 292.96258, -26.44994



    ###################################
    # Here is where we finally obtain flux 
    # and magnitude values
    ###################################

    df_g['mag_inst'] = -2.5 * np.log10(flux_g)
    df_r['mag_inst'] = -2.5 * np.log10(flux_r)
    df_i['mag_inst'] = -2.5 * np.log10(flux_i)

    df_g['mag_inst_err'] = (2.5 / np.log(10)) * (fluxerr_g / flux_g)
    df_r['mag_inst_err'] = (2.5 / np.log(10)) * (fluxerr_r / flux_r)
    df_i['mag_inst_err'] = (2.5 / np.log(10)) * (fluxerr_i / flux_i)

    df_g['mag_apparent'] = zero_point_estimation_g + df_g['mag_inst'] - (extinction_coefficient_g * airmass_g)
    df_r['mag_apparent'] = zero_point_estimation_r + df_r['mag_inst'] - (extinction_coefficient_r * airmass_r)
    df_i['mag_apparent'] = zero_point_estimation_i + df_i['mag_inst'] - (extinction_coefficient_i * airmass_i)

    df_g['mag_apparent_error'] = np.sqrt(zero_point_estimation_error_g**2 + df_g['mag_inst_err']**2 +
                                        (airmass_g * extinction_coefficient_error_g)**2 + 
                                        (extinction_coefficient_g * airmass_error_g)**2)

    df_r['mag_apparent_error'] = np.sqrt(zero_point_estimation_error_r**2 + df_r['mag_inst_err']**2 +
                                        (airmass_r * extinction_coefficient_error_r)**2 + 
                                        (extinction_coefficient_r * airmass_error_r)**2)

    df_i['mag_apparent_error'] = np.sqrt(zero_point_estimation_error_i**2 + df_i['mag_inst_err']**2 +
                                        (airmass_i * extinction_coefficient_error_i)**2 + 
                                        (extinction_coefficient_i * airmass_error_i)**2)


    # Get the World Coordinate System information

    wcs_g = WCS(header_g)
    wcs_r = WCS(header_r)
    wcs_i = WCS(header_i)

    # Convert pixel coordinates to world coordinates (RA, DEC)
    # The second argument is the origin (1 for FITS convention)

    ra_dec_g = wcs_g.wcs_pix2world(np.column_stack((df_g['x'], df_g['y'])), 1)
    df_g['ra'] = ra_dec_g[:, 0]
    df_g['dec'] = ra_dec_g[:, 1]

    ra_dec_r = wcs_r.wcs_pix2world(np.column_stack((df_r['x'], df_r['y'])), 1)
    df_r['ra'] = ra_dec_r[:, 0]
    df_r['dec'] = ra_dec_r[:, 1]

    ra_dec_i = wcs_i.wcs_pix2world(np.column_stack((df_i['x'], df_i['y'])), 1)
    df_i['ra'] = ra_dec_i[:, 0]
    df_i['dec'] = ra_dec_i[:, 1]



    ###################################
    # Here I just did a little zoom to focus 
    # on the region around our candidate's centroid. 
    # This step can of course take place much 
    # earlier or later, to your liking.
    # For Yasone 2 and 3, their centroids are in the paper.
    # These centroids were determined when I was searching for overdensities
    # in Pan-STARRS, and they just happened to be the coordinates of the stellar
    # member with the highest number of stellar neighbours, but of course
    # this centroid can be adjusted. I leave it to your criteria. 
    ###################################

    # daniel: changed hard-coded coordinates to ra0, dec0
    df_cut_g = df_g[(df_g['ra'] >= ra0 - 0.5) & (df_g['ra'] <= ra0 + 0.5) &
                    (df_g['dec'] >= dec0 - 0.5) & (df_g['dec'] <= dec0 + 0.5)].copy()

    df_cut_r = df_r[(df_r['ra'] >= ra0 - 0.5) & (df_r['ra'] <= ra0 + 0.5) &
                    (df_r['dec'] >= dec0 - 0.5) & (df_r['dec'] <= dec0 + 0.5)].copy()

    df_cut_i = df_i[(df_i['ra'] >= ra0 - 0.5) & (df_i['ra'] <= ra0 + 0.5) &
                    (df_i['dec'] >= dec0 - 0.5) & (df_i['dec'] <= dec0 + 0.5)].copy()


    # p1_cut = p1_cut.reset_index(drop = True)
    df_cut_g = df_cut_g.reset_index(drop = True)
    df_cut_r = df_cut_r.reset_index(drop = True)
    df_cut_i = df_cut_i.reset_index(drop = True)



    ###################################
    # Now that we have 3 different dataframes with 
    # a varying number of sources 
    # (depending on the threshold value, the definition of the aperture etc.)
    # this is how I would typically "join them" so that we can evaluate all 
    # the bands information for each source. Since I was especially interested
    # in recreating g vs r-g and r vs r-i plots, I created 2 dataframes, one containing
    # the information from the g and r bands, and another one with r and i.
    ###################################

    observed_coords_g = df_cut_g[['ra', 'dec']].to_numpy()
    observed_coords_r = df_cut_r[['ra', 'dec']].to_numpy()
    observed_coords_i = df_cut_i[['ra', 'dec']].to_numpy()

    # Create KDTree for the bottleneck (the array with least amount of sources). In g-r that's g. In r-i that's r.
    tree_g = KDTree(observed_coords_g)
    tree_r = KDTree(observed_coords_r)

    # Query the KDTree for the nearest neighbor of each observed coordinate. Here is where you specify the secoundary array (r, i)
    # daniel, remove distance upper bound
    distances_gr, indices_gr = tree_g.query(observed_coords_r,)  # 1 arcsecond = 1/3600 degree
    distances_ri, indices_ri = tree_r.query(observed_coords_i,)

    # Filter out matches within 1 arcsecond
    valid_matches_gr = distances_gr < np.inf
    valid_indices_gr = indices_gr[valid_matches_gr]
    valid_distances_gr = distances_gr[valid_matches_gr]

    valid_matches_ri = distances_ri < np.inf
    valid_indices_ri = indices_ri[valid_matches_ri]
    valid_distances_ri = distances_ri[valid_matches_ri]

    # Build the dataframe of matched sources for g-r and r-i
    # ¡¡¡IMPORTANT!!!
    # WHEN EXTRACTING THE ROWS OF INTEREST FROM THE ORIGINAL BOTTLENECK DATAFRAME ==> USE valid_indices
    # WHEN EXTRACTING THE ROWS OF INTEREST FROM THE ORIGINAL LARGER DATAFRAME ==> USE valid_matches
    # daniel: added a few new columns, including distances
    matched_sources_gr = pd.DataFrame({
        'ra_g': df_cut_g.iloc[valid_indices_gr]['ra'].values,
        'dec_g': df_cut_g.iloc[valid_indices_gr]['dec'].values,
        'mag_g': df_cut_g.iloc[valid_indices_gr]['mag_apparent'].values,
        'mag_g_err': df_cut_g.iloc[valid_indices_gr]['mag_apparent_error'].values,
        'ellipticity_g': df_cut_g.iloc[valid_indices_gr]['ellipticity'].values,
        'ra_r': df_cut_r.iloc[valid_matches_gr]['ra'].values,
        'dec_r': df_cut_r.iloc[valid_matches_gr]['dec'].values,
        'mag_r': df_cut_r.iloc[valid_matches_gr]['mag_apparent'].values,
        'mag_r_err': df_cut_r.iloc[valid_matches_gr]['mag_apparent_error'].values,
        'ellipticity_r': df_cut_r.iloc[valid_matches_gr]['ellipticity'].values,
        'distance_gr': valid_distances_gr,
        'id': df_cut_r.iloc[valid_matches_gr]['id'].values,
        'fwhm': df_cut_r.iloc[valid_matches_gr]['fwhm'].values,

    })

    matched_sources_gr['mag_gr'] = matched_sources_gr['mag_g'] - matched_sources_gr['mag_r']
    matched_sources_gr['mag_gr_err'] = np.sqrt(matched_sources_gr['mag_g_err']**2 + matched_sources_gr['mag_r_err']**2)

    matched_sources_ri = pd.DataFrame({
        'ra_r': df_cut_r.iloc[valid_indices_ri]['ra'].values,
        'dec_r': df_cut_r.iloc[valid_indices_ri]['dec'].values,
        'mag_r': df_cut_r.iloc[valid_indices_ri]['mag_apparent'].values,
        'mag_r_err': df_cut_r.iloc[valid_indices_ri]['mag_apparent_error'].values,
        'ellipticity_r': df_cut_r.iloc[valid_indices_ri]['ellipticity'].values,
        'ra_i': df_cut_i.iloc[valid_matches_ri]['ra'].values,
        'dec_i': df_cut_i.iloc[valid_matches_ri]['dec'].values,
        'mag_i': df_cut_i.iloc[valid_matches_ri]['mag_apparent'].values,
        'mag_i_err': df_cut_i.iloc[valid_matches_ri]['mag_apparent_error'].values,
        'ellipticity_i': df_cut_i.iloc[valid_matches_ri]['ellipticity'].values,
        'distance_ri': valid_distances_ri,
        'id': df_cut_r.iloc[valid_indices_ri]['id'].values,

    })

    matched_sources_ri['mag_ri'] = matched_sources_ri['mag_r'] - matched_sources_ri['mag_i']
    matched_sources_ri['mag_ri_err'] = np.sqrt(matched_sources_ri['mag_r_err']**2 + matched_sources_ri['mag_i_err']**2)



    # Selecting sources from the candidate g-r
    # daniel: use ra0, dec0 instead of hard-coded coordinates
    radius_g = np.sqrt((matched_sources_gr['ra_g'] - ra0)**2 + (matched_sources_gr['dec_g'] - dec0)**2)
    radius_r = np.sqrt((matched_sources_gr['ra_r'] - ra0)**2 + (matched_sources_gr['dec_r'] - dec0)**2)
    radius_gr = (radius_g + radius_r)/2

    matched_sources_gr['radius'] = radius_gr


    #############################
    # daniel: My own analysis
    #############################
    matched_sources_both = matched_sources_gr.join(matched_sources_ri.set_index("id")[["ra_i", "dec_i", "mag_i", "mag_i_err", "distance_ri", "mag_ri", "mag_ri_err"]], on="id", )
    matched_sources_both

    matched_sources_both.rename(
       columns= {
        "mag_r": "RMAG",
        "mag_g": "GMAG",
        "mag_i": "IMAG",
        "ra_r": "ALPHA_J2000",
        "dec_r": "DELTA_J2000",
        
        
        }
    )


    # daniel: use variable xmatch radius. Note that the only columns affected
    # are mag_g and mag_i. Use these magnitudes for ...

    matched_sources_both.loc[matched_sources_both["distance_gr"] > max_xmatch_radius/3600, "mag_g"] = np.nan
    matched_sources_both.loc[matched_sources_both["distance_ri"] >
                             max_xmatch_radius/3600, "mag_i"] = np.nan


    # corrections to go to panstarrs magnitude
    matched_sources_both['mag_g_ps'] = (matched_sources_both['mag_g'] - 0.011 - 0.125*(matched_sources_both['mag_gr']) 
                                           - 0.015 * (matched_sources_both['mag_gr']** 2))

    matched_sources_both['mag_r_ps'] = (matched_sources_both['mag_r'] + 0.001 - 0.006*(matched_sources_both['mag_gr'])
                                           - 0.002 * (matched_sources_both['mag_gr'] ** 2))

    matched_sources_both['mag_i_ps'] = (matched_sources_both['mag_i'] + 0.004 - 0.014*(matched_sources_both['mag_gr'])
                                           - 0.001 * (matched_sources_both['mag_gr'] ** 2))

    # Error propagation for mag_g_good

    matched_sources_both['mag_g_ps_err'] = np.sqrt(
        ((0.875 - 0.03 * matched_sources_both['mag_gr']) * matched_sources_both['mag_g_err']) ** 2 +
        ((-0.125 + 0.03 * matched_sources_both['mag_gr']) * matched_sources_both['mag_r_err']) ** 2 +
        (0.006**2) #Associated error of the conversion, Table 6 THE Pan-STARRS1 PHOTOMETRIC SYSTEM Tonry et al. 2012
    )

    # Error propagation for mag_r_good

    matched_sources_both['mag_r_ps_err'] = np.sqrt(
        ((-0.006 - 0.004 * matched_sources_both['mag_gr']) * matched_sources_both['mag_g_err']) ** 2 +
        ((1.006 - 0.004 * matched_sources_both['mag_gr']) * matched_sources_both['mag_r_err']) ** 2 +
        (0.002**2) #Associated error of the conversion, Table 6 THE Pan-STARRS1 PHOTOMETRIC SYSTEM Tonry et al. 2012
    )

    matched_sources_both['mag_i_ps_err'] = np.sqrt(
        ((0.004 - 0.024 * matched_sources_both['mag_gr']) * matched_sources_both['mag_g_err']) ** 2 +
        ((0.004 - 0.024 * matched_sources_both['mag_gr']) * matched_sources_both['mag_r_err']) ** 2 +
        + matched_sources_both['mag_i_err']**2 + 
        (0.003**2) #Associated error of the conversion, Table 6 THE Pan-STARRS1 PHOTOMETRIC SYSTEM Tonry et al. 2012
    )

    matched_sources_both['mag_gr_ps'] = matched_sources_both['mag_g_ps'] - matched_sources_both['mag_r_ps']

    matched_sources_both['mag_gr_ps_err'] = np.sqrt(
        matched_sources_both['mag_g_ps_err']**2 + matched_sources_both['mag_r_ps_err']**2)

    matched_sources_both['mag_ri_ps'] = matched_sources_both['mag_r_ps'] - matched_sources_both['mag_i_ps']

    matched_sources_both['mag_ri_ps_err'] = np.sqrt(
        matched_sources_both['mag_r_ps_err']**2 + matched_sources_both['mag_i_ps_err']**2)



    matched_out = matched_sources_both.rename(
        columns={
        "mag_r": "RMAG",
        "mag_g": "GMAG",
        "mag_i": "IMAG",
        "mag_r_err": "RMAG_ERR",
        "mag_g_err": "GMAG_ERR",
        "mag_i_err": "IMAG_ERR",
        "ra_r": "ALPHA_J2000",
        "dec_r": "DELTA_J2000",
        
        
        }
    )


    tab_out = astropy.table.Table.from_pandas(matched_out)

    return tab_out
