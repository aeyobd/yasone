# Yasone data processing

These notes describe my data processing steps to take the raw GTC/OSIRIS+ data and create reliable and deep photometric catalogues of stars for reliable detection and characterization of possible member stars for each object.  

# Overview

## About the instrument and observations

GTC/OSIRIS+ is a multi-purpose spectrograph and imager. 

We have observations of about 15 exposures in each colour (g, r, and i) for each object. For every observation night, we also have associated calibration frames.



## Initial data reduction

We initially stack bias images, then create flats for each night and create unbiased, flat-fielded and trimmed images for further analysis. 

The noise on the bias images appears to be around 2 ADU/frame (likely near the read noise of the detector). In addition, the flat-field noise is approximately consistent with poisson statistics given the flat-field exposure time.

Flats show increased vignetting near edges, where optical efficiency may drop below 70%. However, we will likely discard the edges regardless in the combined image. 

We calculate the pixel-by-pixel RMS as follows
$$
\delta^2 \textrm{pixel} = \delta^2 \textrm{background}  + ...
$$

## Removing artefacts and bright stars

We mask out pixels known to be bad (from SAUSERO).

We adopt a saturation level of 62,000, as is used in the XXX catalogue. 

Bright stars cause large saturation wings which spread out in multiple directions. To identify such features, we use a combination of binary dilation of saturated pixels, and an algorithm to calculate the convex hull of saturated regions and fill the inner region with pixels. 

To stack images, we use SCAMP to generate warping and alignment corrections and SWARP to create de-warped stacks of each image. Residuals in *Gaia* appear to be well-behaved and small. We use median combining (and median combine the weights from swarp). 

We specially treat the brightest stars affected by large-scale saturation extending more than *16* pixels. After masking the saturated regions, we fit Moffat profiles (with powerlaw index *1.1*) to each star and subtract the model from the image, one at a time starting from the stars with the most saturated pixels. This step reduces the "halo" around the bright stars, however regions within about 30 pixels of the brightest stars retain strong residuals. We *will* then mask these affected regions using a circular aperture. Finally, we automatically mask out connected regiones 3$\sigma$ above the background after subtracting the psf model from each bright star to flag regions affected by diffraction and saturation spikes. 


## Photometric measurements

For source detection, we initially fit the bright-star-subtracted image using source extractor. We then subtract out all model images (psf or galaxy?) and rerun source extractor to find fainter sources. From these detections, we build a catalogue of possible sources from the stacked images. 

We *will* use both aperture and PSF photometry. An empirical PSF is generated with PSFEx which is then used with source extractor. 



## Star galaxy seperation

To identify likely stars, we use:
$$
0.1< I2 - 0.5 I3 < 0.5
$$
in the $r$ band where $r_1, r_2, r_3$ are the xx, xx, and xx apertures respectively.

## Photometric calibration

Using the same photometric and data reduction pipeline as for science images, we measure the photometry of standard stars as provided by OSIRIS. One complication is that the standard stars are unfortunantly saturated.

Finally, we cross match the catalogue against PANSTARRS for each science field and derive the zero-point correction based on the median offset of unsaturated sources near

We additionally fit the magnitude deviation of  saturated stars with a quadratic polynomial, enabling reliable recovering of the photometry of mildy-saturated stars (less than 16 saturated pixels) up to magnitudes 16??!.



# Catalogue & Data processing details

Steps:

Imaging:

1. trim_images.py, trims images based on fits header, creating images in `obsname/trimmed`
2. stack_biases.py creates master biases for each observation date, creating a `stacked_bias.fits` in each calibration date folder
3. subtract_biases.py linearly subtracts the bias from each trimmed frame, creating imagings in `obsname/unbiased`
4. `stack_flats.py` stacks flats in each colour.
5. `flat_field.py` creates flat-fielded images in the directories `obsname/flat_fielded`
6. `make_flags.py` creates the flag maps for flags 4, 8, and 16 (see below).
7. `make_errs.py` calculates the total pixel-by-pixel RMS error for the flat-fielded images.

Photometry:

1. `link_files.py` simply reorganizes the file trees for photometric analysis (now one reduced image per directory)
2. `astrometrize_all.sh` astrometrizes images using astrometry.net to create an initial calibration that SCAMP will improve on
3. `run_se_all.sh` runs an initial high-sigma detection on each image for WCS calibration. This is also where we background-subtract images
4. `run_scamp.sh` Runs scamp for each image group on each folder and outputs updated `nobkg.head` files
5. `make_bkg_err.py` Adds in quadrature the background uncertainty
6. `run_swarp.py` and `run_swarp_flags.py` runs swarp on the images, uncertainties, and flag files, created stacked images in directories like `yasone2/stacked_i` 
7. `analyze_bright_stars.jl` creates a catalogue of very bright stars on a given image, subtracts out a model psf, and identifies pixels affected by flags 1 and 2. 
8. `run_se_stack_detection.sh` runs an initial detection on stacked images
9. `run_iterative_photometry.sh` iteratively detects and creates residual images on stacks
10. `combine_catalogues.py` combines the detection i, g, r catalogues on each source to create a main list of possible sources with identifiers
11. `run_forced_photometry.sh` now uses the combined catalogues to run forced photometry on each source.
12. `calibrate_vs_panstarrs.py` calculates zeropoints based on the panstarrs and notes our calibration results and calibration catalogues
13. `combine_all_photometry.py` combines the photometry for each band into a final uncalibrated catalogue. 

Analysis (TODO):

1. make_cmd_detection_plot.py
2. fit_cmd.py
3. fit_light_profile.py



flag images:

- 1: Possible diffraction spike
- 2: Near bright star
- 4: Individually saturated pixel
- 8: Part of very saturated region near bright star
- 16: bad pixel



Catalogue columns

- 