# Questions

- mean vs median, sigma-clipping outlier removal does not fully address CR biases...
- dealing with (satellite) contamination, bright CR events
- level of data inspection?
- when to stack science images?

# General notes

- bad pixel mask within trim region appears to only remove top rows of images, which appear to do weird things nontheless and the upper corners are heavily vignetted
- 

# Object 1

This object is one of three separate exposures targeting Yasone 1. 



### Bias fiels

- No comments, everything looks okay. Total bias variation is < 5 ADU across detector!
- General notes:
- I assume horizontal banding is read noise?
- Some frindge like behaviour near column 0, but this is well into the overscan region
- There is a sharp change in bias at the right edge but this should region will likely be trimmed during image combination. This begins around absolute pixel number 2030
- No patterns in residuals besides banding and absolute normalization shifts over time.

Flat fields

- No major comments
- Large obscured regions (due to filter wheel) is mostly removed when restricting to TRIMSEC
- some moderate vignetting (30%) near the upper left and right edge are apparent
- A few dust particles, and diagonal bands appear on the G filter
- In the residuals, there is field-level variations. Most prominantly, there are strange ripple-like large scale patterns on the right edge of the detector. However, nearly all variations are less than 0.2% relative to the combined flat-field
- The i-filter flats have a strong (cosmic-ray) vertical group of pixels which is of reduced value in most frames but enhanced counts in one frame
- A few possibly hot pixels (mid-upper left)

Standard star fields:

- Read noise appears to be high
- ereduced images show some vertical bands of pixels which are darker
- Otherwise looks good
- no concerning residuals

Object fields:

- No apparent problems besides the saturation wings of a few bright stars, observations appear to be of excellent quality
- Parachutes (artefacts) appear in middle lower corner: watch for spurious point sources and arcs in reduced image
- smallscale changes in background apparent at top right corner (like in flat fields)
- diffraction spikes slightly different between frames (likely good!)
- Some shifts to overall background (~30 ADU), but up to 200 ADU shift for r-2047
- Strange bands (cosmic rays? or satellites) in top left of i-band images. Appears to move across image over time, extending out to pixel coordinate~500-1500



# Object 2

Bias:

- same comments as above

Flat:

- same comments as above
- 990 appears to have a cosmic ray

Standard:

- same as above

Object:

- Moving artefact begining left of bright star and apparent in all images in g, r, and first i image, probably needs to be treated seperately
- same as above otherwise (similar parashute)
- I-band images suffer strong shifts in background magnitude, likely need to sky-subtract or renormalize before stacking images

# Object X

Bias:

- 

Flat:

- 

Standard:

- 

Object:

- 
