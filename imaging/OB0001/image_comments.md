# Questions

- mean vs median, sigma-clipping outlier removal does not fully address CR biases...
- dealing with (satellite) contamination, bright CR events
- level of data inspection?
- when to stack science images?
- multiple exposures across different observation blocks?

# General notes

- bad pixel mask within trim region appears to only remove top rows of images, which appear to do weird things nontheless and the upper corners are heavily vignetted
- median



Objects 1-3; 4-6; 8-9 share calibration images

Objects 1-3 are for Yasone 1 (17:42:04.730 13:08:05.676), all observed on 20230708

Objects 4-6 are for Yasone 2 (17:29:23.691 06:23:11.275), all observed on 20230709

Objects 7-9 are for Yasone 3 (19:31:50.897 -26:29:11.482)

- Obj 7 observed on 20230723
- Obj 8 and 9 on 20230816 except flats on 20230811

for all science fields, the exposure time is 190.0 seconds

for all flats, they are twilight flats observed around 8:40pm on the same day, with exposures of 0.8 to 3.0 seconds

biases are 0.0 second exposures taken in early evening.

most standard fields are towards 15:30:39.431 05:59:01.480, with 0.5s (or 1.0 for g)  exposures, with 2 r band exposures and 1 in g and i..

However, objects 7 and 8-9 which are  22:41:34.940 01:08:58.156, and 15:32:06.640 06:25:37.482 respectively and with variable numbers of exposures with different exposure lengths.

# Yasone 1 (Obj 1-3)



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

Object fields (1):

- No apparent problems besides the saturation wings of a few bright stars, observations appear to be of excellent quality
- Parachutes (artefacts) appear in middle lower corner: watch for spurious point sources and arcs in reduced image
- smallscale changes in background apparent at top right corner (like in flat fields)
- diffraction spikes slightly different between frames (likely good!)
- Some shifts to overall background (~30 ADU), but up to 200 ADU shift for r-2047
- Strange bands (cosmic rays? or satellites) in top left of i-band images. Appears to move across image over time, extending out to pixel coordinate~500-1500

Object (2):

- Moving artefact begining left of bright star and apparent in all images in g, r, and first i image, probably needs to be treated seperately
- same as above otherwise (similar parashute)
- I-band images suffer strong shifts in background magnitude, likely need to sky-subtract or renormalize before stacking images

Object (3):

- blackground level shifts by over 1000 adu between frames, normalization likely required
- some minor mysterious ghots in bottom right and mid-centre
- one frame has much poorer psf shape than the other in r, r band also shows more severe (several 100s ADU) fluctuations in the right edge



# Yasone 2 (Obj 4-6)

Bias:

- no comments, only a few noticable small cosmic rays

Flat:

- some stars ? in the flat frames?, 2 in g, several feeatures in i

Standard:

- good

Object (4):

- 2236 has significant blur
- background shift in i-band, faint satellite blur

Object (5):

- 2253 aligned poorly (seeing?)
- background shift in i
- g band excellent
- satellite trail in 2255 and 2256
- 2258 poor image

Object (6):

- 2266, 2274 poor image
- background shift (strong) in gri



# Yasone 3 

## (Obj 7)

Bias:

- Very little frame-frame variation

Flat:

- some birght horizontal in r-band and many in i-band
- stronger right-edge ripples

Standard:

- suspistious background variation on all frames, towards top right, possibly problematic flat fields?

Object:

- g band excellent
- r band decent, background shift concerns
- i band appears problematic, oversubtracting unfocused images?? and strong background shifts



## (Obj 8-9)

Bias:

- somme shift ~60 adu in overscan region between stds frames
- overall systematic shift in bias over time

Flat:

- look excellent besides possible stars in the frame (top edge)?

Standard:

- strange viegnetting and missing regions for r and g, only last four images look satesfactor
- nonetheless large background variations, viegnetting towards top left corner and suffering ripples and variations on right edge.

Object (8):

- bright satellite streak in 4816
- 4812, 4814 has spurrious streak and poorer seeing
- mild background variations
- large background variations and spurrious dust in i

Object (9)

- 

