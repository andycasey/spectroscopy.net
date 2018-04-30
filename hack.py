
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from astropy.io import fits
from astropy.table import Table
from itertools import permutations


known_arc_lines = Table.read("lovis_et_al_2007.fits")[:1000]

# OK, let's measure positions for three lines.
thar = fits.open("thar/HARPS.2004-01-01T20:07:02.313_e2ds_A.fits")

# Calculate the wavelength calibration given the coefficients in the header.
thar_intensities = thar[0].data
thar_wavelengths = np.zeros_like(thar_intensities)
for o in range(thar_intensities.shape[0]):

    d = thar[0].header["ESO DRS CAL TH DEG LL"]
    A = thar[0].header["ESO DRS CAL TH COEFF LL"] 

    raise a

thar_index = 0
intensities = thar[0].data[thar_index]

def gaussian(x, *parameters):
    """
    Evaluate a Gaussian profile at x, given the profile parameters.

        y = amplitude * exp(-(x - position)**2 / (2.0 * sigma**2))

    :param x:
        The x-values to evaluate the Gaussian profile at.

    :param parameters:
        The position, sigma, and amplitude of the Gaussian profile.
    """
    position, sigma, amplitude, background = parameters
    return amplitude * np.exp(-(x - position)**2 / (2.0 * sigma**2)) + background


def fit_arc_line(intensities, pixel_min, pixel_max, **kwargs):

    pixels = np.arange(len(intensities))
    region = (pixel_max >= pixels) * (pixels >= pixel_min)

    x, y = (pixels[region], intensities[region])

    hw = 0.5 * (pixel_max - pixel_min)
    x0 = pixel_min + hw

    kwds = dict(
        p0=[x0, 2.0, np.max(y), 0],
        bounds=np.array([
            [pixel_min, 0, 0, 0],
            [pixel_max, hw, 1.5 * np.max(y), np.max(y)]
        ]))
    p_opt, p_cov = op.curve_fit(gaussian, x, y, **kwds)

    assert not np.all(p_opt == kwds["p0"])
    meta = dict(x=x, p0=kwds["p0"])
    return (p_opt, p_cov, meta)


"""
measured_arc_lines = [
    fit_arc_line(intensities, 966, 1005),
    fit_arc_line(intensities, 3210, 3230),
    fit_arc_line(intensities, 3500, 3540)
]
"""

measured_arc_lines = [
    fit_arc_line(intensities, 2097, 2106),
    fit_arc_line(intensities, 2106, 2117),
    fit_arc_line(intensities, 2720, 2750)
]


fig, ax = plt.subplots()
ax.plot(intensities)
ax.set_xlabel("Pixel")
ax.set_ylabel("Intensity")
fig.tight_layout()

for p_opt, p_cov, meta in measured_arc_lines:
    ax.plot(meta["x"], gaussian(meta["x"], *p_opt), c='r')


def measured_pixel_ratios(measured_pixel_peaks, indices=None):

    if indices is None:
        values = np.array(measured_pixel_peaks)
        indices = np.argsort(values)

    else:
        indices = np.sort(indices)
        values = measured_pixel_peaks[indices]


    min_value, max_value = (np.min(values), np.max(values))
    ratios = (values - min_value)/(max_value - min_value)

    offset = min_value
    scale = max_value - min_value
    return (ratios, indices, offset, scale)

measured_pixel_peaks = np.array([po[0] for po, pc, mt in measured_arc_lines])
pixel_ratios, pixel_indices, pixel_offset, pixel_scale = measured_pixel_ratios(measured_pixel_peaks)



# From the line list, go through all permutations starting with the brightest
def arc_line_permutations(known_arc_lines, N, wavelength_label="Lambda",
    intensity_label="Intensity"):

    """
    all_permutations = np.array(list(permutations(np.arange(M), N)))
    sum_of_permutation_intensities = np.sum(
        known_arc_lines[intensity_label][all_permutations], axis=1)
    brightness_indices = np.argsort(sum_of_permutation_intensities)[::-1]
    """

    # Take the ~brightest combinations first.
    brightness_indices = np.argsort(known_arc_lines[intensity_label])[::-1]

    tolerance = 1e-5

    best = []
    lowest = np.inf

    for i, permutation_indices in enumerate(permutations(brightness_indices, N)):

        # Measure the ratios.
        ratios, permutation_indices, offset, scale = measured_pixel_ratios(
            known_arc_lines[wavelength_label], permutation_indices)

        sum_diff = np.sum(np.abs(pixel_ratios - ratios))
        if len(best) == 0 or sum_diff < lowest:
            best = [ratios, permutation_indices]

        if np.allclose(pixel_ratios, ratios, atol=tolerance):

            # Make a proposal.

            #  convert all wavelengths to pixels
            all_ratios = (known_arc_lines[wavelength_label] - offset)/scale

            proposed_pixels = all_ratios * np.ptp(measured_pixel_peaks) + np.min(measured_pixel_peaks)


            # If this match is true, then we should predict the wavelengths of
            # all other lines.

            # calculate ratios for other points.

            # pixel = M * wavelength + B
            # pixel = (pixels/wavelength_ptp) * wavelength + offset
            #foo = np.ptp(permutation_indices)/np.ptp(known_arc_lines[wavelength_label][permutation_indices]) * known_arc_lines[wavelength_label] + known_arc_lines[wavelength_label][permutation_indices][0]



            #foo = (known_arc_lines[wavelength_label] - offset)/scale

            #proposed_pixels = foo * np.ptp(measured_pixel_peaks) \
            #                + np.min(measured_pixel_peaks)


            fig, ax = plt.subplots()
            ax.plot(intensities, c='k')

            ax.scatter(proposed_pixels, np.ones(proposed_pixels.size) * np.array(ax.get_ylim()).mean(), c="r")
            ax.set_xlim(0, len(intensities))

            #    ax.axvline(proposed_pixel, c="#666666", lw=1, linestyle=":")

            raise a


        print("{:.0f} {:.2f} {:.2f} {:.5e}".format(i, ratios[1], pixel_ratios[1],
            np.sum(known_arc_lines[intensity_label][permutation_indices])))


arc_line_permutations(known_arc_lines, len(measured_pixel_peaks))







# For each permutation, see how close the ratio is to the observed one.

# If the ratio is within some error (that we specify in wavelength space),
# then make a proposal.

# Evaluate a proposal.


# https://www.dropbox.com/sh/zw63o3qiqamc1we/AAARKM3b268mCTM_WpBPi6xta?dl=0

# Load in an arc lamp, and let's measure the positions of three lines.

