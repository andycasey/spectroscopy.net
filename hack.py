
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from astropy.io import fits
from astropy.table import Table
from itertools import permutations


known_arc_lines = Table.read("lovis_et_al_2007.fits")

keep = (5000 > known_arc_lines["Lambda"]) \
     * (known_arc_lines["Lambda"] > 4000)

known_arc_lines = known_arc_lines[keep]

# OK, let's measure positions for three lines.
thar = fits.open("thar/HARPS.2004-01-01T20:07:02.313_e2ds_A.fits")

# Calculate the wavelength calibration given the coefficients in the header.
O, P = thar[0].data.shape

D = 1 + thar[0].header["ESO DRS CAL TH DEG X"]
x = np.arange(P) # 1-based indexing for HARPS

coefficients = np.array([
    thar[0].header["ESO DRS CAL TH COEFF LL{:.0f}".format(i)] \
    for i in range(D * O)])

design_matrix = np.array([x**i for i in range(D)]).T
thar_wavelengths = np.dot(design_matrix, coefficients.reshape((-1, D)).T).T
thar_intensities = thar[0].data


def air_to_vacuum(wavelength):
    """
    It converts spectrum's wavelengths (nm) from air to vacuum
    """
    # Following the air to vacuum conversion from VALD3 (computed by N. Piskunov) http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    s_square = np.power(1.e4 / wavelength, 2)
    n2 = 1. + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s_square) + 0.0001599740894897 / (38.92568793293 - s_square)
    return wavelength*n2 # Angstroms


def vacuum_to_air(wavelength):
    """
    It converts spectrum's wavelengths from vacuum to air
    """
    # Following the vacuum to air conversion the formula from Donald Morton (2000, ApJ. Suppl., 130, 403) which is also a IAU standard
    # - More info: http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    s_square = np.power(1.e4 / wavelength, 2)
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s_square) + 0.00015998 / (38.9 - s_square)
    return wavelength/n # Angstroms


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

order_index = 30
wavelengths = air_to_vacuum(thar_wavelengths[order_index])
intensities = thar_intensities[order_index]


# For order_index = 30
measured_arc_lines = [
    fit_arc_line(intensities, 1285, 1300),
    fit_arc_line(intensities, 1798, 1810),
    fit_arc_line(intensities, 3899, 3909)
]

measured_arc_lines = [
    fit_arc_line(intensities, 1535, 1545),
    fit_arc_line(intensities, 1555, 1570),
    fit_arc_line(intensities, 1660, 1671)
]



fig, axes = plt.subplots(2)
axes[0].plot(intensities)
axes[0].set_xlabel("Pixel")
axes[0].set_ylabel("Intensity")
axes[0].set_xlim(0, P)

for p_opt, p_cov, meta in measured_arc_lines:
    axes[0].plot(meta["x"], gaussian(meta["x"], *p_opt), c='r')


axes[1].plot(wavelengths, intensities)

ok = (wavelengths[-1] >= known_arc_lines["Lambda"]) \
   * (known_arc_lines["Lambda"] >= wavelengths[0])

for known_arc_line in known_arc_lines[ok]:
    axes[1].axvline(known_arc_line["Lambda"], c="#666666", linestyle=":", zorder=-1)

axes[1].set_xlim(wavelengths[0], wavelengths[-1])

axes[1].set_xlabel("Wavelenth (vac)")
axes[1].set_ylabel("Intensity")

fig.tight_layout()



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

    #for i, permutation_indices in enumerate(permutations(brightness_indices, N)):
    for i, permutation_indices in enumerate(np.array([[2203, 2204, 2211]])):

        # Measure the ratios.
        ratios, permutation_indices, offset, scale = measured_pixel_ratios(
            known_arc_lines[wavelength_label], permutation_indices)

        sum_diff = np.sum(np.abs(pixel_ratios - ratios))
        if len(best) == 0 or sum_diff < lowest:
            best = [ratios, permutation_indices]

        
        #if np.allclose(pixel_ratios, ratios, atol=tolerance):

        if np.all(np.array([2203, 2204, 2211]) == permutation_indices):


            # Make a proposal.
            arc_wavelengths = known_arc_lines[wavelength_label]
            model_ratios = (arc_wavelengths - arc_wavelengths[permutation_indices[0]]) \
                         / (arc_wavelengths[permutation_indices[-1]] - arc_wavelengths[permutation_indices[0]])




            #  convert all wavelengths to pixels
            #all_ratios = (known_arc_lines[wavelength_label] - offset)/scale

            proposed_pixels = model_ratios * np.ptp(measured_pixel_peaks) + np.min(measured_pixel_peaks)


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

            pixels = np.arange(intensities.size)
            ok = (intensities.size >= proposed_pixels) * (proposed_pixels >= 0)


            ax.scatter(proposed_pixels[ok], intensities[pixels.searchsorted(proposed_pixels[ok])], c='r')



            #ax.scatter(proposed_pixels, intensities[np.arange(intensities.size).searchsorted(proposed_pixels)], c="r")
            ax.set_xlim(0, len(intensities))

            #    ax.axvline(proposed_pixel, c="#666666", lw=1, linestyle=":")
            print(known_arc_lines[wavelength_label][permutation_indices])
            raise a


        print("{:.0f} {:.2f} {:.2f} {:.5e} {}".format(i, ratios[1], pixel_ratios[1],
            np.sum(known_arc_lines[intensity_label][permutation_indices]), permutation_indices))


arc_line_permutations(known_arc_lines, len(measured_pixel_peaks))







# For each permutation, see how close the ratio is to the observed one.

# If the ratio is within some error (that we specify in wavelength space),
# then make a proposal.

# Evaluate a proposal.


# https://www.dropbox.com/sh/zw63o3qiqamc1we/AAARKM3b268mCTM_WpBPi6xta?dl=0

# Load in an arc lamp, and let's measure the positions of three lines.

