import numpy as np
import bottleneck
import math

#These are slightly modified functions from Courty's code
#primarily from: https://github.com/lrntct/pxr/blob/master/ev_fit.py


iscalar = np.int32
fscalar = np.float32
shape_param = None
shape_param = -0.114

# Euler-Mascheroni constant
EM = fscalar(0.577215664901532860606512090082)

def fit_gev(arr):
    n_obs = len(arr)
    ax_year = 0
    # rank samples
    rank = bottleneck.nanrankdata(arr, axis=ax_year).astype(fscalar)
    # fit distribution. ev_params is a tuple of ndarrays.
    ecdf = comp_ecdf(rank, n_obs)
    #print(ecdf)
    ev_params = gev_pwm(arr, ecdf, n_obs, ax_year, shape=shape_param)
    # Convert the output to a tuple of floats
    ev_params_floats = (float(ev_params[0]), float(ev_params[1]), float(ev_params[2]))
    print(ev_params_floats)
    return ev_params_floats

def comp_ecdf(rank, n_obs):
    """Return the ECDF
    Recommended as an unbiased estimator for PWM by:
    Hosking, J. R. M., and J. R. Wallis. 1995.
    “A Comparison of Unbiased and Plotting-Position Estimators of L Moments.”
    Water Resources Research 31 (8): 2019–25.
    https://doi.org/10.1029/95WR01230.
    """
    return rank / n_obs

def gev_pwm(ams, ecdf, n_obs, ax_year, shape=None):
    """Fit the GEV using the Method of Probability-Weighted Moments.
    EV type II when shape<0.

    Hosking, J. R. M., Wallis, J. R., & Wood, E. F. (1985).
    Estimation of the Generalized Extreme-Value Distribution
    by the Method of Probability-Weighted Moments.
    Technometrics, 27(3), 251–261.
    https://doi.org/10.1080/00401706.1985.10488049
    """
    # b0 = np.mean(ams, axis=ax_year)
    b0 = gen_bvalue(ecdf, ams, n_obs, iscalar(0), ax_year)
    b1 = gen_bvalue(ecdf, ams, n_obs, iscalar(1), ax_year)
    l1 = b0
    l2 = iscalar(2) * b1 - b0
    if shape is not None:
        arr_shape = np.atleast_1d(np.full_like(b0, shape))
    else:
        b2 = gen_bvalue(ecdf, ams, n_obs, iscalar(2), ax_year)
        arr_shape = gev_shape(b0, b1, b2)
    arr_scale = np.where(arr_shape == 0,
                         gumbel_scale(l2),
                         gev_scale(l2, arr_shape))
    arr_loc = np.where(arr_shape == 0,
                       gumbel_loc(l1, arr_scale),
                       gev_loc(l1, arr_scale, arr_shape))
    return arr_loc, arr_scale, arr_shape

def gen_bvalue(ecdf, ams, n_obs, order, axis):
    """Estimation of bvalue not depending on xarray
    """
    pr_sum = (ecdf**order * ams).sum(axis=axis)
    bvalue = pr_sum / n_obs
    return bvalue

def gev_shape(b0, b1, b2):
    """Shape parameter of the GEV
    EV type II when shape<0.
    Hosking, J. R. M., Wallis, J. R., & Wood, E. F. (1985).
    Estimation of the Generalized Extreme-Value Distribution
    by the Method of Probability-Weighted Moments.
    Technometrics, 27(3), 251–261.
    https://doi.org/10.1080/00401706.1985.10488049
    """
    c = (iscalar(2)*b1 - b0) / (iscalar(3)*b2 - b0) - fscalar(0.63093)  # ln(2) / ln(3)
    return fscalar(7.859) * c + fscalar(2.9554) * c**iscalar(2)

def gumbel_scale(l2):
    return l2 / fscalar(0.6931)  # ln(2)

def gev_scale(l2, shape):
    """Scale parameter of the GEV
    EV type II when shape<0.
    Hosking, J. R. M., & Wallis, J. R. (1997).
    Appendix: L-moments for some specific distributions.
    In Regional Frequency Analysis (pp. 191–209).
    Cambridge: Cambridge University Press.
    https://doi.org/10.1017/CBO9780511529443.012
    """
    return (l2 * shape) / (gamma(iscalar(1) + shape) * (iscalar(1) - iscalar(2)**-shape))

def gamma(x):
    return math.gamma(x)

def gumbel_loc(l1, scale):
    return l1 - EM * scale

def gev_loc(l1, scale, shape):
    """Location parameter of the GEV
    EV type II when shape<0.
    Hosking, J. R. M., & Wallis, J. R. (1997).
    Appendix: L-moments for some specific distributions.
    In Regional Frequency Analysis (pp. 191–209).
    Cambridge: Cambridge University Press.
    https://doi.org/10.1017/CBO9780511529443.012
    """
    return l1 - scale * ((iscalar(1) - gamma(iscalar(1) + shape)) / shape)