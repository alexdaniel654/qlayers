import numpy as np

from numpy.random import Generator, PCG64
from scipy import stats
from scipy.optimize import curve_fit, fsolve


def logistic(x, L, x0, k):
    """
    Logistic function.

    Parameters
    ----------
    x : numpy.ndarray
        The input array.
    L : float
        The curve's maximum value.
    x0 : float
        The x-value of the sigmoid's midpoint.
    k : float
        The logistic growth rate or steepness of the curve.

    Returns
    -------
    numpy.ndarray
        The logistic function applied to x.
    """
    return L / (1 + np.exp(-k * (x - x0)))


def gaussian(x, a, b, c):
    """
    Gaussian function.

    Parameters
    ----------
    x : numpy.ndarray
        The input array.
    a : float
        The height of the curve's peak.
    b : float
        The position of the center of the peak.
    c : float
        Controls the width of the curve.

    Returns
    -------
    numpy.ndarray
        The Gaussian function applied to x.
    """
    return a * np.exp(-((x - b) ** 2) / (2 * c ** 2))


def estimate_logistic_params(x, y):
    """
    Estimates the parameters of a logistic function based on some sample data.

    Parameters
    ----------
    x : numpy.ndarray
        The x data.
    y : numpy.ndarray
        The y data.

    Returns
    -------
    params : numpy.ndarray
        The estimated parameters [L, x0, k].
    err : numpy.ndarray
        The uncertainty in the estimated parameters.
    """
    # Provide some initial guess values for L, x0, k
    initial_guess = [max(y), np.median(x), -1]

    # Use curve_fit to estimate the logistic function parameters
    params, params_covariance = curve_fit(logistic, x, y, p0=initial_guess,
                                          bounds=([0, 0, -np.inf], [np.inf, 500, 0]))
    err = np.sqrt(np.diag(params_covariance))
    return params, err


def estimate_gaussian_params(x, y):
    """
    Estimates the parameters of a Gaussian function based on some sample data.

    Parameters
    ----------
    x : numpy.ndarray
        The x data.
    y : numpy.ndarray
        The y data.

    Returns
    -------
    params : numpy.ndarray
        The estimated parameters [a, b, c].
    err : numpy.ndarray
        The uncertainty in the estimated parameters.
    """
    # Provide some initial guess values for a, b, c
    initial_guess = [max(y), np.median(x), np.std(x)]

    # Use curve_fit to estimate the Gaussian function parameters
    params, params_covariance = curve_fit(gaussian, x, y, p0=initial_guess,
                                          bounds=([0, 0, 0], [np.inf, 500, np.inf]))
    err = np.sqrt(np.diag(params_covariance))
    return params, err


def equation_system(variable, L, x0, k, a, b, c):
    """
    Equation system for logistic and Gaussian functions.

    Parameters
    ----------
    variable : float
        The variable to solve for.
    L : float
        The curve's maximum value for logistic function.
    k : float
        The logistic growth rate or steepness of the curve for logistic function.
    x0 : float
        The x-value of the sigmoid's midpoint for logistic function.
    a : float
        The height of the curve's peak for Gaussian function.
    b : float
        The position of the center of the peak for Gaussian function.
    c : float
        Controls the width of the curve for Gaussian function.

    Returns
    -------
    float
        The difference between logistic and Gaussian functions.
    """
    x = variable
    eq = logistic(x, L, x0, k) - gaussian(x, a, b, c)
    return eq


def cortical_thickness(qlayers, est_error=False):
    """
    Computes the cortical thickness of the kidneys.

    Parameters
    ----------
    qlayers : object
        The QLayers object.
    est_error : bool, optional
        Default False
        If True, estimate the error in the cortical thickness, this is done using a stochastic approach. Lots of samples
        are drawn from the estimated parameters and the cortical thickness is calculated for each sample.
        The mean and standard deviation of the cortical thickness is returned. Estimating errors is therefore slow and
        off by default.

    Returns
    -------
    thickness : float
        The cortical depth.
    thickness_err : float
        The uncertainty in the cortical depth.
    """
    if qlayers.space != "layers":
        raise ValueError(
            "Cortical thickness can only be computed if the "
            "QLayers object is in layers space"
        )
    df = qlayers.get_df("wide")
    df = df.dropna()
    df = df[df["depth"] > 0]

    if "tissue" not in df.columns:
        raise ValueError(
            "Cortical thickness can only be computed if tissue "
            "labels have been added to the QLayers object"
        )
    if not sorted(list(map(str.lower, df["tissue"].unique()))) == ["cortex", "medulla"]:
        raise ValueError(
            "Cortical thickness can only be computed if tissue "
            'labels are "Cortex" and "Medulla"'
        )

    # Convert samples into a distribution that the curves can be fit to
    bin_width = 0.5
    bins = np.arange(0, df["depth"].max() + bin_width, bin_width)
    density_cortex, bins = np.histogram(
        df.loc[df["tissue"] == "Cortex", "depth"], bins=bins
    )
    density_medulla, bins = np.histogram(
        df.loc[df["tissue"] == "Medulla", "depth"], bins=bins
    )

    # Fit a logistic function to the cortex data and a Gaussian function to the
    # medulla data
    x = (bins + bin_width / 2)[:-1]
    params_cortex, err_cortex = estimate_logistic_params(x, density_cortex)
    params_medulla, err_medulla = estimate_gaussian_params(x, density_medulla)

    # Estimate the intersection of the two curves, this is the point at
    # which more voxels are medulla than cortex and can be taken as the
    # cortical thickness
    cortical_depth = fsolve(
        equation_system, [5.0], args=(*params_cortex, *params_medulla)
    )
    if est_error:
        n_samp = 1000
        rng = Generator(PCG64(seed=0))
        cortex_samples = rng.normal(params_cortex, err_cortex, (n_samp, len(params_cortex)))
        medulla_samples = rng.normal(params_medulla, err_medulla, (n_samp, len(params_medulla)))
        cortical_depth_samples = np.zeros(n_samp)
        for i in range(n_samp):
            roots = fsolve(equation_system, np.linspace(0, df['depth'].max(), 10),
                           args=(*cortex_samples[i], *medulla_samples[i])
                           )
            roots = roots[0 < roots]
            roots = roots[roots < df["depth"].max()]
            if len(roots) > 0:
                cortical_depth_samples[i] = stats.mode(roots, keepdims=False)[0]
            else:
                cortical_depth_samples[i] = np.nan
        thickness_err = np.nanstd(cortical_depth_samples)
        thickness_mean = np.nanmean(cortical_depth_samples)
        return thickness_mean, thickness_err
    else:
        return cortical_depth[0]
