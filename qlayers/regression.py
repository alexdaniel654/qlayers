import numpy as np
import pandas as pd

from scipy.stats import linregress


def slope(qlayers, maps="all", outer=5.0, inner=15.0, unit="mm", agg=np.nanmedian):
    """
    Explore the profiles produced by the layers algorithm. Calculate summary
    statistics for the outer layers, inner layers and the gradient between the
    two, these are similar to cortex, medulla and cortical-medullary
    difference respectively.

    Parameters
    ----------
    qlayers : QLayers
        The QLayers object containing the depth and quantitative map data.
    maps : str or list of str, optional
        Default "all"
        The names of the maps to calculate the slope for.
        If "all", all maps in the QLayers object are used.
    outer : float, optional
        Default 5.0
        The outer boundary of the depth range to calculate the slope within.
        This should be approximately the depth of the cortex. It can be
        specified in mm, percent of the maximum depth (i.e. between 0
        and 100) or as a proportion of the maximum depth (i.e. between 0
        and 1).

    inner : float, optional
        Default 15.0
        The inner boundary of the depth range to calculate the slope
        within. It can be specified in mm, percent of the maximum depth
        (i.e. between 0 and 100) or as a proportion of the maximum depth
        (i.e. between 0 and 1).
    unit : str, optional
        The unit of the outer and inner parameters. Can be "mm", "percent",
         or "prop". Default is "mm".
    agg : function, optional
        The aggregation function to use when calculating the summary
        value for the outer and inner regions. Default is np.median.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the calculated slopes and other statistics for
        each map.
    """
    if qlayers.space != "layers":
        raise ValueError(
            "Cortical thickness can only be computed if the "
            "QLayers object is in layers space"
        )
    if unit not in ["mm", "percent", "prop"]:
        raise ValueError("unit must be one of 'mm', 'percent' or 'prop'")

    df = qlayers.get_df("long")
    df = df.dropna()

    if maps == "all":
        maps = df['measurement'].unique().tolist()

    if unit == "percent":
        outer = (outer / 100) * df["depth"].max()
        inner = (inner / 100) * df["depth"].max()
    elif unit == "prop":
        outer = outer * df["depth"].max()
        inner = inner * df["depth"].max()

    if type(maps) is str:
        maps = [maps]
    results_df = pd.DataFrame(
        index=maps,
        columns=["inner", "outer", "grad", "inner_std", "outer_std", "grad_se"],
    )
    for m in maps:
        if m not in df['measurement'].unique():
            raise ValueError(f"{m} is not a valid map")

        results_df.loc[m, "outer"] = agg(df.loc[(df["depth"] < outer) &
                                                (df["measurement"] == m),
                                                'value'])
        results_df.loc[m, "inner"] = agg(df.loc[(df["depth"] > inner) &
                                                (df["measurement"] == m),
                                                'value'])
        reg = linregress(
            df.loc[(df["depth"] > outer) &
                   (df["depth"] < inner) &
                   (df["measurement"] == m), "depth"],
            df.loc[(df["depth"] > outer) &
                   (df["depth"] < inner) &
                   (df["measurement"] == m), 'value'],
        )
        results_df.loc[m, "grad"] = reg.slope
        results_df.loc[m, "outer_std"] = (
            np.std(df.loc[(df["depth"] < outer) &
                          (df["measurement"] == m),
                          'value']))
        results_df.loc[m, "inner_std"] = (
            np.std(df.loc[(df["depth"] > inner) &
                          (df["measurement"] == m),
                          'value']))
        results_df.loc[m, "grad_se"] = reg.stderr

    return results_df
