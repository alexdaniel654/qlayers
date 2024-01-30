import pytest
from qlayers.regression import slope
from qlayers import QLayers
import numpy as np
import pandas as pd

np.random.seed(0)


class TestSlope:
    class MockQLayers:
        def __init__(self):
            self.space = "layers"
            self.maps = ["t1", "t2"]

        def get_df(self, format="wide"):
            if format == "long":
                df = pd.DataFrame(
                    {
                        "depth": np.tile(np.linspace(0, 20, 50), 2),
                        "measurement": np.repeat(["t1", "t2"], 50),
                        "value": np.concatenate((np.linspace(1000, 2000, 50)
                                                 + np.random.randn(50) * 10,
                                                 np.linspace(100, 200,50)
                                                 + np.random.randn(50))),
                    }
                )
            elif format == "wide":
                df = pd.DataFrame(
                {
                    "depth": np.linspace(0, 20, 50),
                    "t1": np.linspace(1000, 2000, 50)
                    + np.random.randn(50) * 10,
                    "t2": np.linspace(100, 200, 50) + np.random.randn(50),
                }
            )
            return df

    def test_slope_with_invalid_space(self):
        qlayers = self.MockQLayers()
        qlayers.space = "map"

        with pytest.raises(ValueError):
            slope(qlayers)

    def test_slope_with_invalid_unit(self):
        with pytest.raises(ValueError):
            slope(self.MockQLayers(), unit="inches")

    def test_slope_with_invalid_map(self):
        with pytest.raises(ValueError):
            slope(self.MockQLayers(), maps="invalid_map")

    def test_slope_with_all_maps(self):
        result = slope(self.MockQLayers(), maps="all")
        assert "t1" in result.index
        assert "t2" in result.index

    def test_slope_with_specific_map(self):
        result = slope(self.MockQLayers(), "t1")
        assert "t1" in result.index
        assert "t2" not in result.index

    def test_slope_output_values(self):
        result = slope(self.MockQLayers(), "t1")
        assert np.isclose(result.loc["t1", "inner"], 1865.41024)
        assert np.isclose(result.loc["t1", "outer"], 1109.10639)
        assert np.isclose(result.loc["t1", "grad"], 49.62542)
        assert np.isclose(result.loc["t1", "inner_std"], 75.91713)
        assert np.isclose(result.loc["t1", "outer_std"], 76.70788)
        assert np.isclose(result.loc["t1", "grad_se"], 0.73361)

    def test_slope_boundaries(self):
        result = slope(self.MockQLayers(), "t1", outer=2.0, inner=18.0)
        assert np.isclose(result.loc["t1", "inner"], 1949.53755)
        assert np.isclose(result.loc["t1", "outer"], 1048.48295)
        assert np.isclose(result.loc["t1", "grad"], 49.81742)
        assert np.isclose(result.loc["t1", "inner_std"], 27.54807)
        assert np.isclose(result.loc["t1", "outer_std"], 29.13423)
        assert np.isclose(result.loc["t1", "grad_se"], 0.3219759)

    def test_slope_with_unit_percent(self):
        result = slope(self.MockQLayers(), unit="percent", outer=20, inner=80)
        assert np.isclose(result.loc["t1", "inner"], 1904.02085)
        assert np.isclose(result.loc["t1", "outer"], 1087.19313)
        assert np.isclose(result.loc["t1", "grad"], 49.68695)
        assert np.isclose(result.loc["t1", "inner_std"], 59.09864)
        assert np.isclose(result.loc["t1", "outer_std"], 56.31960)
        assert np.isclose(result.loc["t1", "grad_se"], 0.5826736)

    def test_slope_with_unit_prop(self):
        result = slope(self.MockQLayers(), unit="prop", outer=0.15, inner=0.85)
        assert np.isclose(result.loc["t1", "inner"], 1932.19128)
        assert np.isclose(result.loc["t1", "outer"], 1065.58927)
        assert np.isclose(result.loc["t1", "grad"], 50.49024)
        assert np.isclose(result.loc["t1", "inner_std"], 50.44155)
        assert np.isclose(result.loc["t1", "outer_std"], 46.29691)
        assert np.isclose(result.loc["t1", "grad_se"], 0.3447094)

    def test_slope_different_agg(self):
        result = slope(self.MockQLayers(), agg=np.mean)
        assert np.isclose(result.loc["t1", "inner"], 1874.41801)
        assert np.isclose(result.loc["t1", "outer"], 1123.73675)
        assert np.isclose(result.loc["t1", "grad"], 49.87365)
        assert np.isclose(result.loc["t1", "inner_std"], 79.41030)
        assert np.isclose(result.loc["t1", "outer_std"], 79.36174)
        assert np.isclose(result.loc["t1", "grad_se"], 0.7652260)
