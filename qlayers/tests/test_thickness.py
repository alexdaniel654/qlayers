import numpy as np
import pandas as pd
import pytest

from numpy.random import Generator, PCG64
from qlayers.thickness import (
    logistic,
    gaussian,
    estimate_logistic_params,
    estimate_gaussian_params,
    equation_system,
    cortical_thickness,
)


class TestLogistic:
    def test_logistic_function_returns_expected_values(self):
        assert np.isclose(logistic(0, 1, 0, 1), 0.5)
        assert np.isclose(logistic(1, 1, 0, 1), 0.7310585786300049)


class TestGaussian:
    def test_gaussian_function_returns_expected_values(self):
        assert np.isclose(gaussian(0, 1, 0, 1), 1)
        assert np.isclose(gaussian(1, 1, 0, 1), 0.6065306597126334)


class TestEstimateLogisticParams:
    def test_estimate_params_for_logistic_function(self):
        x = np.array([0, 1, 2, 3, 4, 5])
        y = logistic(x, 1000, 5, -1)
        params, err = estimate_logistic_params(x, y)
        assert np.allclose(params, [1000, 5, -1])


class TestEstimateGaussianParams:
    def test_estimate_params_for_gaussian_function(self):
        x = np.array([0, 1, 2, 3, 4, 5])
        y = gaussian(x, 1000, 10, 2)
        params, err = estimate_gaussian_params(x, y)
        assert np.allclose(params, [1000, 10, 2])


class TestEquationSystem:
    def test_equation_system_returns_zero_for_same_params(self):
        assert equation_system(0, 1, 0, 1, 1, 0, 1) == -0.5


class TestCorticalThickness:
    def test_cortical_thickness_raises_error_for_wrong_space(self):
        class MockQLayers:
            def __init__(self):
                self.space = "map"

        with pytest.raises(ValueError):
            cortical_thickness(MockQLayers())

    def test_cortical_thickness_raises_error_for_missing_tissue_column(self):
        class MockQLayers:
            def __init__(self):
                self.space = "layers"

            def get_df(self, _):
                return pd.DataFrame({"depth": [0, 1, 2, 3, 4, 5]})

        with pytest.raises(ValueError):
            cortical_thickness(MockQLayers())

    def test_cortical_thickness_raises_error_for_wrong_tissue_column(self):
        class MockQLayers:
            def __init__(self):
                self.space = "layers"

            def get_df(self, _):
                return pd.DataFrame(
                    {
                        "depth": [0, 1, 2, 3, 4, 5],
                        "tissue": [
                            "Left",
                            "Left",
                            "Left",
                            "Right",
                            "Right",
                            "Right",
                        ],
                    }
                )

        with pytest.raises(ValueError):
            cortical_thickness(MockQLayers())

    def test_cortical_thickness_returns_expected_value_no_error(self):
        class MockQLayers:
            def __init__(self):
                self.space = "layers"

            def get_df(self, _):
                # Range of depths
                x = np.linspace(0, 20, 50)

                # Distributions to draw from
                cortex_dist = logistic(x, 500, 10, -0.4)
                medulla_dist = gaussian(x, 300, 10, 4)

                # Number of samples from each tissue type
                n = 1000

                # Draw samples
                rng = Generator(PCG64(seed=0))
                cortex_depths = rng.choice(
                    x, size=n, p=cortex_dist / cortex_dist.sum()
                )
                medulla_depths = rng.choice(
                    x, size=n, p=medulla_dist / medulla_dist.sum()
                )
                df_wide = pd.DataFrame(
                    {
                        "depth": np.concatenate(
                            (cortex_depths, medulla_depths)
                        ),
                        "tissue": np.repeat(["Cortex", "Medulla"], n),
                    }
                )
                return df_wide

        thickness = cortical_thickness(MockQLayers(), est_error=False)
        assert np.isclose(thickness, 6.96579)

    def test_cortical_thickness_returns_expected_value_with_error(self):
        class MockQLayers:
            def __init__(self):
                self.space = "layers"

            def get_df(self, _):
                # Range of depths
                x = np.linspace(0, 20, 1000)

                # Distributions to draw from
                cortex_dist = logistic(x, 500, 10, -0.4)
                medulla_dist = gaussian(x, 300, 10, 4)

                # Number of samples from each tissue type
                n = 10000

                # Draw samples
                rng = Generator(PCG64(seed=0))
                cortex_depths = rng.choice(
                    x, size=n, p=cortex_dist / cortex_dist.sum()
                )
                medulla_depths = rng.choice(
                    x, size=n, p=medulla_dist / medulla_dist.sum()
                )
                df_wide = pd.DataFrame(
                    {
                        "depth": np.concatenate(
                            (cortex_depths, medulla_depths)
                        ),
                        "tissue": np.repeat(["Cortex", "Medulla"], n),
                    }
                )
                return df_wide

        thickness, thickness_err = cortical_thickness(MockQLayers(), est_error=True)
        print(f"Thickness: {thickness}, Error: {thickness_err}")
        assert np.isclose(thickness, 7.10116, atol=1e-4)
        assert np.isclose(thickness_err, 0.11867, atol=1e-4)
