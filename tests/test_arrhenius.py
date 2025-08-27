"""
Tests for Arrhenius analysis functionality.
"""

import pytest
import numpy as np
from unittest.mock import Mock

from garnet_analysis.core.arrhenius import ArrheniusAnalyzer, ArrheniusDataPoint, ArrheniusResult


class TestArrheniusDataPoint:
    """Test ArrheniusDataPoint dataclass."""
    
    def test_valid_data_point(self):
        """Test creating a valid data point."""
        point = ArrheniusDataPoint(temperature=1000, diffusivity=1e-6, std_error=1e-7)
        
        assert point.temperature == 1000
        assert point.diffusivity == 1e-6
        assert point.std_error == 1e-7
    
    def test_negative_temperature_raises_error(self):
        """Test that negative temperature raises error."""
        with pytest.raises(ValueError, match="Temperature must be positive"):
            ArrheniusDataPoint(temperature=-100, diffusivity=1e-6, std_error=1e-7)
    
    def test_negative_diffusivity_raises_error(self):
        """Test that negative diffusivity raises error."""
        with pytest.raises(ValueError, match="Diffusivity cannot be negative"):
            ArrheniusDataPoint(temperature=1000, diffusivity=-1e-6, std_error=1e-7)


class TestArrheniusAnalyzer:
    """Test ArrheniusAnalyzer class."""
    
    def test_analyzer_initialization(self):
        """Test basic analyzer initialization."""
        analyzer = ArrheniusAnalyzer()
        
        assert len(analyzer.data_points) == 0
        assert analyzer.fit_result is None
    
    def test_add_data_point(self):
        """Test adding data points."""
        analyzer = ArrheniusAnalyzer()
        analyzer.add_data_point(1000, 1e-6, 1e-7)
        
        assert len(analyzer.data_points) == 1
        assert analyzer.data_points[0].temperature == 1000
    
    def test_add_data_from_dict(self):
        """Test adding data from dictionary format."""
        analyzer = ArrheniusAnalyzer()
        data_dict = {
            1000: {'diff': 1e-6, 'std': 1e-7},
            800: {'diff': 5e-7, 'std': 5e-8}
        }
        
        analyzer.add_data_from_dict(data_dict)
        
        assert len(analyzer.data_points) == 2
        temperatures = [point.temperature for point in analyzer.data_points]
        assert 1000 in temperatures
        assert 800 in temperatures
    
    def test_fit_insufficient_data(self):
        """Test that fitting fails with insufficient data."""
        analyzer = ArrheniusAnalyzer()
        analyzer.add_data_point(1000, 1e-6, 1e-7)
        
        with pytest.raises(ValueError, match="Need at least 2 data points"):
            analyzer.fit()
    
    def test_fit_with_valid_data(self):
        """Test fitting with valid data."""
        analyzer = ArrheniusAnalyzer()
        
        # Add some synthetic data that should fit well
        temperatures = [1000, 900, 800, 700, 600]
        # Generate data following Arrhenius behavior
        ea = 0.3  # eV
        d0 = 1e-3
        kb = 8.617e-5  # eV/K
        
        for temp in temperatures:
            diff = d0 * np.exp(-ea / (kb * temp))
            analyzer.add_data_point(temp, diff, diff * 0.1)  # 10% error
        
        result = analyzer.fit()
        
        assert result is not None
        assert isinstance(result, ArrheniusResult)
        assert result.n_points == 5
        assert result.r_squared > 0.9  # Should be a good fit
        assert 0.2 < result.activation_energy < 0.4  # Should recover original Ea
    
    def test_extrapolate_without_fit(self):
        """Test that extrapolation fails without fitting."""
        analyzer = ArrheniusAnalyzer()
        
        with pytest.raises(ValueError, match="Must fit data before extrapolation"):
            analyzer.extrapolate_to_temperature(298)
    
    def test_extrapolate_after_fit(self):
        """Test extrapolation after successful fit."""
        analyzer = ArrheniusAnalyzer()
        
        # Add data and fit
        temperatures = [1000, 800, 600]
        for temp in temperatures:
            diff = 1e-6 * np.exp(-3000 / temp)  # Simple exponential
            analyzer.add_data_point(temp, diff, diff * 0.1)
        
        result = analyzer.fit()
        d_298k, d_min, d_max = analyzer.extrapolate_to_temperature(298)
        
        assert d_298k > 0
        assert d_min < d_298k < d_max
        assert d_min > 0
    
    def test_get_data_arrays(self):
        """Test getting data as arrays."""
        analyzer = ArrheniusAnalyzer()
        analyzer.add_data_point(1000, 1e-6, 1e-7)
        analyzer.add_data_point(800, 5e-7, 5e-8)
        
        x_data, y_data, y_errors = analyzer.get_data_arrays()
        
        assert len(x_data) == 2
        assert len(y_data) == 2
        assert len(y_errors) == 2
        
        # Check that x_data is 1000/T
        expected_x = [1000/1000, 1000/800]
        np.testing.assert_array_almost_equal(sorted(x_data), sorted(expected_x))


class TestArrheniusResult:
    """Test ArrheniusResult functionality."""
    
    def create_mock_result(self):
        """Create a mock ArrheniusResult for testing."""
        return ArrheniusResult(
            slope=-3000,
            intercept=2,
            slope_error=100,
            intercept_error=0.1,
            activation_energy=0.3,
            activation_energy_error=0.01,
            r_squared=0.95,
            n_points=5,
            temperature_range=(600, 1000)
        )
    
    def test_extrapolation_method(self):
        """Test the extrapolation method of ArrheniusResult."""
        result = self.create_mock_result()
        
        d_298k, d_min, d_max = result.extrapolate_to_temperature(298)
        
        assert d_298k > 0
        assert d_min < d_298k < d_max
        assert d_min > 0


if __name__ == '__main__':
    pytest.main([__file__])