"""
Tests for diffusivity analysis functionality.
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch

from garnet_analysis.core.diffusivity import DiffusivityAnalyzer, DiffusivityResult, ElementMSDData


class TestDiffusivityResult:
    """Test DiffusivityResult dataclass."""
    
    def test_valid_result_creation(self):
        """Test creating a valid DiffusivityResult."""
        result = DiffusivityResult(
            diffusivity=1e-6,
            std_error=1e-7,
            msd_slope=6.0,
            time_range=(10.0, 90.0),
            temperature=1000,
            species='Li',
            n_data_points=100,
            r_squared=0.95
        )
        
        assert result.diffusivity == 1e-6
        assert result.temperature == 1000
        assert result.species == 'Li'
    
    def test_negative_diffusivity_raises_error(self):
        """Test that negative diffusivity raises ValueError."""
        with pytest.raises(ValueError, match="Diffusivity cannot be negative"):
            DiffusivityResult(
                diffusivity=-1e-6,
                std_error=1e-7,
                msd_slope=6.0,
                time_range=(10.0, 90.0),
                temperature=1000,
                species='Li',
                n_data_points=100,
                r_squared=0.95
            )
    
    def test_zero_temperature_raises_error(self):
        """Test that zero temperature raises ValueError."""
        with pytest.raises(ValueError, match="Temperature must be positive"):
            DiffusivityResult(
                diffusivity=1e-6,
                std_error=1e-7,
                msd_slope=6.0,
                time_range=(10.0, 90.0),
                temperature=0,
                species='Li',
                n_data_points=100,
                r_squared=0.95
            )


class TestDiffusivityAnalyzer:
    """Test DiffusivityAnalyzer class."""
    
    def test_analyzer_initialization(self):
        """Test basic analyzer initialization."""
        analyzer = DiffusivityAnalyzer(time_step=1.0, step_skip=500)
        
        assert analyzer.time_step == 1.0
        assert analyzer.step_skip == 500
        assert len(analyzer.structures) == 0
    
    def test_fallback_calculation_no_trajectory(self):
        """Test that calculation fails without loaded trajectory."""
        analyzer = DiffusivityAnalyzer()
        
        with pytest.raises(ValueError, match="No trajectory data loaded"):
            analyzer.calculate_diffusivity()
    
    @patch('garnet_analysis.core.diffusivity.HAS_PYMATGEN', False)
    def test_lammps_loading_without_pymatgen(self):
        """Test that LAMMPS loading fails without PyMatGen."""
        with pytest.raises(ImportError, match="PyMatGen required"):
            DiffusivityAnalyzer.from_lammps_dump("dummy.dump", "dummy.data")
    
    def test_element_msd_data_creation(self):
        """Test ElementMSDData creation."""
        data = ElementMSDData(
            element='La',
            diffusivity=1e-10,
            msd_start=0.1,
            msd_middle=1.5,
            msd_final=3.2
        )
        
        assert data.element == 'La'
        assert data.diffusivity == 1e-10
        assert data.msd_final == 3.2


class TestFallbackDiffusivity:
    """Test fallback diffusivity calculation methods."""
    
    def test_fallback_with_mock_structures(self):
        """Test fallback calculation with mock structures."""
        # Create mock structures
        mock_site = Mock()
        mock_site.species_string = 'Li'
        mock_site.coords = np.array([0.0, 0.0, 0.0])
        
        mock_structure = Mock()
        mock_structure.__iter__ = Mock(return_value=iter([mock_site] * 5))
        mock_structure.__len__ = Mock(return_value=5)
        
        analyzer = DiffusivityAnalyzer()
        analyzer.structures = [mock_structure] * 50
        analyzer.time_points = list(range(50))
        
        # Mock the AIMD package to be unavailable
        with patch('garnet_analysis.core.diffusivity.HAS_AIMD', False):
            # This would normally work with proper mock structure setup
            # For now, just test that the method exists and handles the Li species check
            with pytest.raises(NotImplementedError, match="Fallback method only supports Li analysis"):
                analyzer._calculate_fallback('La', 1000)


if __name__ == '__main__':
    pytest.main([__file__])