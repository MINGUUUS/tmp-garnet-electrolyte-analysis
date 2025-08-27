"""
Arrhenius analysis module for temperature-dependent diffusivity analysis.

This module provides functionality to fit Arrhenius equations to diffusivity data
and extrapolate to room temperature conditions.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from scipy.optimize import curve_fit
import warnings

try:
    from pymatgen.core.structure import Structure
    from pymatgen.core.periodic_table import Specie
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False


@dataclass
class ArrheniusDataPoint:
    """Single temperature-diffusivity data point."""
    temperature: float  # K
    diffusivity: float  # cm²/s
    std_error: float   # cm²/s
    
    def __post_init__(self):
        """Validate data point."""
        if self.temperature <= 0:
            raise ValueError("Temperature must be positive")
        if self.diffusivity < 0:
            raise ValueError("Diffusivity cannot be negative") 
        if self.std_error < 0:
            raise ValueError("Standard error cannot be negative")


@dataclass 
class ArrheniusResult:
    """Results from Arrhenius fitting."""
    slope: float                    # K
    intercept: float               # log10(D)
    slope_error: float             # K  
    intercept_error: float         # log10(D)
    activation_energy: float       # eV
    activation_energy_error: float # eV
    r_squared: float              # Fit quality
    n_points: int                 # Number of data points
    temperature_range: Tuple[float, float]  # K
    
    def extrapolate_to_temperature(self, temperature: float) -> Tuple[float, float, float]:
        """
        Extrapolate diffusivity to specified temperature.
        
        Args:
            temperature: Target temperature in K
            
        Returns:
            Tuple of (diffusivity, lower_bound, upper_bound) in cm²/s
        """
        x_pred = 1000.0 / temperature
        log_d_pred = self.slope * x_pred + self.intercept
        
        # Error propagation
        log_d_error = np.sqrt((self.slope_error * x_pred)**2 + (self.intercept_error)**2)
        
        d_pred = 10**log_d_pred
        d_min = 10**(log_d_pred - log_d_error)
        d_max = 10**(log_d_pred + log_d_error)
        
        return d_pred, d_min, d_max


class ArrheniusAnalyzer:
    """
    Analyze temperature dependence of ionic diffusivity using Arrhenius equation.
    
    The Arrhenius equation: D = D0 * exp(-Ea / (kB * T))
    In log form: log10(D) = log10(D0) - Ea/(kB * ln(10)) * (1000/T)
    """
    
    def __init__(self):
        """Initialize Arrhenius analyzer."""
        self.data_points: List[ArrheniusDataPoint] = []
        self.fit_result: Optional[ArrheniusResult] = None
        
    def add_data_point(self, temperature: float, diffusivity: float, 
                      std_error: float = 0.0) -> None:
        """
        Add a temperature-diffusivity data point.
        
        Args:
            temperature: Temperature in K
            diffusivity: Diffusivity in cm²/s
            std_error: Standard error in cm²/s
        """
        point = ArrheniusDataPoint(temperature, diffusivity, std_error)
        self.data_points.append(point)
        
    def add_data_from_dict(self, data_dict: Dict[float, Dict[str, float]]) -> None:
        """
        Add data from dictionary format.
        
        Args:
            data_dict: {temperature: {'diff': value, 'std': error}}
        """
        for temp, data in data_dict.items():
            self.add_data_point(temp, data['diff'], data.get('std', 0.0))
            
    def fit(self, temperatures_to_use: Optional[List[float]] = None) -> ArrheniusResult:
        """
        Fit Arrhenius equation to data points.
        
        Args:
            temperatures_to_use: Specific temperatures to include in fit.
                               If None, uses all data points.
                               
        Returns:
            ArrheniusResult containing fit parameters and statistics
        """
        if len(self.data_points) < 2:
            raise ValueError("Need at least 2 data points for Arrhenius fit")
            
        # Filter data points if specific temperatures requested
        if temperatures_to_use is not None:
            points_to_use = [
                point for point in self.data_points 
                if point.temperature in temperatures_to_use
            ]
        else:
            points_to_use = self.data_points
            
        if len(points_to_use) < 2:
            raise ValueError("Not enough valid data points for fit")
            
        # Sort by temperature (highest to lowest for consistency)
        points_to_use.sort(key=lambda p: p.temperature, reverse=True)
        
        # Prepare data for fitting
        x_data = np.array([1000.0 / point.temperature for point in points_to_use])
        y_data = np.array([np.log10(point.diffusivity) for point in points_to_use])
        
        # Calculate weights from standard errors
        y_errors = np.array([
            np.log10(np.e) * point.std_error / point.diffusivity
            if point.std_error > 0 else 1.0
            for point in points_to_use
        ])
        
        # Perform weighted linear fit
        try:
            (slope, intercept), cov_matrix = curve_fit(
                self._linear_function, x_data, y_data,
                sigma=y_errors, 
                absolute_sigma=True
            )
        except Exception as e:
            raise RuntimeError(f"Arrhenius fit failed: {e}")
            
        # Extract parameter errors
        slope_error = np.sqrt(cov_matrix[0, 0])
        intercept_error = np.sqrt(cov_matrix[1, 1])
        
        # Calculate activation energy (eV)
        # Ea = -slope * kB * ln(10) * 1000, where kB = 8.617e-5 eV/K
        boltzmann_ev = 8.617e-5  # eV/K
        ea = -slope * boltzmann_ev * np.log(10) * 1000
        ea_error = slope_error * boltzmann_ev * np.log(10) * 1000
        
        # Calculate R² 
        y_pred = slope * x_data + intercept
        ss_res = np.sum((y_data - y_pred)**2)
        ss_tot = np.sum((y_data - np.mean(y_data))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        # Create result
        self.fit_result = ArrheniusResult(
            slope=slope,
            intercept=intercept,
            slope_error=slope_error,
            intercept_error=intercept_error,
            activation_energy=ea,
            activation_energy_error=ea_error,
            r_squared=r_squared,
            n_points=len(points_to_use),
            temperature_range=(
                min(point.temperature for point in points_to_use),
                max(point.temperature for point in points_to_use)
            )
        )
        
        return self.fit_result
        
    def extrapolate_to_temperature(self, temperature: float) -> Tuple[float, float, float]:
        """
        Extrapolate diffusivity to specified temperature.
        
        Args:
            temperature: Target temperature in K
            
        Returns:
            Tuple of (diffusivity, lower_bound, upper_bound) in cm²/s
        """
        if self.fit_result is None:
            raise ValueError("Must fit data before extrapolation")
            
        return self.fit_result.extrapolate_to_temperature(temperature)
        
    def get_conductivity_at_temperature(self, temperature: float, structure: Structure,
                                      species: str = 'Li+') -> Tuple[float, float, float]:
        """
        Calculate ionic conductivity at specified temperature using Nernst-Einstein equation.
        
        Args:
            temperature: Temperature in K
            structure: Crystal structure 
            species: Ionic species (must include charge, e.g., 'Li+')
            
        Returns:
            Tuple of (conductivity, lower_bound, upper_bound) in mS/cm
        """
        if not HAS_PYMATGEN:
            raise ImportError("PyMatGen required for conductivity calculation")
            
        # Get diffusivity at temperature
        d_pred, d_min, d_max = self.extrapolate_to_temperature(temperature)
        
        # Calculate conversion factor
        conversion_factor = self._get_conductivity_conversion_factor(
            structure, species, temperature
        )
        
        return (
            conversion_factor * d_pred,
            conversion_factor * d_min, 
            conversion_factor * d_max
        )
        
    def _linear_function(self, x: np.ndarray, slope: float, intercept: float) -> np.ndarray:
        """Linear function for curve fitting."""
        return slope * x + intercept
        
    def _get_conductivity_conversion_factor(self, structure: Structure, 
                                          species: str, temperature: float) -> float:
        """
        Calculate conversion factor from diffusivity to conductivity.
        
        Uses Nernst-Einstein equation:
        σ = z²F²c/(RT) * D
        where c = n/(V*NA)
        
        Args:
            structure: Crystal structure
            species: Species string with oxidation (e.g., 'Li+')  
            temperature: Temperature in K
            
        Returns:
            Conversion factor (conductivity = factor * diffusivity)
        """
        try:
            specie = Specie.from_str(species)
        except:
            raise ValueError(f"Invalid species format: {species}")
            
        z = specie.oxi_state  # Charge
        
        # Count number of species in structure
        if isinstance(list(structure.composition.items())[0][0], Specie):
            # Oxidation-decorated structure
            n = structure.composition[species]
        else:
            # Non-decorated structure
            n = structure.composition[str(specie.element)]
            
        if n == 0:
            raise ValueError(f"No {species} found in structure")
            
        # Physical constants
        vol_cm3 = structure.volume * 1e-24  # Convert Å³ to cm³
        N_A = 6.022140857e23  # Avogadro's number
        e = 1.6021766208e-19  # Elementary charge (C)
        R = 8.3144598  # Gas constant (J/mol/K)
        
        # Conversion factor (mS/cm per cm²/s)
        factor = (1000 * n / (vol_cm3 * N_A) * 
                 z**2 * (N_A * e)**2 / (R * temperature))
        
        return factor
        
    def get_fit_line_data(self, n_points: int = 100, 
                         extend_to_298k: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get x,y data for plotting fitted Arrhenius line.
        
        Args:
            n_points: Number of points for smooth line
            extend_to_298k: Whether to extend line to 298K
            
        Returns:
            Tuple of (x_values, y_values) where x = 1000/T and y = log10(D)
        """
        if self.fit_result is None:
            raise ValueError("Must fit data before getting line data")
            
        temp_min, temp_max = self.fit_result.temperature_range
        
        if extend_to_298k:
            temp_min = min(temp_min, 298)
            
        x_min = 1000 / temp_max
        x_max = 1000 / temp_min
        
        x_line = np.linspace(x_min, x_max, n_points)
        y_line = self.fit_result.slope * x_line + self.fit_result.intercept
        
        return x_line, y_line
        
    def get_data_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get experimental data as arrays for plotting.
        
        Returns:
            Tuple of (x_values, y_values, y_errors) where x = 1000/T
        """
        if not self.data_points:
            return np.array([]), np.array([]), np.array([])
            
        x_data = np.array([1000.0 / point.temperature for point in self.data_points])
        y_data = np.array([point.diffusivity for point in self.data_points])
        y_errors = np.array([point.std_error for point in self.data_points])
        
        return x_data, y_data, y_errors


def create_reference_analyzers() -> Dict[str, ArrheniusAnalyzer]:
    """
    Create Arrhenius analyzers for LLZO reference materials.
    
    Returns:
        Dictionary with 'tLLZO' and 'cLLZO' analyzers
    """
    # Reference data from experimental/computational studies
    tllzo_data = {
        1000: {'diff': 1.448e-06, 'std': 1.335e-07},
        900:  {'diff': 7.694e-07, 'std': 9.250e-08},
        800:  {'diff': 4.525e-07, 'std': 6.393e-08},
        700:  {'diff': 1.346e-07, 'std': 2.935e-08},
        600:  {'diff': 8.079e-08, 'std': 2.192e-08}
    }
    
    cllzo_data = {
        1000: {'diff': 8.463e-06, 'std': 5.399e-07},
        900:  {'diff': 5.244e-06, 'std': 3.601e-07},
        800:  {'diff': 3.690e-06, 'std': 2.764e-07},
        700:  {'diff': 3.683e-06, 'std': 2.776e-07},
        600:  {'diff': 1.558e-06, 'std': 1.462e-07},
        500:  {'diff': 7.688e-07, 'std': 9.115e-08}
    }
    
    analyzers = {}
    
    # Create tLLZO analyzer
    tllzo_analyzer = ArrheniusAnalyzer()
    tllzo_analyzer.add_data_from_dict(tllzo_data)
    tllzo_analyzer.fit()
    analyzers['tLLZO'] = tllzo_analyzer
    
    # Create cLLZO analyzer  
    cllzo_analyzer = ArrheniusAnalyzer()
    cllzo_analyzer.add_data_from_dict(cllzo_data)
    cllzo_analyzer.fit()
    analyzers['cLLZO'] = cllzo_analyzer
    
    return analyzers