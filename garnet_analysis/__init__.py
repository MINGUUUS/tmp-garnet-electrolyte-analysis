"""
Garnet Electrolyte Analysis Suite

A comprehensive Python toolkit for computational analysis of garnet-structured solid electrolytes.
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import main classes and functions for easy access
from .core.diffusivity import DiffusivityAnalyzer, DiffusivityResult
from .core.arrhenius import ArrheniusAnalyzer, ArrheniusResult
# from .core.energy import EnergyAnalyzer, FormationEnergyCalculator
# from .core.structure import StructureHandler, CifProcessor

from .plotting.diffusivity_plots import plot_msd, plot_diffusivity_screening
from .plotting.arrhenius_plots import plot_arrhenius, plot_temperature_dependence
# from .plotting.screening_plots import plot_formation_energy_vs_diffusivity

from .utils.constants import (
    LLZO_REFERENCE_DIFFUSIVITY, 
    CLLZO_REFERENCE_DIFFUSIVITY,
    TEMPERATURE_CONVERSION_FACTORS
)

# Define what gets imported with "from garnet_analysis import *"
__all__ = [
    # Core classes
    'DiffusivityAnalyzer', 'DiffusivityResult',
    'ArrheniusAnalyzer', 'ArrheniusResult', 
    'EnergyAnalyzer', 'FormationEnergyCalculator',
    'StructureHandler', 'CifProcessor',
    
    # Plotting functions
    'plot_msd', 'plot_diffusivity_screening',
    'plot_arrhenius', 'plot_temperature_dependence',
    'plot_formation_energy_vs_diffusivity',
    
    # Constants
    'LLZO_REFERENCE_DIFFUSIVITY',
    'CLLZO_REFERENCE_DIFFUSIVITY',
    'TEMPERATURE_CONVERSION_FACTORS'
]