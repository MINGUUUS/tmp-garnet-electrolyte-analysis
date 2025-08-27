"""
Physical constants and reference data for garnet electrolyte analysis.
"""

import numpy as np
from typing import Dict, Any

# Physical constants
BOLTZMANN_CONSTANT_EV = 8.617333262145e-5  # eV/K
AVOGADRO_NUMBER = 6.02214076e23  # mol⁻¹
ELEMENTARY_CHARGE = 1.602176634e-19  # C
GAS_CONSTANT = 8.314462618  # J/(mol·K)

# Unit conversion factors
ANGSTROM_TO_CM = 1e-8
PS_TO_S = 1e-12
EV_TO_J = 1.602176634e-19

# Temperature conversion factors
TEMPERATURE_CONVERSION_FACTORS = {
    'kelvin_to_celsius': lambda K: K - 273.15,
    'celsius_to_kelvin': lambda C: C + 273.15,
    'inverse_temperature_factor': 1000.0  # For 1000/T in Arrhenius plots
}

# Reference diffusivity data (cm²/s)
LLZO_REFERENCE_DIFFUSIVITY = {
    'tetragonal': {
        1000: {'value': 1.448e-06, 'std': 1.335e-07},
        900:  {'value': 7.694e-07, 'std': 9.250e-08},
        800:  {'value': 4.525e-07, 'std': 6.393e-08},
        700:  {'value': 1.346e-07, 'std': 2.935e-08},
        600:  {'value': 8.079e-08, 'std': 2.192e-08}
    },
    'cubic': {
        1000: {'value': 8.463e-06, 'std': 5.399e-07},
        900:  {'value': 5.244e-06, 'std': 3.601e-07},
        800:  {'value': 3.690e-06, 'std': 2.764e-07},
        700:  {'value': 3.683e-06, 'std': 2.776e-07},
        600:  {'value': 1.558e-06, 'std': 1.462e-07},
        500:  {'value': 7.688e-07, 'std': 9.115e-08}
    }
}

# Convenience constants for backward compatibility
TLLZO_REFERENCE_DIFFUSIVITY = LLZO_REFERENCE_DIFFUSIVITY['tetragonal']
CLLZO_REFERENCE_DIFFUSIVITY = LLZO_REFERENCE_DIFFUSIVITY['cubic']

# Extrapolated 298K values (calculated from Arrhenius fit)
LLZO_298K_DIFFUSIVITY = {
    'tetragonal': {
        'value': 2.089e-11,  # cm²/s
        'lower_bound': 1.123e-11,
        'upper_bound': 3.885e-11
    },
    'cubic': {
        'value': 3.635e-08,  # cm²/s 
        'lower_bound': 1.879e-08,
        'upper_bound': 7.025e-08
    }
}

# Analysis parameters
DEFAULT_ANALYSIS_PARAMS = {
    'diffusivity': {
        'time_step': 1.0,  # fs
        'step_skip': 500,
        'li_analysis_params': {
            'lower_bound': 4.5,  # Å²
            'upper_bound': 0.99,  # fraction
            'minimum_msd_diff': 4.5,  # Å²
            'total_sim_time_limit': 100  # ps
        },
        'other_element_params': {
            'lower_bound': 0.1,  # Å²
            'upper_bound': 0.99,  # fraction
            'minimum_msd_diff': 0.01,  # Å²
            'total_sim_time_limit': 0.001  # ps
        }
    },
    'screening': {
        'melting_threshold': 20.0,  # Å² for MSD
        'melting_diffusivity_threshold': 6.5e-06,  # cm²/s
        'min_diffusivity_threshold': 1e-7,  # cm²/s (failed calculations)
        'llzo_fraction_threshold': 0.9  # 90% of LLZO diffusivity
    }
}

# Atomic masses for common garnet elements (amu)
ATOMIC_MASSES = {
    'Li': 6.941, 'La': 138.906, 'Zr': 91.224, 'O': 15.999,
    'Yb': 173.045, 'Al': 26.982, 'Ga': 69.723, 'Ta': 180.948,
    'W': 183.84, 'Fe': 55.845, 'Co': 58.933, 'Mn': 54.938,
    'Cr': 51.996, 'V': 50.942, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'As': 74.922, 'Se': 78.971, 'Te': 127.60,
    'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
    'Au': 196.967, 'Pd': 106.42, 'Ru': 101.07, 'Tc': 98.0,
    'Mo': 95.95, 'Nb': 92.906, 'Ti': 47.867, 'Sc': 44.956,
    'Y': 88.906, 'Lu': 174.967, 'Hf': 178.49, 'Ce': 140.116,
    'Pr': 140.908, 'Nd': 144.242, 'Pm': 145.0, 'Sm': 150.36,
    'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.925, 'Dy': 162.500,
    'Ho': 164.930, 'Er': 167.259, 'Tm': 168.934, 'Ge': 72.630,
    'Sn': 118.710, 'Pb': 207.2, 'Bi': 208.980, 'Ac': 227.0,
    'Th': 232.038, 'Pa': 231.036, 'U': 238.029, 'B': 10.811,
    'C': 12.011, 'N': 14.007, 'F': 18.998, 'Cl': 35.45,
    'Br': 79.904, 'I': 126.904, 'In': 114.818, 'Tl': 204.383
}

# Common garnet compositions and templates
GARNET_TEMPLATES = {
    'LLZO_type': {
        'formula': 'Li7La3Zr2O12',
        'space_group': 'I-43d',  # cubic
        'prototype': 'garnet'
    },
    'A7B3C2D12_type': {
        'formula': 'A7B3C2D12',
        'description': 'General garnet formula where A=Li, B=La/RE, C=Zr/transition metals, D=O'
    }
}

# Color schemes for plotting
PLOT_COLOR_SCHEMES = {
    'default': {
        'tLLZO': 'green',
        'cLLZO': 'blue', 
        'screened_pass': 'red',
        'screened_fail': 'lightgray',
        'melted': 'orange',
        'failed': 'black'
    },
    'templates': {
        'LiLaZrO': '#1f77b4',
        'LiCrSiO': '#ff7f0e',
        'LiVSiO': '#2ca02c',
        'LiCoMnO': '#d62728',
        'HZnAsO': '#9467bd'
    }
}

# Plotting parameters
PLOT_PARAMS = {
    'figure_size': (8, 6),
    'dpi': 300,
    'font_size': 12,
    'title_size': 14,
    'label_size': 12,
    'legend_size': 10,
    'marker_size': 6,
    'line_width': 1.5,
    'error_bar_width': 1,
    'cap_size': 3
}

# File format specifications
FILE_FORMATS = {
    'lammps_data': {
        'extension': '.data',
        'atom_style': 'atomic'
    },
    'lammps_dump': {
        'extensions': ['.dump', '.traj'],
        'timestep_pattern': r'TIMESTEP\n(\d+)'
    },
    'cif': {
        'extension': '.cif',
        'required_fields': ['_cell_length_a', '_atom_site_label']
    },
    'poscar': {
        'patterns': ['POSCAR*', 'CONTCAR*'],
        'coordinate_types': ['Direct', 'Cartesian']
    }
}

def get_element_mass(element_symbol: str) -> float:
    """
    Get atomic mass for element symbol.
    
    Args:
        element_symbol: Chemical element symbol (e.g., 'Li', 'La')
        
    Returns:
        Atomic mass in amu
        
    Raises:
        KeyError: If element not found in database
    """
    if element_symbol in ATOMIC_MASSES:
        return ATOMIC_MASSES[element_symbol]
    else:
        raise KeyError(f"Atomic mass not found for element: {element_symbol}")


def get_reference_diffusivity(phase: str = 'tetragonal', temperature: float = None) -> Dict[str, Any]:
    """
    Get reference LLZO diffusivity data.
    
    Args:
        phase: 'tetragonal' or 'cubic'
        temperature: Specific temperature (K). If None, returns all temperatures.
        
    Returns:
        Dictionary with diffusivity data
    """
    if phase not in LLZO_REFERENCE_DIFFUSIVITY:
        raise ValueError(f"Unknown phase: {phase}. Use 'tetragonal' or 'cubic'")
        
    data = LLZO_REFERENCE_DIFFUSIVITY[phase]
    
    if temperature is not None:
        if temperature in data:
            return {temperature: data[temperature]}
        else:
            raise KeyError(f"No reference data for {temperature}K in {phase} LLZO")
    
    return data


def calculate_conversion_factor_diffusivity_to_conductivity(
    n_carriers: float, volume_cm3: float, charge: float, temperature: float
) -> float:
    """
    Calculate conversion factor from diffusivity (cm²/s) to conductivity (mS/cm).
    
    Uses Nernst-Einstein equation: σ = (z²F²c)/(RT) * D
    where c = n/(V*N_A) is the concentration
    
    Args:
        n_carriers: Number of charge carriers
        volume_cm3: Unit cell volume in cm³
        charge: Charge state of carrier (e.g., 1 for Li⁺)
        temperature: Temperature in K
        
    Returns:
        Conversion factor (multiply diffusivity by this to get conductivity)
    """
    concentration = n_carriers / (volume_cm3 * AVOGADRO_NUMBER)  # mol/cm³
    
    factor = (1000 * concentration * charge**2 * 
             (AVOGADRO_NUMBER * ELEMENTARY_CHARGE)**2 / 
             (GAS_CONSTANT * temperature))
    
    return factor