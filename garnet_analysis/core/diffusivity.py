"""
Diffusivity analysis module for garnet electrolyte screening.

This module provides classes and functions for calculating ionic diffusivity
from molecular dynamics trajectories, with specific optimizations for 
Li-ion conducting garnets.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
import warnings

try:
    from pymatgen.core.structure import Structure
    from pymatgen.io.lammps.outputs import parse_lammps_dumps
    from pymatgen.io.lammps.data import LammpsData
    from pymatgen.core.periodic_table import Element, Specie
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False
    warnings.warn("PyMatGen not available. Some functionality will be limited.")

try:
    import aimd
    from aimd.diffusion import DiffusivityAnalyzer as AIMDAnalyzer
    from aimd.diffusion import ErrorAnalysisFromDiffusivityAnalyzer
    HAS_AIMD = True
except ImportError:
    HAS_AIMD = False
    warnings.warn("AIMD package not available. Will use fallback methods.")


@dataclass
class DiffusivityResult:
    """Container for diffusivity calculation results."""
    diffusivity: float  # cm²/s
    std_error: float    # cm²/s
    msd_slope: float    # Å²/ps
    time_range: Tuple[float, float]  # ps
    temperature: float  # K
    species: str
    n_data_points: int
    r_squared: float
    
    def __post_init__(self):
        """Validate result data."""
        if self.diffusivity < 0:
            raise ValueError("Diffusivity cannot be negative")
        if self.temperature <= 0:
            raise ValueError("Temperature must be positive")


@dataclass
class ElementMSDData:
    """Container for element-specific MSD data."""
    element: str
    diffusivity: float  # cm²/s
    msd_start: float    # Å²
    msd_middle: float   # Å²
    msd_final: float    # Å²


class DiffusivityAnalyzer:
    """
    Analyze ionic diffusivity from molecular dynamics trajectories.
    
    This class provides methods to calculate diffusivity from LAMMPS dumps,
    POSCAR trajectories, or pre-calculated MSD data.
    """
    
    def __init__(self, time_step: float = 1.0, step_skip: int = 500):
        """
        Initialize the analyzer.
        
        Args:
            time_step: Simulation time step in fs
            step_skip: Interval between saved trajectory frames
        """
        self.time_step = time_step
        self.step_skip = step_skip
        self.structures = []
        self.time_points = []
        
    @classmethod
    def from_lammps_dump(cls, dump_path: str, data_path: str, 
                        time_step: float = 1.0, step_skip: int = 500,
                        supercell: bool = False) -> 'DiffusivityAnalyzer':
        """
        Create analyzer from LAMMPS dump file.
        
        Args:
            dump_path: Path to LAMMPS dump file
            data_path: Path to LAMMPS data file
            time_step: Time step in fs
            step_skip: Frame interval
            supercell: Whether to create 2x2x2 supercell
            
        Returns:
            Initialized DiffusivityAnalyzer
        """
        if not HAS_PYMATGEN:
            raise ImportError("PyMatGen required for LAMMPS dump parsing")
            
        analyzer = cls(time_step=time_step, step_skip=step_skip)
        analyzer._load_lammps_trajectory(dump_path, data_path, supercell)
        return analyzer
        
    @classmethod 
    def from_poscar_directory(cls, poscar_dir: str, time_step: float = 1.0,
                            step_skip: int = 500) -> 'DiffusivityAnalyzer':
        """
        Create analyzer from directory of POSCAR files.
        
        Args:
            poscar_dir: Directory containing POSCAR-* files
            time_step: Time step in fs
            step_skip: Frame interval
            
        Returns:
            Initialized DiffusivityAnalyzer
        """
        if not HAS_PYMATGEN:
            raise ImportError("PyMatGen required for POSCAR parsing")
            
        analyzer = cls(time_step=time_step, step_skip=step_skip)
        analyzer._load_poscar_trajectory(poscar_dir)
        return analyzer
        
    def _load_lammps_trajectory(self, dump_path: str, data_path: str, 
                               supercell: bool = False) -> None:
        """Load trajectory from LAMMPS dump file."""
        # Load LAMMPS data file for structure template
        ld = LammpsData.from_file(data_path, atom_style='atomic')
        if supercell:
            template_structure = ld.structure.make_supercell([2, 2, 2])
        else:
            template_structure = ld.structure
            
        # Parse dump file
        dumps = list(parse_lammps_dumps(dump_path))
        
        # Convert dumps to structures
        self.structures = []
        self.time_points = []
        
        for dump in dumps:
            structure = self._dump_to_structure(dump, template_structure, ld.masses)
            self.structures.append(structure)
            self.time_points.append(dump.timestep * self.time_step / 1000)  # Convert to ps
            
    def _dump_to_structure(self, dump, template_structure: Structure, 
                          masses: np.ndarray) -> Structure:
        """Convert LAMMPS dump to PyMatGen Structure."""
        # Get element symbols from masses
        def get_element_by_mass(mass: float, tolerance: float = 0.1) -> str:
            for el in Element:
                if abs(el.atomic_mass - mass) <= tolerance:
                    return el.symbol
            return "Unknown"
            
        atomic_masses = np.array(masses).flatten().tolist()
        element_order = [get_element_by_mass(mass) for mass in atomic_masses]
        

        
        
        # Create new structure with updated coordinates
        new_structure = template_structure.copy()
        new_structure.sort(key=lambda site: element_order.index(site.species_string))
        
        
        # Sort dump data by atom type
        dump.data = dump.data.sort_values('type')
        
        if len(new_structure) != len(dump.data):
            raise ValueError("Structure and dump have different number of atoms")
            
        
        # Update coordinates
        new_coords = np.array(dump.data[['x', 'y', 'z']])
        for i, site in enumerate(new_structure):
            site.coords = new_coords[i]
            
        return new_structure
        
    def _load_poscar_trajectory(self, poscar_dir: str) -> None:
        """Load trajectory from directory of POSCAR files."""
        import os
        from pathlib import Path
        
        poscar_files = sorted([f for f in os.listdir(poscar_dir) if f.startswith('POSCAR-')])
        
        self.structures = []
        self.time_points = []
        
        for poscar_file in poscar_files:
            timestep = int(poscar_file.split('-')[1])
            structure = Structure.from_file(Path(poscar_dir) / poscar_file)
            
            self.structures.append(structure)
            self.time_points.append(timestep * self.time_step / 1000)  # Convert to ps
            
    def calculate_diffusivity(self, species: str = 'Li', temperature: float = 1000,
                            lower_bound: float = 4.5, upper_bound: float = 0.99,
                            min_msd_diff: float = 4.5, max_sim_time: float = 100) -> DiffusivityResult:
        """
        Calculate Li-ion diffusivity from loaded trajectory.
        
        Args:
            species: Ion species to analyze (default 'Li')
            temperature: Simulation temperature in K
            lower_bound: Lower bound for linear fit region (Å²) 
            upper_bound: Upper bound fraction for linear fit
            min_msd_diff: Minimum MSD difference for valid fit
            max_sim_time: Maximum simulation time to consider
            
        Returns:
            DiffusivityResult containing calculated diffusivity and metadata
        """
        if not self.structures:
            raise ValueError("No trajectory data loaded")
            
        if HAS_AIMD:
            return self._calculate_with_aimd(species, temperature, lower_bound, 
                                           upper_bound, min_msd_diff, max_sim_time)
        else:
            return self._calculate_fallback(species, temperature)
            
    def _calculate_with_aimd(self, species: str, temperature: float,
                           lower_bound: float, upper_bound: float,
                           min_msd_diff: float, max_sim_time: float) -> DiffusivityResult:
        """Calculate diffusivity using AIMD package."""
        # Set up species
        if species == 'Li':
            specie = Specie.from_str('Li+')
        else:
            specie = Specie(species, 0)
        # Create AIMD analyzer
        difs = AIMDAnalyzer.from_structures(
            structures=self.structures,
            specie=str(specie.element),
            temperature=temperature,
            time_step=self.time_step,
            step_skip=self.step_skip,
            time_intervals_number=len(self.structures),
            spec_dict={
                'lower_bound': lower_bound,
                'upper_bound': upper_bound,
                'minimum_msd_diff': min_msd_diff,
                'total_sim_time_limit': max_sim_time
            }
        )
        
        # Error analysis
        ea = ErrorAnalysisFromDiffusivityAnalyzer(difs)
        if species == 'Li':
            summary = ea.get_summary_dict(oxidized_specie='Li+')
        else:
            summary = ea.get_summary_dict(oxidized_specie=species)
        # Calculate time range
        time_ps = np.array(difs.dt) / 1000
        time_range = (
            float(time_ps[difs.lower_bound_index-1]),
            float(time_ps[difs.upper_bound_index])
        )
        
        # Calculate R²
        fit_indices = slice(difs.lower_bound_index, difs.upper_bound_index)
        x_fit = time_ps[fit_indices] 
        y_fit = difs.msd[fit_indices]
        
        if len(x_fit) > 1:
            coeffs = np.polyfit(x_fit, y_fit, 1)
            y_pred = np.polyval(coeffs, x_fit)
            ss_res = np.sum((y_fit - y_pred) ** 2)
            ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        else:
            r_squared = 0
            
        return DiffusivityResult(
            diffusivity=difs.diffusivity,
            std_error=summary['diffusivity_standard_deviation'],
            msd_slope=coeffs[0] if len(x_fit) > 1 else 0,
            time_range=time_range,
            temperature=temperature,
            species=species,
            n_data_points=len(x_fit),
            r_squared=r_squared
        )
        
    def _calculate_fallback(self, species: str, temperature: float) -> DiffusivityResult:
        """Fallback diffusivity calculation without AIMD package."""
        # Simple MSD calculation
        if species != 'Li':
            raise NotImplementedError("Fallback method only supports Li analysis")
        
        # Extract Li positions over time
        li_trajectories = []
        for structure in self.structures:
            li_positions = []
            for site in structure:
                if site.species_string == 'Li':
                    li_positions.append(site.coords)
            li_trajectories.append(np.array(li_positions))
        
        # Calculate MSD
        n_li = len(li_trajectories[0])
        n_times = len(li_trajectories)
        msd_values = []
        
        initial_positions = li_trajectories[0]
        for t_idx in range(n_times):
            current_positions = li_trajectories[t_idx]
            displacements = current_positions - initial_positions
            
            # Apply periodic boundary conditions if needed
            # (This is simplified - real implementation would need cell parameters)
            
            msd = np.mean(np.sum(displacements**2, axis=1))
            msd_values.append(msd)
        
        # Fit linear region (simplified)
        time_ps = np.array(self.time_points)
        msd_array = np.array(msd_values)
        
        # Use middle 50% for fit
        start_idx = len(msd_values) // 4
        end_idx = 3 * len(msd_values) // 4
        
        coeffs = np.polyfit(time_ps[start_idx:end_idx], msd_array[start_idx:end_idx], 1)
        slope = coeffs[0]  # Å²/ps
        
        # Convert to diffusivity: D = slope / 6 * 1E-4 (cm²/s)
        diffusivity = slope / 6 * 1e-4
        
        # Estimate error (simplified)
        y_fit = np.polyval(coeffs, time_ps[start_idx:end_idx])
        residuals = msd_array[start_idx:end_idx] - y_fit
        std_error = np.std(residuals) / 6 * 1e-4
        
        # Calculate R²
        y_mean = np.mean(msd_array[start_idx:end_idx])
        ss_tot = np.sum((msd_array[start_idx:end_idx] - y_mean)**2)
        ss_res = np.sum(residuals**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return DiffusivityResult(
            diffusivity=diffusivity,
            std_error=std_error, 
            msd_slope=slope,
            time_range=(time_ps[start_idx], time_ps[end_idx]),
            temperature=temperature,
            species=species,
            n_data_points=end_idx - start_idx,
            r_squared=r_squared
        )
        
    def analyze_all_elements(self, temperature: float = 1000) -> Dict[str, ElementMSDData]:
        """
        Analyze MSD for all non-Li elements in the structure.
        
        Args:
            temperature: Simulation temperature in K
            
        Returns:
            Dictionary mapping element symbols to ElementMSDData
        """
        if not self.structures:
            raise ValueError("No trajectory data loaded")
            
        first_structure = self.structures[0]
        elements = list({site.specie.symbol for site in first_structure if site.specie.symbol != 'Li'})
        
        results = {}
        
        for element in elements:
            if HAS_AIMD:
                
                # Use AIMD for non-Li elements
                specie = Specie(element, 0)
                difs = AIMDAnalyzer.from_structures(
                    structures=self.structures,
                    specie=str(specie.element),
                    temperature=temperature,
                    time_step=self.time_step,
                    step_skip=self.step_skip,
                    time_intervals_number=len(self.structures),
                    spec_dict={
                        'lower_bound': 0.1,
                        'upper_bound': 0.99,
                        'minimum_msd_diff': 0.01,
                        'total_sim_time_limit': 0.001
                    }
                )
                
                
                time_ps = np.array(difs.dt) / 1000
                n_points = len(time_ps)
                
                results[element] = ElementMSDData(
                    element=element,
                    diffusivity=difs.diffusivity,
                    msd_start=difs.msd[0],
                    msd_middle=difs.msd[n_points//2], 
                    msd_final=difs.msd[-1]
                )
            else:
                
                # Fallback method for other elements
                results[element] = ElementMSDData(
                    element=element,
                    diffusivity=0.0,  # Not calculated in fallback
                    msd_start=0.0,
                    msd_middle=0.0,
                    msd_final=0.0
                )
                
        return results