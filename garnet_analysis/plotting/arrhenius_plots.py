"""
Plotting functions for Arrhenius analysis and temperature-dependent studies.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Any
import pandas as pd

from ..core.arrhenius import ArrheniusAnalyzer, ArrheniusResult
from ..utils.constants import PLOT_PARAMS, PLOT_COLOR_SCHEMES


def plot_arrhenius(analyzers: Dict[str, ArrheniusAnalyzer],
                  reference_analyzers: Optional[Dict[str, ArrheniusAnalyzer]] = None,
                  extend_to_298k: bool = True,
                  show_error_bars: bool = True,
                  color_scheme: str = 'default',
                  title: str = "Arrhenius Plot",
                  save_path: Optional[str] = None,
                  figsize: Tuple[float, float] = (8, 8)) -> plt.Figure:
    """
    Plot Arrhenius analysis for multiple materials.
    
    Args:
        analyzers: Dictionary of {name: ArrheniusAnalyzer} for custom materials
        reference_analyzers: Dictionary of reference materials (tLLZO, cLLZO)
        extend_to_298k: Whether to extend fit lines to 298K
        show_error_bars: Whether to show experimental error bars
        color_scheme: Color scheme name from PLOT_COLOR_SCHEMES
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = PLOT_COLOR_SCHEMES[color_scheme]
    
    # Plot reference materials first
    if reference_analyzers:
        for name, analyzer in reference_analyzers.items():
            if analyzer.fit_result is None:
                continue
                
            # Get data points
            x_data, y_data, y_errors = analyzer.get_data_arrays()
            x_data_inv_temp = x_data  # Already 1000/T
            
            # Plot data points with error bars
            if show_error_bars and np.any(y_errors > 0):
                ax.errorbar(x_data_inv_temp, y_data, yerr=y_errors,
                           fmt='o', color=colors.get(name, 'blue'), 
                           markersize=PLOT_PARAMS['marker_size'],
                           markeredgecolor='black',
                           capsize=PLOT_PARAMS['cap_size'],
                           capthick=2, elinewidth=2,
                           label=name, zorder=3)
            else:
                ax.scatter(x_data_inv_temp, y_data, 
                          c=colors.get(name, 'blue'),
                          s=PLOT_PARAMS['marker_size']**2,
                          edgecolor='black', 
                          label=name, zorder=3)
            
            # Plot fit line
            x_line, y_line = analyzer.get_fit_line_data(extend_to_298k=extend_to_298k)
            y_line_actual = 10**y_line  # Convert from log10(D) to D
            
            ax.plot(x_line, y_line_actual,
                   color=colors.get(name, 'blue'),
                   linewidth=PLOT_PARAMS['line_width'],
                   zorder=2)
            
            # Plot 298K extrapolation point
            if extend_to_298k:
                x_298k = 1000 / 298
                d_298k, d_min, d_max = analyzer.extrapolate_to_temperature(298)
                
                ax.errorbar(x_298k, d_298k, 
                           yerr=[[d_298k - d_min], [d_max - d_298k]],
                           fmt='o', color=colors.get(name, 'blue'),
                           markersize=PLOT_PARAMS['marker_size'] + 2,
                           markeredgecolor='black',
                           capsize=PLOT_PARAMS['cap_size'] + 2,
                           capthick=2, elinewidth=2, zorder=4)
    
    # Plot custom materials
    for name, analyzer in analyzers.items():
        if analyzer.fit_result is None:
            continue
            
        # Determine color based on screening criteria or default
        color = _get_material_color(name, analyzer)
        
        # Get data points
        x_data, y_data, y_errors = analyzer.get_data_arrays()
        x_data_inv_temp = x_data  # Already 1000/T
        
        # Plot fit line
        x_line, y_line = analyzer.get_fit_line_data(extend_to_298k=extend_to_298k)
        y_line_actual = 10**y_line
        
        ax.plot(x_line, y_line_actual, '--', color=color,
               linewidth=0.8, alpha=0.6, zorder=1)
        
        # Plot data points
        ax.scatter(x_data_inv_temp, y_data, c=color,
                  s=PLOT_PARAMS['marker_size']**2 - 10,
                  alpha=0.6, zorder=1)
        
        # Plot error bars if requested
        if show_error_bars and np.any(y_errors > 0):
            ax.errorbar(x_data_inv_temp, y_data, yerr=y_errors,
                       fmt='none', ecolor=color, alpha=0.4,
                       capsize=PLOT_PARAMS['cap_size'],
                       capthick=1, elinewidth=1, zorder=1)
        
        # Plot 298K extrapolation point
        if extend_to_298k:
            x_298k = 1000 / 298
            d_298k, d_min, d_max = analyzer.extrapolate_to_temperature(298)
            
            ax.errorbar(x_298k, d_298k,
                       yerr=[[d_298k - d_min], [d_max - d_298k]],
                       fmt='o', color=color, alpha=0.6,
                       markersize=PLOT_PARAMS['marker_size'],
                       capsize=PLOT_PARAMS['cap_size'],
                       capthick=1, elinewidth=1, zorder=1)
    
    # Formatting
    ax.set_yscale('log')
    ax.set_xlabel("1000/T (K⁻¹)", fontsize=PLOT_PARAMS['label_size'])
    ax.set_ylabel("D (cm²/s)", fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    ax.grid(False)
    
    if reference_analyzers:
        ax.legend(loc='lower left', fontsize=PLOT_PARAMS['legend_size'])
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def plot_activation_energy_distribution(analyzers: Dict[str, ArrheniusAnalyzer],
                                       bins: int = 20,
                                       title: str = "Activation Energy Distribution",
                                       save_path: Optional[str] = None,
                                       figsize: Tuple[float, float] = (8, 6)) -> plt.Figure:
    """
    Plot histogram of activation energies.
    
    Args:
        analyzers: Dictionary of ArrheniusAnalyzer objects
        bins: Number of histogram bins
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Extract activation energies
    activation_energies = []
    for name, analyzer in analyzers.items():
        if analyzer.fit_result is not None:
            activation_energies.append(analyzer.fit_result.activation_energy)
    
    if not activation_energies:
        ax.text(0.5, 0.5, 'No valid activation energy data', 
               ha='center', va='center', transform=ax.transAxes)
        return fig
    
    activation_energies = np.array(activation_energies)
    
    # Create histogram
    ax.hist(activation_energies, bins=bins, alpha=0.7, 
           color='skyblue', edgecolor='black')
    
    # Add statistics
    mean_ea = np.mean(activation_energies)
    std_ea = np.std(activation_energies)
    
    ax.axvline(mean_ea, color='red', linestyle='--', linewidth=2,
              label=f'Mean: {mean_ea:.3f} eV')
    ax.fill_between([mean_ea - std_ea, mean_ea + std_ea], 
                    0, ax.get_ylim()[1], alpha=0.2, color='red',
                    label=f'±1σ: {std_ea:.3f} eV')
    
    ax.set_xlabel('Activation Energy (eV)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_ylabel('Count', fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def plot_298k_extrapolation_comparison(analyzers: Dict[str, ArrheniusAnalyzer],
                                     reference_values: Optional[Dict[str, float]] = None,
                                     error_bars: bool = True,
                                     title: str = "298K Diffusivity Comparison",
                                     save_path: Optional[str] = None,
                                     figsize: Tuple[float, float] = (10, 6)) -> plt.Figure:
    """
    Compare 298K extrapolated diffusivities across materials.
    
    Args:
        analyzers: Dictionary of ArrheniusAnalyzer objects
        reference_values: Reference diffusivity values for comparison
        error_bars: Whether to show error bars
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Extract 298K diffusivities
    materials = []
    diffusivities = []
    lower_bounds = []
    upper_bounds = []
    
    for name, analyzer in analyzers.items():
        if analyzer.fit_result is None:
            continue
            
        d_298k, d_min, d_max = analyzer.extrapolate_to_temperature(298)
        
        materials.append(name)
        diffusivities.append(d_298k)
        lower_bounds.append(d_298k - d_min)
        upper_bounds.append(d_max - d_298k)
    
    if not materials:
        ax.text(0.5, 0.5, 'No valid 298K extrapolation data',
               ha='center', va='center', transform=ax.transAxes)
        return fig
    
    # Sort by diffusivity
    sorted_indices = np.argsort(diffusivities)[::-1]  # Descending order
    materials = [materials[i] for i in sorted_indices]
    diffusivities = [diffusivities[i] for i in sorted_indices]
    lower_bounds = [lower_bounds[i] for i in sorted_indices]
    upper_bounds = [upper_bounds[i] for i in sorted_indices]
    
    x_positions = np.arange(len(materials))
    
    # Plot bars with error bars
    if error_bars:
        ax.errorbar(x_positions, diffusivities,
                   yerr=[lower_bounds, upper_bounds],
                   fmt='o', markersize=PLOT_PARAMS['marker_size'],
                   capsize=PLOT_PARAMS['cap_size'],
                   color='steelblue', markeredgecolor='black')
    else:
        ax.scatter(x_positions, diffusivities,
                  s=PLOT_PARAMS['marker_size']**2,
                  color='steelblue', edgecolor='black')
    
    # Add reference lines if provided
    if reference_values:
        for ref_name, ref_value in reference_values.items():
            ax.axhline(y=ref_value, linestyle='--', 
                      label=f'{ref_name}: {ref_value:.2e}')
    
    ax.set_yscale('log')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(materials, rotation=45, ha='right')
    ax.set_ylabel('298K Diffusivity (cm²/s)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    
    if reference_values:
        ax.legend()
    
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def plot_temperature_dependence(analyzer: ArrheniusAnalyzer,
                               temperature_range: Tuple[float, float] = (250, 1200),
                               highlight_temperatures: Optional[List[float]] = None,
                               title: str = "Temperature Dependence",
                               save_path: Optional[str] = None,
                               figsize: Tuple[float, float] = (8, 6)) -> plt.Figure:
    """
    Plot diffusivity vs temperature for a single material.
    
    Args:
        analyzer: ArrheniusAnalyzer object with fitted data
        temperature_range: Temperature range to plot (K)
        highlight_temperatures: Specific temperatures to highlight
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    if analyzer.fit_result is None:
        raise ValueError("Analyzer must be fitted before plotting")
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Generate temperature array
    temps = np.linspace(temperature_range[0], temperature_range[1], 200)
    diffusivities = []
    
    for temp in temps:
        d, _, _ = analyzer.extrapolate_to_temperature(temp)
        diffusivities.append(d)
    
    diffusivities = np.array(diffusivities)
    
    # Plot main curve
    ax.plot(temps, diffusivities, 'b-', linewidth=PLOT_PARAMS['line_width'],
           label='Arrhenius fit')
    
    # Plot experimental data points
    x_data, y_data, y_errors = analyzer.get_data_arrays()
    exp_temps = 1000.0 / x_data  # Convert back to temperature
    
    ax.errorbar(exp_temps, y_data, yerr=y_errors,
               fmt='ro', markersize=PLOT_PARAMS['marker_size'],
               capsize=PLOT_PARAMS['cap_size'],
               label='Experimental data')
    
    # Highlight specific temperatures
    if highlight_temperatures:
        for temp in highlight_temperatures:
            if temperature_range[0] <= temp <= temperature_range[1]:
                d, d_min, d_max = analyzer.extrapolate_to_temperature(temp)
                ax.errorbar(temp, d, yerr=[[d - d_min], [d_max - d]],
                           fmt='gs', markersize=PLOT_PARAMS['marker_size'] + 2,
                           capsize=PLOT_PARAMS['cap_size'] + 2,
                           label=f'{temp}K extrapolation')
    
    ax.set_yscale('log')
    ax.set_xlabel('Temperature (K)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_ylabel('Diffusivity (cm²/s)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def _get_material_color(name: str, analyzer: ArrheniusAnalyzer) -> str:
    """Determine color for material based on screening criteria."""
    # Simple heuristic - could be made more sophisticated
    if analyzer.fit_result is None:
        return 'gray'
    
    # Color based on 298K diffusivity relative to tLLZO
    d_298k, _, _ = analyzer.extrapolate_to_temperature(298)
    tllzo_298k = 2.089e-11  # Reference value
    
    if d_298k >= tllzo_298k:
        return PLOT_COLOR_SCHEMES['default']['screened_pass']
    else:
        return PLOT_COLOR_SCHEMES['default']['screened_fail']