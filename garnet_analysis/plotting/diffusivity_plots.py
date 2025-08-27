"""
Plotting functions for diffusivity analysis and screening results.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Any
import pandas as pd

from ..core.diffusivity import DiffusivityResult, ElementMSDData
from ..utils.constants import PLOT_PARAMS, PLOT_COLOR_SCHEMES, LLZO_298K_DIFFUSIVITY


def plot_msd(time_ps: np.ndarray, msd: np.ndarray, 
             diffusivity_result: Optional[DiffusivityResult] = None,
             other_elements: Optional[Dict[str, ElementMSDData]] = None,
             title: str = "Mean Square Displacement",
             save_path: Optional[str] = None,
             figsize: Tuple[float, float] = (8, 6)) -> plt.Figure:
    """
    Plot mean square displacement with diffusivity fit.
    
    Args:
        time_ps: Time points in ps
        msd: MSD values in Å²
        diffusivity_result: Li diffusivity analysis results
        other_elements: MSD data for other elements
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size (width, height)
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot Li MSD
    ax.scatter(time_ps, msd, c='black', s=5, alpha=0.7, label='Li', zorder=3)
    
    # Plot fit line if available
    if diffusivity_result is not None:
        t_start, t_end = diffusivity_result.time_range
        fit_mask = (time_ps >= t_start) & (time_ps <= t_end)
        
        if np.any(fit_mask):
            t_fit = time_ps[fit_mask]
            # Reconstruct fit line
            slope = diffusivity_result.msd_slope
            intercept = msd[fit_mask][0] - slope * t_fit[0]  # Calculate intercept
            msd_fit = slope * t_fit + intercept
            
            ax.plot(t_fit, msd_fit, 'gray', linewidth=2, alpha=0.8, zorder=2)
            
        # Add diffusivity to title
        title += f"\nD = {diffusivity_result.diffusivity:.3e} cm²/s"
    
    # Plot other elements if available
    if other_elements:
        for element, data in other_elements.items():
            # Simple scatter plot for other elements (would need time series for full plot)
            ax.scatter([], [], marker='o', s=3, alpha=0.3, label=element)
    
    ax.set_xlabel('Time (ps)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_ylabel('MSD (Å²)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def plot_diffusivity_screening(diffusivity_data: pd.DataFrame,
                              formation_energy_col: str = 'formation_energy',
                              diffusivity_col: str = 'diffusivity',
                              std_col: Optional[str] = 'diffusivity_std',
                              color_by: Optional[str] = None,
                              reference_lines: bool = True,
                              screening_criteria: Optional[Dict[str, Any]] = None,
                              title: str = "Diffusivity Screening Results",
                              save_path: Optional[str] = None,
                              figsize: Tuple[float, float] = (10, 6)) -> plt.Figure:
    """
    Plot diffusivity vs formation energy screening results.
    
    Args:
        diffusivity_data: DataFrame with screening results
        formation_energy_col: Column name for formation energy
        diffusivity_col: Column name for diffusivity
        std_col: Column name for standard error (optional)
        color_by: Column to use for color coding
        reference_lines: Whether to show LLZO reference lines
        screening_criteria: Criteria for pass/fail coloring
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Prepare data
    x_data = diffusivity_data[formation_energy_col]
    y_data = diffusivity_data[diffusivity_col]
    
    # Handle error bars
    if std_col and std_col in diffusivity_data.columns:
        y_err = diffusivity_data[std_col]
        ax.errorbar(x_data, y_data, yerr=y_err, 
                   fmt='none', ecolor='gray', alpha=0.5, 
                   capsize=PLOT_PARAMS['cap_size'], zorder=1)
    
    # Color coding
    if color_by and color_by in diffusivity_data.columns:
        scatter = ax.scatter(x_data, y_data, 
                           c=diffusivity_data[color_by], 
                           s=PLOT_PARAMS['marker_size']**2,
                           edgecolor='black', alpha=0.8, zorder=2)
        plt.colorbar(scatter, label=color_by)
    elif screening_criteria:
        # Apply screening criteria coloring
        colors = []
        for _, row in diffusivity_data.iterrows():
            if _passes_screening_criteria(row, screening_criteria):
                colors.append(PLOT_COLOR_SCHEMES['default']['screened_pass'])
            else:
                colors.append(PLOT_COLOR_SCHEMES['default']['screened_fail'])
                
        ax.scatter(x_data, y_data, c=colors, 
                  s=PLOT_PARAMS['marker_size']**2,
                  edgecolor='black', alpha=0.8, zorder=2)
    else:
        ax.scatter(x_data, y_data, 
                  c=PLOT_COLOR_SCHEMES['default']['screened_pass'],
                  s=PLOT_PARAMS['marker_size']**2,
                  edgecolor='black', alpha=0.8, zorder=2)
    
    # Reference lines
    if reference_lines:
        tllzo_298k = LLZO_298K_DIFFUSIVITY['tetragonal']['value']
        cllzo_298k = LLZO_298K_DIFFUSIVITY['cubic']['value']
        
        ax.axhline(y=tllzo_298k, color=PLOT_COLOR_SCHEMES['default']['tLLZO'],
                  linestyle='--', linewidth=PLOT_PARAMS['line_width'],
                  label='tLLZO (298K)', zorder=0)
        ax.axhline(y=cllzo_298k, color=PLOT_COLOR_SCHEMES['default']['cLLZO'],
                  linestyle='-', linewidth=PLOT_PARAMS['line_width'],
                  label='cLLZO (298K)', zorder=0)
    
    ax.set_yscale('log')
    ax.set_xlabel('|Formation Energy| (meV/atom)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_ylabel('Diffusivity (cm²/s)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    
    if reference_lines:
        ax.legend(fontsize=PLOT_PARAMS['legend_size'])
    
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def plot_diffusivity_histogram(diffusivity_data: pd.DataFrame,
                              diffusivity_col: str = 'diffusivity',
                              formation_energy_col: str = 'formation_energy',
                              bins: int = 50,
                              screening_criteria: Optional[Dict[str, Any]] = None,
                              title: str = "Diffusivity Distribution",
                              save_path: Optional[str] = None,
                              figsize: Tuple[float, float] = (10, 8)) -> plt.Figure:
    """
    Plot histogram of diffusivity values with formation energy distribution.
    
    Args:
        diffusivity_data: DataFrame with results
        diffusivity_col: Column name for diffusivity
        formation_energy_col: Column name for formation energy  
        bins: Number of histogram bins
        screening_criteria: Criteria for highlighting passed structures
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    from matplotlib import gridspec
    
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4], hspace=0.1)
    
    # Top panel: Formation energy histogram
    ax_top = plt.subplot(gs[0])
    fe_data = diffusivity_data[formation_energy_col]
    
    bin_width = (fe_data.max() - fe_data.min()) / bins
    fe_bins = np.arange(fe_data.min() - bin_width/2, fe_data.max() + bin_width, bin_width)
    
    ax_top.hist(fe_data, bins=fe_bins, color='lightblue', alpha=0.7, edgecolor='black')
    ax_top.set_xlim(fe_data.min() - 5, fe_data.max() + 5)
    ax_top.set_ylabel('Count')
    ax_top.set_title('Formation Energy Distribution')
    
    # Bottom panel: Diffusivity scatter plot
    ax_main = plt.subplot(gs[1], sharex=ax_top)
    
    x_data = diffusivity_data[formation_energy_col]
    y_data = diffusivity_data[diffusivity_col]
    
    if screening_criteria:
        colors = []
        for _, row in diffusivity_data.iterrows():
            if _passes_screening_criteria(row, screening_criteria):
                colors.append(PLOT_COLOR_SCHEMES['default']['screened_pass'])
            else:
                colors.append(PLOT_COLOR_SCHEMES['default']['screened_fail'])
    else:
        colors = PLOT_COLOR_SCHEMES['default']['screened_pass']
    
    ax_main.scatter(x_data, y_data, c=colors, 
                   s=PLOT_PARAMS['marker_size']**2,
                   edgecolor='black', alpha=0.8)
    
    ax_main.set_yscale('log')
    ax_main.set_xlabel('|Formation Energy| (meV/atom)', fontsize=PLOT_PARAMS['label_size'])
    ax_main.set_ylabel('Diffusivity (cm²/s)', fontsize=PLOT_PARAMS['label_size'])
    ax_main.grid(True, axis='y', alpha=0.3)
    
    plt.setp(ax_top.get_xticklabels(), visible=False)
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def plot_template_comparison(diffusivity_data: pd.DataFrame,
                           template_col: str = 'template',
                           diffusivity_col: str = 'diffusivity',
                           std_col: Optional[str] = 'diffusivity_std',
                           title: str = "Diffusivity by Template",
                           save_path: Optional[str] = None,
                           figsize: Tuple[float, float] = (10, 6)) -> plt.Figure:
    """
    Compare diffusivity distributions across different templates.
    
    Args:
        diffusivity_data: DataFrame with results
        template_col: Column name for template classification
        diffusivity_col: Column name for diffusivity
        std_col: Column name for standard error
        title: Plot title
        save_path: Path to save figure
        figsize: Figure size
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    templates = diffusivity_data[template_col].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(templates)))
    
    for i, template in enumerate(templates):
        template_data = diffusivity_data[diffusivity_data[template_col] == template]
        
        x_data = np.full(len(template_data), i)  # x position for template
        y_data = template_data[diffusivity_col]
        
        # Add some jitter to x positions
        x_jitter = x_data + np.random.normal(0, 0.1, len(x_data))
        
        ax.scatter(x_jitter, y_data, c=[colors[i]], 
                  s=PLOT_PARAMS['marker_size']**2,
                  alpha=0.7, edgecolor='black',
                  label=template)
        
        # Add error bars if available
        if std_col and std_col in template_data.columns:
            y_err = template_data[std_col]
            ax.errorbar(x_jitter, y_data, yerr=y_err,
                       fmt='none', ecolor=colors[i], alpha=0.5)
    
    ax.set_yscale('log')
    ax.set_xticks(range(len(templates)))
    ax.set_xticklabels(templates, rotation=45, ha='right')
    ax.set_ylabel('Diffusivity (cm²/s)', fontsize=PLOT_PARAMS['label_size'])
    ax.set_title(title, fontsize=PLOT_PARAMS['title_size'])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
        
    return fig


def _passes_screening_criteria(row: pd.Series, criteria: Dict[str, Any]) -> bool:
    """Check if a structure passes screening criteria."""
    # Example criteria checking - customize as needed
    diffusivity_threshold = criteria.get('min_diffusivity', 1e-10)
    max_msd_threshold = criteria.get('max_msd', 20.0)
    
    passes = row.get('diffusivity', 0) >= diffusivity_threshold
    
    # Check MSD criteria if available
    for col in ['msd1', 'msd2', 'msd3']:
        if col in row and row[col] > max_msd_threshold:
            passes = False
            break
    
    return passes