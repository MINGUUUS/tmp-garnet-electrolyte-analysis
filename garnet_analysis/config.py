"""
Configuration management for garnet electrolyte analysis.

This module provides configuration settings that can be customized
for different analysis workflows.
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional
import json

# Default configuration
DEFAULT_CONFIG = {
    'analysis': {
        'diffusivity': {
            'time_step': 1.0,  # fs
            'step_skip': 500,
            'li_params': {
                'lower_bound': 4.5,
                'upper_bound': 0.99,
                'minimum_msd_diff': 4.5,
                'total_sim_time_limit': 100
            },
            'other_element_params': {
                'lower_bound': 0.1,
                'upper_bound': 0.99,
                'minimum_msd_diff': 0.01,
                'total_sim_time_limit': 0.001
            }
        },
        'screening': {
            'melting_msd_threshold': 20.0,
            'melting_diffusivity_threshold': 6.5e-06,
            'min_diffusivity_threshold': 1e-7,
            'llzo_fraction_threshold': 0.9
        },
        'arrhenius': {
            'default_temperatures': [600, 700, 800, 900, 1000],
            'extrapolation_temperature': 298,
            'min_data_points': 2
        }
    },
    'plotting': {
        'figure_size': [8, 6],
        'dpi': 300,
        'font_size': 12,
        'save_format': 'png',
        'color_scheme': 'default'
    },
    'io': {
        'output_precision': 6,
        'csv_separator': ',',
        'default_encoding': 'utf-8'
    }
}


class Config:
    """Configuration manager for garnet analysis."""
    
    def __init__(self, config_dict: Optional[Dict[str, Any]] = None):
        """
        Initialize configuration.
        
        Args:
            config_dict: Custom configuration dictionary. If None, uses defaults.
        """
        self._config = DEFAULT_CONFIG.copy()
        if config_dict:
            self._config.update(config_dict)
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value using dot notation.
        
        Args:
            key: Configuration key (e.g., 'analysis.diffusivity.time_step')
            default: Default value if key not found
            
        Returns:
            Configuration value
        """
        keys = key.split('.')
        value = self._config
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def set(self, key: str, value: Any) -> None:
        """
        Set configuration value using dot notation.
        
        Args:
            key: Configuration key (e.g., 'analysis.diffusivity.time_step')
            value: Value to set
        """
        keys = key.split('.')
        config = self._config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
    
    def update(self, config_dict: Dict[str, Any]) -> None:
        """
        Update configuration with new values.
        
        Args:
            config_dict: Dictionary of configuration updates
        """
        self._deep_update(self._config, config_dict)
    
    def save(self, file_path: str) -> None:
        """
        Save configuration to JSON file.
        
        Args:
            file_path: Path to save configuration file
        """
        with open(file_path, 'w') as f:
            json.dump(self._config, f, indent=2)
    
    @classmethod
    def load(cls, file_path: str) -> 'Config':
        """
        Load configuration from JSON file.
        
        Args:
            file_path: Path to configuration file
            
        Returns:
            Config object loaded from file
        """
        with open(file_path, 'r') as f:
            config_dict = json.load(f)
        
        return cls(config_dict)
    
    def _deep_update(self, base_dict: Dict, update_dict: Dict) -> None:
        """Recursively update nested dictionary."""
        for key, value in update_dict.items():
            if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
                self._deep_update(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def __getitem__(self, key: str) -> Any:
        """Allow dictionary-style access."""
        return self.get(key)
    
    def __setitem__(self, key: str, value: Any) -> None:
        """Allow dictionary-style setting."""
        self.set(key, value)


# Global configuration instance
_global_config = Config()


def get_config() -> Config:
    """Get the global configuration instance."""
    return _global_config


def set_config(config: Config) -> None:
    """Set the global configuration instance."""
    global _global_config
    _global_config = config


def load_config(file_path: str) -> None:
    """Load configuration from file and set as global config."""
    config = Config.load(file_path)
    set_config(config)


def save_config(file_path: str) -> None:
    """Save current global configuration to file."""
    _global_config.save(file_path)


# Convenience functions for common configuration access
def get_analysis_params(analysis_type: str = 'diffusivity') -> Dict[str, Any]:
    """Get analysis parameters for specified type."""
    return _global_config.get(f'analysis.{analysis_type}', {})


def get_plotting_params() -> Dict[str, Any]:
    """Get plotting parameters."""
    return _global_config.get('plotting', {})


def get_screening_criteria() -> Dict[str, Any]:
    """Get screening criteria."""
    return _global_config.get('analysis.screening', {})