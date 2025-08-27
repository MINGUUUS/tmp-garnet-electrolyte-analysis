# Garnet Electrolyte Analysis Suite

A comprehensive Python toolkit for computational analysis of garnet-structured solid electrolytes, developed for high-throughput screening and property prediction of Li-ion conducting materials.

## Overview

This repository contains modular, well-documented code for:
- **Diffusivity Analysis**: Calculate Li-ion diffusivity from molecular dynamics trajectories
- **Arrhenius Analysis**: Perform temperature-dependent analysis and extrapolation
- **Energy Analysis**: Formation energy and energy hull calculations
- **Screening Workflows**: High-throughput material screening pipelines
- **Visualization**: Professional plotting functions for publication-quality figures

## Features

- **Modular Design**: Clean separation of functionality into focused modules
- **Publication Ready**: Code used in peer-reviewed journal articles
- **Extensive Documentation**: Comprehensive API documentation and examples
- **Error Handling**: Robust input validation and error reporting
- **Reproducible**: Configuration management for consistent results

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/garnet-electrolyte-analysis.git
cd garnet-electrolyte-analysis

# Create conda environment
conda env create -f environment.yml
conda activate garnet-analysis

# Install package
pip install -e .
```

## Quick Start

```python
import garnet_analysis as ga

# Load MD trajectory and calculate diffusivity
analyzer = ga.DiffusivityAnalyzer.from_lammps_dump('trajectory.dump', 'structure.data')
diffusivity = analyzer.calculate_diffusivity(species='Li', temperature=1000)

# Create Arrhenius analysis
arrhenius = ga.ArrheniusAnalyzer()
arrhenius.add_data_point(1000, diffusivity.value, diffusivity.std)
arrhenius.fit()
D_298K = arrhenius.extrapolate_to_temperature(298)

# Plot results
ga.plot_arrhenius(arrhenius, save_path='arrhenius_plot.png')
```

## Documentation

- [API Documentation](docs/api/)
- [Tutorial Notebooks](examples/)
- [Configuration Guide](docs/configuration.md)

## Citation

If you use this code in your research, please cite:

```
[Your journal article citation here]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Acknowledgments

- Materials Project for structure databases
- AIMD community for trajectory analysis methods
- Reviewers and collaborators for valuable feedback