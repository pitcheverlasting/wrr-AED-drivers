# wrr-AED-drivers

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li><a href="#about the project">About the Project</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

Atmospheric evaporative demand (AED) is an important variable linking climate with the terrestrial water cycle and the biosphere. Understanding the dynamics of AED has substantial economic, ecological, and social implications. However, how AED varies at different time scales and the drivers of variability remain elusive. This study used spectral coherence analysis to analyze the relationships between observed and modeled AED and climate drivers across multiple time scales at 228 Chinese stations and explored the cross-scale effects of climate forcings on AED. 

### Highlights 

- Vapor pressure deficit (VPD) is essential to predicting atmospheric evaporative demand in both energy-limited and water-limited regions, therefore models that do not incorporate VPD or underestimate the relative importance of VPD have relatively lower skill in predicting AED. Model predictions for AED and associated hydrologic impacts may not be valid in a changing climate when the key controls on AED and their relative importance are not appropriately represented.
- Short-term forcing variability has potential impacts on the long-term AED changes through tempera- ture and associated land-atmosphere feedbacks.

This code repo is a collection of the final scripts for our study.

### Reference
Peng, L., Li, D., &amp; Sheffield, J. (2018). Drivers of Variability in Atmospheric Evaporative Demand: Multiscale Spectral Analysis Based on Observations and Physically Based Modeling. Water Resources Research, 54(5), 3510–3529. http://doi.org/10.1029/2017WR022104


<!-- USAGE EXAMPLES -->
## Usage
* Main analysis: 
  * Ep_obs_analysis.py
  * Ep_mod_analysis.py
  * Experiment.py
* Data processing: Examine data quality, detect step changes
  * Quality_control.py
  * Data_quality.py
* Library:
  * IO.py
  * FFT.py
  * PET_Library.py
  * Plotting.py


<!-- LICENSE -->
## License
Distributed under the Princeton Terrestrial Hydrology Group License. 


<!-- CONTACT -->
## Contact
Project Link: [https://github.com/pitcheverlasting/wrr-AED-drivers](https://github.com/pitcheverlasting/wrr-AED-drivers)

