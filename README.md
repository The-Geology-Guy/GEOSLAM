# SEISMO-LINEAMENT ANALYSIS METHOD (SLAM)

Main Project Resources: [PAJER, Luke](mailto:luke.pajer@gmail.com); [CRONIN, Vincent](mailto:vince_cronin@baylor.edu)

_Last Updated: October 2020_

[![croninprojects.org home](https://img.shields.io/badge/croninprojects.org-home-F78C26.svg)](http://croninprojects.org/)
[![License](https://img.shields.io/badge/LICENSE-mit-43B02A.svg)](/LICENSE)
[![jupyterlab](https://img.shields.io/badge/jupyterlab-0.35.4-F37821.svg)](https://jupyterlab.readthedocs.io/en/stable/)
[![python](https://img.shields.io/badge/python-3.6.5-yellow.svg)](https://jupyterlab.readthedocs.io/en/stable/)

-----

# PROJECT OVERVIEW

The Seismo-Lineament Analysis Method (SLAM) Python package is a simple translation of the mathematica code developed by Vince Cronin. The purpose of this package is to make SLAM available to those who are interested in SLAM, but are more comfortable using Python.

_From [Vince Cronin's website](https://croninprojects.org/Vince/SLAM/index.htm):_

> The seismo-lineament analysis method is a tool to spatially correlate a shallow-focus earthquake to the surface trace of the fault that generated it. SLAM is the intellectual property and work product of Vince Cronin, and has been developed with assistance from Brandon Rasaka, Victoria Worrell, Jeremy Ashburn, Brian Bayliss, Chris Breed, Bruce Byars, Ryan Campbell, David Cleveland, Jon Cook, Kelly Cronin, Jordan Dickinson, Daniel Lancaster, Ryan Lindsay, Mark Millard, Shane Prochnow, Tyler Reed, Stephen Secrest, Lauren Seidman Robinson, Keith Sverdrup and Lisa Zygo, with funds from AAPG, Baylor University, Colorado Scientific Society, Ellis Exploration, Ft. Worth Geological Society, Geological Society of America, GCAGS, Samson Resources, Roy Shlemon Scholarship Fund, Sigma Xi, and SIPES.

If there are any issues or concerns with the python package, please reach out to [Luke Pajer](mailto:luke.pajer@gmail.com). For any questions regarding SLAM, please reach out to [Vince Cronin](mailto:vince_cronin@baylor.edu).

## CONTRIBUTORS

This project is an open project, and contributions are welcome from any individual. All contributors to this project are bound by a [code of conduct](/CODE_OF_CONDUCT.md). Please review and follow this code of conduct as part of your contribution.

#### Contributions to the SLAM Python Package
- [Luke Pajer](mailto:luke.pajer@gmail.com) [![orcid](https://img.shields.io/badge/orcid-0000--0002--5218--7650-brightgreen.svg)](https://orcid.org/0000-0002-5218-7650)

#### SLAM Author/Developer
- [Vince Cronin](mailto:vince_cronin@baylor.edu) [![orcid](https://img.shields.io/badge/orcid-0000--0002--3069--6470-brightgreen.svg)](https://orcid.org/0000-0002-3069-6470)

In addition, there are a few thesis projects used when developing the SLAM package. These were instrumental in developing and troubleshooting the package. Below are the Thesis projects referenced during development:

- Victoria E. Worrell: ["The Seismo-Lineament Analysis Method (SLAM) Applied to the South Napa Earthquake and Antecedent Events"](https://baylor-ir.tdl.org/bitstream/handle/2104/9796/WORRELL-THESIS-2016.pdf?sequence=1&isAllowed=y) 
- Jeremy A. Ashburn: ["Investigation of a Lineament that Might Mark the Ground-Surface Trace of the Dog Valley Fault, Truckee Area, Northern California"](https://croninprojects.org/Vince/AshburnBSThesis2015.pdf) 
- Brandon M. Rasaka: ["Correlation of Selected Earthquakes with Seismogenic Faults, Central Oklahoma"](https://croninprojects.org/Rasaka/Rasaka-MS-Thesis-2016.pdf)

### Tips for Contributing

Issues and bug reports are always welcome.  Code clean-up, and feature additions can be done either through pull requests to [project forks]() or branches.

All products of the SLAM project are licensed under an [MIT License](LICENSE) unless otherwise noted.

-----

## HOW TO USE THIS REPOSITORY

This repository is available to be 

Base overview for SLAM map generation (_see the [SLAM Wiki](https://github.com/The-Geology-Guy/SLAM/wiki) for more information_):
1. Get a DEM file for a specific area (_max lat/lon and min lat/lon must be set via variables_)  
2. Instantiate qFaults shapefile for quaternary faults and folds to be plotted on map  
3. Gather Event data (_you can either use the SLAM function to query or enter data yourself_)  
4. Calculate Grid North Adjustment  
5. Produce computations relative to the 1st Nodal Plane  
6. Determine errors relative to the 1st Nodal Plane  
7. Calculate light direction for hill-shade map  
8. Calculate the 7 swaths using the `swath_calc` function  
9. Determine Seismo-Lineament Boundaries using the `get_shade` function (_there also may be a need to fill in some gaps in the boundary area, to do this, specify where the gaps are on the map via the 'corners_fill' attribute -- see the [SLAM Wiki](https://github.com/The-Geology-Guy/SLAM/wiki) for more information_)  
10. Plot the focal mechanism and Seismo-Lineament boundaries on a map using any of the map types available -- (_see the [SLAM Wiki](https://github.com/The-Geology-Guy/SLAM/wiki) for more information_)

Once again, this is a simple overview of a typical SLAM task. This is in no way the limit of what can be done. See the [SLAM Wiki](https://github.com/The-Geology-Guy/SLAM/wiki) for more information.

### System Requirements

This project is developed using Python. There should be no issues with these projects running on Mac, Windows, or Linux. If there are any issues, please submit an issue and it will be investigated.

### Data Resources used in SLAM

#### Data Sources

#### Physical Maps and Digital Elevation Model (DEM) Sources

- [OpenTopography](https://opentopography.org/) This work is based on [data, processing] services provided by the OpenTopography Facility with support from the National Science Foundation under NSF Award Numbers 1948997, 1948994 & 1948857.

OpenTopography also hosts and makes available some global data, _all of which are accessible via the OpenTopography services and made available with functions found in the GEOSLAM code_. Below are the citations for the Global Data: 

_GMRT:_ ([Terms of Use](https://www.marine-geo.org/about/terms_of_use.php))
> Ryan, W.B.F., S.M. Carbotte, J.O. Coplan, S. O'Hara, A. Melkonian, R. Arko, R.A. Weissel, V. Ferrini, A. Goodwillie, F. Nitsche, J. Bonczkowski, and R. Zemsky (2009), Global Multi-Resolution Topography synthesis, Geochem. Geophys. Geosyst., 10, Q03014. https://doi.org/10.1029/2008GC002332  

_ALOS World 3D:_
> J. Takaku, T. Tadono, K. Tsutsui : Generation of High Resolution Global DSM from ALOS PRISM, The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, pp.243-248, Vol. XL-4, ISPRS TC IV Symposium, Suzhou, China, 2014. 
> Here is a link to a [PDF file](https://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XL-4/243/2014/isprsarchives-XL-4-243-2014.pdf) of this publication.

_SRTM:_
> Farr, T.G., Rosen, P.A., Caro, E., Crippen, R., Duren, R., Hensley, S., Kobrick, M., Paller, M., Rodriguez, E., Roth, L., Seal, D., Shaffer, S., Shimada, J., Umland, J., Werner, M., Oskin, M., Burbank, D., Alsdorf, D. (2007), The Shuttle Radar Topography Mission, Rev. Geophys., 45, RG2004. https://doi.org/10.1029/2005RG000183

- [Stamen Map Tile Sets](http://maps.stamen.com/#watercolor/12/37.7706/-122.3782) are used to generate the physical maps in this package. The Stamen map tile sets are copyright Stamen Design, under a Creative Commons Attribution (CC BY 3.0) license.

### Key Outputs

This project generates (an API, some log files, what?)

## REFERENCES AND NOTABLE PACKAGES USED


- https://croninprojects.org/Rasaka/Rasaka-MS-Thesis-2016.pdf  
- https://baylor-ir.tdl.org/bitstream/handle/2104/9796/WORRELL-THESIS-2016.pdf?sequence=1&isAllowed=y  
- https://croninprojects.org/Vince/SLAM/CurrentBaseCodes.html  
- http://serc.carleton.edu/files/NAGTWorkshops/structure04/Focal_mechanism_primer.pdf  
- https://croninprojects.org/Vince/SLAM/SLAMWorkflow.html  
- https://portal.opentopography.org/apidocs/#/Public/getGlobalDem  
- https://rasterio.readthedocs.io/en/latest/  
- https://earthquake.usgs.gov/fdsnws/event/1/  
- https://docs.obspy.org/  

-----

