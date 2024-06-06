# pybvalue
This repository contains spatial mapping of b-value of earthquakes at Salton Trough and Imperial valley, Southern California. 
Workflow for spatail mapping of b-value of earthquakes includes the following steps: 
## 1. Specify the area of interest and gather follwoing information prior to access data: ##

To define the region of interest for example, the Salton Trough/Imperial Valley area, you need to specify the following parameters:
#### Geographic Boundaries: ####
- Minimum Latitude = 32.00 
- Maximum Latitude = 36.00
- Minimum Longitude = -120.0
- Maximum Longitude = -111.0
#### Earthquake Characteristics: ####
- Minimum Magnitude = -1.0
- Maximum Magnitude = 9.0
- Minimum Depth  = -5.0
- Maximum Depth  = 30.0
#### Time Boundaries: ####
- Start Time = 1977:01:01:00
- End Time   = 2024:06:07:00 (Today)

use above information when accessing earthquake catalog data. you can excess catalog data for example from Southern California Earthquake Data Center (SCEDC) from [here](https://service.scedc.caltech.edu/eq-catalogs/date_mag_loc.php) After accessing data you can save under the data directory. 

## 2. Prepare catalog for analysis ##
b-value estimation include the following steps: 
- analysze magntidue distribution the follow Gutenberg Richeter relationship: $\log(N)= a - bM$
- estimate magnitude of completeness of the distribution
- fit GR distribution

for spatial b-value mapping we will implement approach by Schorlemmer et. al., 2004


### Use the following references when using this code: ###

Aki, K., 1965, Maximum likelihood estimate of b in the formula log N = a - bM and its confidence limits: Bull. Earthquake Res. Inst., Tokyo Univ., v. 43, p. 237–239.

Clauset, A., Shalizi, C.R., and Newmann, M.E.J., 2009, Power-law distributions in empirical data: SIAM review, v. 51, no. 4, p. 661–703.

Schorlemmer, D., Wiemer, S., Wyss, M., & Jackson, D. D. (2004). Earthquake statistics at parkfield: 2. probabilistic forecasting and testing. Journal of Geophysical Research: Solid Earth, 109(12), 1–12. (Publisher: Blackwell Publishing Ltd) doi: 10.1029/2004JB003235