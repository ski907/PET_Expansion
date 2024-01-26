The first script to run is gee_ERA_download.py

Enter the lat/lon and Location info in the "locations" list, this can accomodate several stations

This will download hourly total precipitation data from the ERA5-Land dataset 
from Google Earth Engine for the period 1950-2023 and stores it in a local CSV file
https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY

The script requires a Google Earth Engine account:
https://earthengine.google.com/

And the google earth engine python library
https://developers.google.com/earth-engine/guides/python_install



Once the script generates the 'ERA5_land_hourly_total_precipitation_{station_name}.csv' output file, 
a second script is used to conduct the Generalized Extreme Value fit to the data.

ERA5_GEV_toDSS_GEE.py

The script uses methods extracted from the PXR source code 
primarily from: https://github.com/lrntct/pxr/blob/master/ev_fit.py

The user must specify the "installation", "duration" and "T" (return period) of the desired output hyetograph
These parameters are near the top of the main() function
The script will process the ERA5-Land data with that installation name and create an alternating block hyetograph
for the duration and return period specified, then output that to a DSS file.

This script requires the pydsstools library which can be tricky to install. See the Dependencies section of 
the PXRExtract README for more info: https://github.com/ski907/PXRextract

