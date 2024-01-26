#-------------------------------------------------------------------------------
# Name          Global Precip Tool
# Description:  Tool to obtain precip intensity and temporal distribution
#               globally from GEV fits to ERA5 Land Precip Data from
#               Google Earth Engine 
#               https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY
#               Then output to HEC-DSS format for use in HEC-RAS modeling
#               GEV fit methods follow those from PXR study from https://github.com/lrntct/pxr
# Author:       Chandler Engel
#               US Army Corps of Engineers
#               Cold Regions Research and Engineering Laboratory (CRREL)
#               Chandler.S.Engel@usace.army.mil
# Created:      26 January 2024
# Updated:      26 January 2024
#               
#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
from courty_gev_fit_methods import *

def gev_fit_ERA5data(ERA5_filename='ERA5_land_hourly_total_precipitation.csv'):
    df = pd.read_csv(ERA5_filename)
    # Convert the 'date' column to datetime and set it as the index
    df['date'] = pd.to_datetime(df['date'])
    df.set_index('date', inplace=True)

    ###compare with Courty's POR:
    #comment out to use whole POR in csv
    #start_date = '1979-01-01'
    #end_date = '2018-12-31'
    #df = df.loc[start_date:end_date]
    ###############


    # Dictionary to store the GEV parameters for each interval
    gev_parameters_by_hour = {}

    durations = [1,2,3,4,6,8,10,12,18,24,48,72,96,120,144,192,240,288,360]

    # Loop over the hourly intervals
    for duration in durations:  # 1 to 360 hours
        # Calculate the rolling sum for the current interval
        column_name = f'rolling_{duration}h_precip'
        df[column_name] = df['total_precipitation_hourly'].rolling(f'{duration}H').sum()

        # Group by year and get the maximum for each year
        # also convert to mm
        max_values_by_year = df[column_name].resample('Y').max() * 1000

        max_intensity_by_year = max_values_by_year/duration

        # Fit the GEV distribution to the max values
        gev_params = fit_gev(max_intensity_by_year.values)  # assuming fit_gev takes a numpy array and returns a tuple of parameters

        # Store the parameters in the dictionary
        gev_parameters_by_hour[duration] = gev_params

    # Convert the dictionary to a DataFrame for easier manipulation and viewing
    gev_parameters_df = pd.DataFrame.from_dict(gev_parameters_by_hour, orient='index', columns=['Location', 'Scale', 'Shape'])

    #print(gev_parameters_df)
    
    return gev_parameters_df

def lnt(T):
    return -np.log(1 - 1/T)

def y_gev_nonzero(T, shape):
    """
    Lu, L.-H., & Stedinger, J. R. (1992).
    Variance of two- and three-parameter GEV/PWM quantile estimators:
    formulae, confidence intervals, and a comparison.
    Journal of Hydrology, 138(1–2), 247–267.
    http://doi.org/10.1016/0022-1694(92)90167-T
    """
    return (1 - lnt(T)**shape) / shape

def calculate_gee_duration_intensities(gev_parameters_df,T):
    """Extract all duration intensities for a given recurrance interval"""
    
    y = y_gev_nonzero(T,-0.114)
    
    intensities = gev_parameters_df.Location.values + gev_parameters_df.Scale.values * y
    
    return intensities

def export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,dt,output_file):
    
    #from datetime import datetime
    from pydsstools.heclib.dss import HecDss
    from pydsstools.core import TimeSeriesContainer,UNDEFINED

    precip = np.diff(cum_rainfall_dist) #convert cumulative rainfall to incremental
    precip = np.insert(precip,0,0)
    
    if dt < 60:
        time_string = str(dt)+"MIN"
    elif dt == 60:
        time_string = "1HOUR"

    dss_file = output_file
    pathname = r"/"+installation+r"/"+installation+r"/PRECIP_INC//"+time_string+r"/"+Part_F+r"/"
    tsc = TimeSeriesContainer()
    tsc.pathname = pathname
    tsc.startDateTime = "01JAN2000 00:00:00"
    tsc.numberValues = len(precip)
    tsc.units = "mm"
    tsc.type = "PER-CUM"
    tsc.interval = 1
    tsc.values = precip




    fid = HecDss.Open(dss_file,version=6)
    fid.deletePathname(tsc.pathname)
    fid.put_ts(tsc)
    ts = fid.read_ts(pathname)
    fid.close()

def create_24_hr_alternating_block_hyetograph(intensities):
    
    durations = [1,2,3,4,6,8,10,12,18,24,48,72,96,120,144,192,240,288,360]
    
    precip = np.zeros(24)
    
    dt = np.diff(durations)
    
    cum_depth = intensities*durations
    inc_depth_irr = np.diff(cum_depth)/dt
    
    inc_depth = np.array(cum_depth[0])
    
    for i,d in enumerate(dt):
        inc_depth = np.append(inc_depth,np.ones(d)*inc_depth_irr[i])
  
    alt_block_index = [12,13,11,14,10,15,9,16,8,17,7,18,6,19,5,20,4,21,3,22,2,23,1,0]    
        
    for i,ind in enumerate(alt_block_index):
        precip[ind] = inc_depth[i]
    #test below
    precip = np.append(precip,precip[-1])       
    
    return precip
 
def generate_alt_block_index(duration):
    # Calculate the starting point
    start_index = duration // 2

    # Initialize the alt_block_index with the starting index
    alt_block_index = [start_index]

    # Alternate incrementing and decrementing
    for i in range(1, start_index + 1):
        # Add next higher index if within duration
        if start_index + i < duration:
            alt_block_index.append(start_index + i)

        # Add next lower index if it's non-negative
        if start_index - i >= 0:
            alt_block_index.append(start_index - i)

    return alt_block_index

def create_alternating_block_hyetograph(intensities,total_duration):
    
    durations = [1,2,3,4,6,8,10,12,18,24,48,72,96,120,144,192,240,288,360]
    
    precip = np.zeros(total_duration)
    
    dt = np.diff(durations)
    
    cum_depth = intensities*durations
    inc_depth_irr = np.diff(cum_depth)/dt
    
    inc_depth = np.array(cum_depth[0])
    
    for i,d in enumerate(dt):
        inc_depth = np.append(inc_depth,np.ones(d)*inc_depth_irr[i])
  
    alt_block_index = generate_alt_block_index(total_duration)
        
    for i,ind in enumerate(alt_block_index):
        precip[ind] = inc_depth[i]
    #test below
    precip = np.append(precip,precip[-1])       
    
    return precip

def main():
    
    installation = 'Freefall'
    data_dir = r'C:\Users\RDCRLCSE\Documents\Python Scripts\PET Expansion'
    ERA_data_file = f'ERA5_land_hourly_total_precipitation_{installation}.csv'

    duration = 24 #storm duration, hours for selecting overall average storm intensity
    T = 100 #return period, years

    Tp=duration/2   #time to peak, hours
    total_duration = duration #storm duration, sets length of storm in output series. Usually the same as duration
 
    
    output_file = 'Global_Precip_'+installation+"_metric.dss"
    

    #this CSV file must be pre-generated ideally by using the "gee_ERA_download.py" file
    #this requires a GEE account, and will pull hourly data from 1950-2023 at a specified point
    #the hourly data are then saved in the CSV format for processing with this script
    gev_parameters_df = gev_fit_ERA5data(ERA_data_file)
    
    ####PXR2 alternating block#######
    intensities = calculate_gee_duration_intensities(gev_parameters_df,T)
    #precip = create_24_hr_alternating_block_hyetograph(intensities)
    precip = create_alternating_block_hyetograph(intensities,total_duration)
    
    cum_rainfall_dist = np.cumsum(precip)
    Part_F = f'{T}-yr {total_duration}-hr ERA5 GEV'
        
    export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,60,output_file) #60 minute blocks
    
    print(intensities)
    
    ###############PXR4 chicago method########
    #V=temporal_dist_PXR4(ds_pxr4,site_coord,total_duration,T,Tp)
    #        
    #cum_rainfall_dist = V
    #Part_F = '100-yr PXR4'
    #    
    #export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,6,output_file) #six minute blocks
    #
    ##############PXR2 depth PXR4 distribution###########
    #intensity_PXR2 = extract_intensity_from_PXR2(ds_pxr2,site_coord,duration,T)
    #V_unit=create_unit_rainfall_dist(V)
    #PXR2_PXR4_dist = gen_cum_precip(V_unit,intensity_PXR2.values[0]*24)
    #Part_F = '100-yr PXR2 Depth PXR4 Dist'
    
    #export_rainfall_toDSS(installation,Part_F,PXR2_PXR4_dist,6,output_file)

    
if __name__ == "__main__":
    main()