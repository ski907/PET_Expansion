import ee
import pandas as pd

# Trigger the authentication flow.
ee.Authenticate()

# Initialize the library.
ee.Initialize()

# Define your list of locations and station names
locations = [#(44.564624, -72.761110, 'Freefall'),
             (40.588934, -111.626780, "Sunnyside")
             # Add more tuples of (latitude, longitude, station_name)
            ]

years = range(1950, 2024)
variables = ['total_precipitation_hourly']

# Loop over each location
for u_lat, u_lon, station_name in locations:
    print(f'Processing data for station: {station_name}')

    u_poi = ee.Geometry.Point(u_lon, u_lat)
    
    # Container for yearly dataframes
    d_temp = []

    # Loop for each year
    for year in years:
        print(f'Getting year {year}...')
        i_date = f'{year}-01-01'
        f_date = f'{year+1}-01-01'
        era = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
        era = era.select(variables).filterDate(i_date, f_date)
        scale = 1
        era_u_poi = era.getRegion(u_poi, scale).getInfo()

        df = pd.DataFrame(era_u_poi[1:], columns=era_u_poi[0])
        df['date'] = pd.to_datetime(df.id)
        df = df.set_index('date')
        
        df['year'] = df.index.year
        df['doy'] = df.index.dayofyear
        df['water_year'] = df.index.year.where(df.index.month < 10, df.index.year + 1)
        
        d_temp.append(df)
        print(f'{year} complete!')

    # Combine all the yearly dataframes
    df_era5_hr = pd.concat(d_temp)

    # Export the dataframe
    output_filename = f'ERA5_land_hourly_total_precipitation_{station_name}.csv'
    df_era5_hr.to_csv(output_filename)
    print(f'Data for {station_name} saved to {output_filename}')