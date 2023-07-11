import json
import pandas as pd
import numpy as np
import refet
import math
from functools import cached_property

class ETRCalculator:
    def __init__(self, config_file_path):
        with open(config_file_path, 'r') as f:
            config = json.load(f)

        self.csv_file_path = config['csv_file_path']
        self.elevation = config['elevation']
        self.zw = config['zw']
        self.lon = config['lon']
        self.lat = config['lat']
        self.df = None

    @cached_property
    def df_merged(self):
        self.df = pd.read_csv(self.csv_file_path)

        # Convert the time column to DateTime type
        self.df['time'] = pd.to_datetime(self.df['time'])

        # Set the time column as the index
        self.df.set_index('time', inplace=True)

        # Resample the DataFrame to daily frequency
        df_daily_max = self.df.resample('D').max()
        df_daily_min = self.df.resample('D').min()

        # Repeat the daily values for each hour of the day
        df_daily_max = df_daily_max.reindex(self.df.index, method='ffill')
        df_daily_min = df_daily_min.reindex(self.df.index, method='ffill')
        
        # Rename the columns to reflect the max and min values
        max_columns = [f'{col}_max' for col in self.df.columns]
        min_columns = [f'{col}_min' for col in self.df.columns]
        df_daily_max.columns = max_columns
        df_daily_min.columns = min_columns
        
        # Concatenate the max and min DataFrames
        df_daily = pd.concat([df_daily_max, df_daily_min], axis=1)
        
        # Merge the daily DataFrame back with the original hourly DataFrame
        df_merged = pd.concat([self.df, df_daily], axis=1)

        # Calculate relative humidity for each row
        df_merged["Relative_Humidity"] = df_merged.apply(self.calculate_relative_humidity, axis=1)

        df_merged["Actual_Vapor_Pressure"] = df_merged.apply(self.calculate_actual_vapor_pressure, axis=1)

        df_merged["day_of_year"] = df_merged.index.dayofyear

        df_merged["utc_hour"] = df_merged.index.tz_localize(None).hour

        return df_merged

    def calculate_relative_humidity(self, row):
        specific_humidity = row["SPFH_2maboveground"]
        temperature = row["TMP_2maboveground"] - 273.15
        pressure = row["PRES_surface"]/100
        p_sat = 6.112 * math.exp((17.67 * temperature) / (temperature + 243.5))
        relative_humidity = (specific_humidity / (1 - specific_humidity)) * (p_sat / pressure) * 100
        return relative_humidity
    
                        
    def calculate_actual_vapor_pressure(self, row):
        specific_humidity = row["SPFH_2maboveground"]
        temperature_kelvin = row["TMP_2maboveground"]
        temperature_celsius = temperature_kelvin - 273.15  # Convert from Kelvin to Celsius
        pressure = row["PRES_surface"]
        relative_humidity = row["Relative_Humidity"]

        # Calculate saturation vapor pressure (e_s) using the Clausius-Clapeyron equation
        e_s = 6.112 * math.exp((17.67 * temperature_celsius) / (temperature_celsius + 243.5))

        # Calculate actual vapor pressure (e)
        actual_vapor_pressure = relative_humidity/100 * e_s

        return actual_vapor_pressure

    
    
    def calculate_etr(self):
        rs = self.df_merged["DSWRF_surface"] * 0.086
        temperature_c = self.df_merged["TMP_2maboveground"] - 273.15
        temperature_f = (self.df_merged["TMP_2maboveground"] - 273.15) * 9/5 + 32
        uz = self.df_merged["VGRD_10maboveground"]
        day_of_year = self.df_merged.index.dayofyear
        time = self.df_merged["utc_hour"]
        ea = self.df_merged["Actual_Vapor_Pressure"] / 1000
        wind_speed = np.sqrt(self.df_merged["VGRD_10maboveground"] **2 +  self.df_merged["UGRD_10maboveground"] **2)

        etr = refet.Hourly(
            tmean=temperature_f, ea=ea, rs=rs, uz=wind_speed, zw=self.zw, elev=self.elevation,
            lat=self.lat, lon=self.lon, doy=day_of_year, time=time, method='asce',
            input_units={'tmean': 'F', 'rs': 'Langleys', 'uz': 'mph', 'lat': 'deg'}
        ).etr()
        
        etr = etr * 0.001  #mm/hr to m/hr
        etr = np.clip(etr, 0, None)
        
        return etr
