from eto import ETo
import json
import pandas as pd
import numpy as np
from functools import cached_property
import os
import math

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
    def new_df(self):
        self.df = pd.read_csv(self.csv_file_path)

        # Convert the time column to DateTime type
        self.df['time'] = pd.to_datetime(self.df['time'])

        # Set the time column as the index
        self.df.set_index('time', inplace=True)
        

        # Calculate relative humidity for each row
        self.df["Relative_Humidity"] = self.df.apply(self.calculate_relative_humidity, axis=1)

        self.df["Actual_Vapor_Pressure"] = self.df.apply(self.calculate_actual_vapor_pressure, axis=1)

        self.df["day_of_year"] = self.df.index.dayofyear

        self.df["utc_hour"] = self.df.index.tz_localize(None).hour
        
        data = {
            'R_s': self.df["DSWRF_surface"] * 0.0036,
            'T_mean': self.df["TMP_2maboveground"] - 273.15,
            'e_a': self.df["Actual_Vapor_Pressure"] / 1000,
            'RH_mean': self.df["Relative_Humidity"],
            'date': self.df.index
        }

        new_df = pd.DataFrame(data)
        new_df.index = pd.to_datetime(new_df.index)
        
        return new_df
    

    
    def calculate_relative_humidity(self, row):
        specific_humidity = row["SPFH_2maboveground"]
        temperature = row["TMP_2maboveground"] - 273.15
        pressure = row["PRES_surface"]/100
        p_sat = 6.112 * math.exp((17.67 * temperature) / (temperature + 243.5))
        relative_humidity = (specific_humidity / (1 - specific_humidity)) * (p_sat / pressure) * 100
        return relative_humidity    
    
            
    
    def calculate_actual_vapor_pressure(self, row):
        temperature = row["TMP_2maboveground"] - 273.15  # Temperature in Celsius
        specific_humidity = row["SPFH_2maboveground"]  # Specific humidity
        pressure = row["PRES_surface"]  # Atmospheric pressure in pascals

        actual_vapor_pressure = specific_humidity * pressure / (0.622 + 0.378 * specific_humidity)
        return actual_vapor_pressure    
    
  
    
    def calculate_etr(self):
       
        etr = ETo(self.new_df, 'H', z_msl=self.elevation, lat=self.lat, lon=self.lon, TZ_lon=self.lon).eto_fao() 
        etr = etr.fillna(0)
        etr = etr * 0.001  #mm/hr to m/hr  
                  
        return etr

    
   