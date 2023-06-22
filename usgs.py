#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    This implements the nwis data retrieval. 
    The data retrieval uses the already exiting python library from https://github.com/USGS-python/dataretrieval. 
    The original code is modified to collect hourly measurements and interact with the data retriever BMI. A validity check is added to inform BMI if there is a measurement. 
    This code is useful for extracting observation from Nextgen framework.
    
    contact
    ----------
    fwolkeba@crimson.ua.edu
    butlerz@oregonstate.edu
    
    Required inputs
    ----------
    site: USGS site ID e.g.  01123000
    start start of the requested data e.g. 2017-05-15
    parameterCd:"00060", parameterCd is constant to collect only flow measurements 
    end: End of the requested e.g. 2017-05-15

    outputs
    ----------
    Flow: Hourly flow measurement of the requested station
    validity: Binary value indicating availability of measurement (0= no measurement, 1 there is a valid measurement)
    Attributes

    References
    ----------
    https://github.com/USGS-python/dataretrieval

    """


import dataretrieval.nwis as nwis
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None


class USGS:

    def __init__(self):
        super(USGS, self).__init__()

    def run_usgs(self, u):

        # -------------------------------------------------------------------

        sites = u.sites
        start = u.start
        end = u.end
        site = nwis.get_record(sites=u.sites, parameterCd='00060',
                                start=u.start, end=u.end)

        # -------------------------------------------------------------------
        # change the extacted values to data frame to do computation on
        # -------------------------------------------------------------------

        #site = pd.DataFrame(site_data[0])
        site.reset_index(inplace=True)                                        # reset index to grab station date
        site['datetime'] = pd.to_datetime(site['datetime'], utc=True,
                format='%Y-%m-%d %H:%M:%S')                                        # transfer to utc so same time throughout
        site_copy = site.copy()
        site_copy['datetime'] = pd.to_datetime(site_copy['datetime'],
                utc=True, format='%Y-%m-%d %H:%M:%S')                                        # convert datetime to average every hour
        site_copy.index = site_copy['datetime']                                         # index so can pull date time in resample
        site_avg = site_copy.select_dtypes(include=['float64', 'int64']).resample('H').mean()
        #site_avg = site_copy.resample('H').mean()                                        # Average every hour based on datetime
        site_avg.reset_index(inplace=True)                                        # reset index again to have datetime
        site_avg.columns = ['Date', 'Flow']

        # check validity of extracted data

        site_avg.loc[site_avg['Flow'] >= 0, 'validity'] = 1                                        # if value positive, consider
        site_avg.loc[site_avg['Flow'] < 0, 'validity'] = 0                                        # if less than zero, not realistic
        site_avg.loc[site_avg['Flow'].isnull() == True, 'validity'] = 0                                        # if NaN not availible

        # Output results to csv file

        #site_avg.to_csv('USGS_' + str(sites) + '_obs_streamflow.csv',
        #                index=False)

        # to check if the code runs on the framework

        u.flow = site_avg['Flow']
        u.validity = site_avg['validity']

        return
