import xarray as xr
import numpy as np
from xclim.indices import tx10p, tx90p, tn10p, tn90p, warm_spell_duration_index, cold_spell_duration_index, frost_days
from xclim.core.calendar import percentile_doy
import warnings

# Helper function
def first_run(data, window, dim='time'):
    """
    Identify and count consecutive sequences (runs) of days that satisfy a condition.

    :param data: The xarray data containing the variable to check.
    :type data: xarray.DataArray
    :param window: The minimum length of a run to be counted.
    :type window: int
    :param dim: The dimension along which to find runs (default is 'time').
    :type dim: str, optional

    :return: An xarray DataArray containing the length of runs that satisfy the condition.
    :rtype: xarray.DataArray

    .. note::
       The function assumes that the condition to check is already applied to the data passed.
       
    :Example:
    >>> data = xr.DataArray([True, True, False, True, True, True], dims=['time'])
    >>> first_run(data, window=2)
    """
    runs = []
    counter = 0
    for day in data[dim].values:
        if data.sel({dim: day}).any():
            counter += 1
        else:
            if counter >= window:
                runs.append(counter)
            counter = 0
    if counter >= window:
        runs.append(counter)
    return xr.DataArray(runs, dims=[dim])

# Helper function to calculate the length of the longest run that satisfies a condition
def longest_run(arr, condition_func):
    """
    Calculate the length of the longest run in an array that satisfies a condition.

    :param arr: The array containing the elements to check.
    :type arr: iterable
    :param condition_func: The function that defines the condition to check.
    :type condition_func: callable

    :return: The length of the longest run that satisfies the condition.
    :rtype: int
    
    :Example:
    >>> longest_run([True, True, False, True, True, True], lambda x: x)
    """
    longest = 0
    current_run = 0
    for val in arr:
        if condition_func(val):
            current_run += 1
            longest = max(longest, current_run)
        else:
            current_run = 0
    return longest

class ClimateIndex:
    """
    A class to handle climate index calculations.
    
    :param datafile: The path to the climate data file in netCDF format.
    :type datafile: str
    
    .. note::
       The class is designed to load climate data and provide various index calculations.
       After initializing, use the `pre_process` method to preprocess the data before calculating any indices.
       
    :Example:
    >>> climate_index = ClimateIndex("path/to/datafile.nc")
    
    :ivar datafile: The path to the climate data file.
    :vartype datafile: str
    :ivar data: The loaded climate data, initialized as None and populated using the `pre_process` method.
    :vartype data: xarray.Dataset or None
    """
    def __init__(self, datafile):
        """
        Initialize the ClimateIndex class with the data file path.
        
        :param datafile: The path to the climate data file in netCDF format.
        :type datafile: str
        
        :return: None
        :rtype: None
        """
        self.datafile = datafile
        self.data = None

    def pre_process(self, time_range=None, fill_missing=None, resample_freq=None, months=None, season=None, resample_method='mean'):
        """
        Pre-process the climate data for further analysis.
        
        :param time_range: A tuple specifying the start and end time for data selection. Example: ('1961-01-01', '1990-12-31'). Default is None.
        :type time_range: tuple, optional
        :param fill_missing: The method for filling missing values, either a numerical value or 'interpolate'. Example: 0. Default is None.
        :type fill_missing: int, float, or str, optional
        :param resample_freq: The frequency for resampling the data. Example: 'MS'. Default is None.
        :type resample_freq: str, optional
        :param months: A list of integers specifying the months for data selection. Example: [1, 2, 3]. Default is None.
        :type months: list of int, optional
        :param season: A string specifying the season for data selection ('DJF', 'MAM', 'JJA', 'SON'). Example: 'JJA'. Default is None.
        :type season: str, optional
        :param resample_method: The method for resampling, should be a method of xarray's DataArray resampler. Default is 'mean'.
        :type resample_method: str, optional
        
        :return: None
        :rtype: None
        
        .. note::
        The function updates the instance's data attribute with the pre-processed data.
        - **Units Check**: The function supports 'K' for temperature and 'm' for precipitation.
        - **Special Rule**: If season is 'DJF', and `resample_freq` is not specified, it defaults to 'QS-DEC'.
        
        :Example:
        >>> climate_index.pre_process(time_range=('2000-01-01', '2020-12-31'), resample_freq='M', fill_missing='interpolate')
        
        :Warnings:
            An exception will be caught and a warning will be issued if the pre-processing cannot be performed.
        """
        try:
            print("Loading data...")
            self.data = xr.open_dataset(self.datafile)
            print("Data loaded successfully.")
            
            if 'time' not in self.data.dims:
                raise ValueError("Time dimension not found in the dataset.")
            
            print("Time axis found.")
            
            # Initialize data_type as None
            self.data_type = None
            
            # Infer data type based on units
            for var in self.data.data_vars:
                units = self.data[var].attrs.get('units', None)
                if units:
                    print(f"Data detected in variable {var}. Units are in {units}.")
                    if units == 'K':
                        self.data_type = 'temperature'
                    elif units == 'm':
                        self.data_type = 'precipitation'
                    else:
                        raise ValueError(f"Unsupported units {units}. Expected units are 'K' for temperature or 'm' for precipitation.")
            
            if self.data_type is None:
                raise ValueError("Could not infer data type. Please check the units in the dataset.")
            
            if time_range:
                start_time, end_time = time_range
                self.data = self.data.sel(time=slice(start_time, end_time))
                print("Time selection completed.")
                
            if months:
                self.data = self.data.sel(time=self.data['time.month'].isin(months))
                print("Month selection completed.")
                
            if season:
                self.data = self.data.sel(time=self.data['time.season'] == season)
                print("Season selection completed.")
                
                # Special resampling for DJF season
                if season == 'DJF' and resample_freq is None:
                    resample_freq = 'QS-DEC'
                
            if resample_freq:
                resampler = self.data.resample(time=resample_freq)
                if resample_method in dir(resampler):
                    self.data = getattr(resampler, resample_method)()
                else:
                    raise ValueError(f"Unknown resampling method: {resample_method}.")
                print("Resampling completed.")
                
            if fill_missing:
                if isinstance(fill_missing, (int, float)):
                    print(f"Filling missing values with specified value: {fill_missing}...")
                    self.data = self.data.fillna(fill_missing)
                elif fill_missing == 'interpolate':
                    print("Filling missing values using interpolation...")
                    self.data = self.data.interpolate_na(dim='time', method='linear')
                else:
                    print(f"Unknown method to fill missing values: {fill_missing}. Skipping this step.")
                print("Missing values handling completed.")
                
        except Exception as e:
            print(f"An error occurred while processing the data: {str(e)}")


    def calculate_tx10p(self, temp_var='t2m'):
        """
        Calculate the number of days when maximum temperature exceeds the 10th percentile (TX10P).

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. Default is 't2m'.
        :type temp_var: str, optional

        :return: An xarray Dataset containing the TX10P index.
        :rtype: xarray.Dataset

        .. note::
        The 10th percentile is calculated over the period from 1961 to 1990.

        :Example:
        >>> tx10p_index = climate_index.calculate_tx10p()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        tasmax_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=10, window=5).sel(percentiles=10)
        result = self.data[temp_var].where(self.data[temp_var] < tasmax_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tx10p"
        return xr.Dataset({"tx10p": percentage_days})

    def calculate_tx90p(self, temp_var='t2m'):
        """
        Calculate the number of days when maximum temperature exceeds the 90th percentile (TX90P).

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. Default is 't2m'.
        :type temp_var: str, optional

        :return: An xarray Dataset containing the TX90P index.
        :rtype: xarray.Dataset

        .. note::
        The 90th percentile is calculated over the period from 1961 to 1990.

        :Example:
        >>> tx90p_index = climate_index.calculate_tx90p()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        tasmax_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=90, window=5).sel(percentiles=90)
        result = self.data[temp_var].where(self.data[temp_var] > tasmax_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tx90p"
        return xr.Dataset({"tx90p": percentage_days})

    def calculate_tn10p(self, temp_var='t2m'):
        """
        Calculate the number of days when minimum temperature is below the 10th percentile (TN10P).

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. Default is 't2m'.
        :type temp_var: str, optional

        :return: An xarray Dataset containing the TN10P index.
        :rtype: xarray.Dataset

        .. note::
        The 10th percentile is calculated over the period from 1961 to 1990.

        :Example:
        >>> tn10p_index = climate_index.calculate_tn10p()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        tasmin_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=10, window=5).sel(percentiles=10)
        result = self.data[temp_var].where(self.data[temp_var] < tasmin_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tn10p"
        return xr.Dataset({"tn10p": percentage_days})

    def calculate_tn90p(self, temp_var='t2m'):
        """
        Calculate the number of days when minimum temperature exceeds the 90th percentile (TN90P).

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. Default is 't2m'.
        :type temp_var: str, optional

        :return: An xarray Dataset containing the TN90P index.
        :rtype: xarray.Dataset

        .. note::
        The 90th percentile is calculated over the period from 1961 to 1990.

        :Example:
        >>> tn90p_index = climate_index.calculate_tn90p()

        :Warnings:
        An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        tasmin_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=90, window=5).sel(percentiles=90)
        result = self.data[temp_var].where(self.data[temp_var] > tasmin_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tn90p"
        return xr.Dataset({"tn90p": percentage_days})


    def calculate_frost_days(self, temp_var='t2m'):
        """
        Count the number of days when the minimum temperature falls below 0°C (Frost Days).

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. Default is 't2m'.
        :type temp_var: str, optional

        :return: An xarray Dataset containing the count of Frost Days.
        :rtype: xarray.Dataset

        :Example:
        >>> frost_days_index = climate_index.calculate_frost_days()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        frost_days_dat = frost_days(self.data[temp_var], thresh='273.15 K')
        frost_days_dat.name = 'frost_days_dat'
        return xr.Dataset({'frost_days_dat': frost_days_dat})

    def calculate_warm_spell(self, temp_var='t2m'):
        """
        Calculate the Warm Spell Duration Index (WSDI).

        This method measures the duration of warm spells when the daily maximum temperature exceeds the 90th percentile for consecutive days.

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter. The 90th percentile is calculated over the period from 1961 to 1990.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the Warm Spell Duration Index.
        :rtype: xarray.Dataset

        :Example:

        >>> warm_spell_index = climate_index.calculate_warm_spell()

        :Warnings:
        An exception will be caught and a warning will be issued if the calculation cannot be performed.

        """
        tasmax_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=90, window=5).sel(percentiles=90)
        warm_days = self.data[temp_var] > tasmax_per
        warm_spell_days = first_run(warm_days, window=6, dim='time')
        warm_spell_days.name = 'wsdi'
        return xr.Dataset({'wsdi': warm_spell_days})

    def calculate_cold_spell(self, temp_var='t2m'):
        """
        Calculate the Cold Spell Duration Index (CSDI).

        This method measures the duration of cold spells when the daily minimum temperature falls below the 10th percentile for consecutive days.

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter. The 10th percentile is calculated over the period from 1961 to 1990.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the Cold Spell Duration Index.
        :rtype: xarray.Dataset

        :Example:

        >>> cold_spell_index = climate_index.calculate_cold_spell()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.

        """
        tasmin_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=10, window=5).sel(percentiles=10)
        cold_days = self.data[temp_var] < tasmin_per
        cold_spell_days = first_run(cold_days, window=6, dim='time')
        cold_spell_days.name = 'csdi'
        return xr.Dataset({'csdi': cold_spell_days})

    
    def calculate_dtr(self, temp_var='t2m'):
        """
        Calculate the Mean Monthly Diurnal Temperature Range (DTR).

        This method calculates the mean monthly diurnal temperature range, which is the difference between the daily maximum and minimum temperature averaged over each month.

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the mean monthly diurnal temperature range.
        :rtype: xarray.Dataset

        :Example:

        >>> dtr_index = climate_index.calculate_dtr()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.

        """
        daily_max = self.data[temp_var].resample(time='D').max(dim='time')
        daily_min = self.data[temp_var].resample(time='D').min(dim='time')
        dtr = daily_max - daily_min
        dtr_monthly_mean = dtr.resample(time='M').mean(dim='time')
        dtr_monthly_mean.name = 'dtr'
        return xr.Dataset({'dtr': dtr_monthly_mean})

    def calculate_su(self, temp_var='t2m'):
        """
        Calculate the Number of Summer Days (SU) in a year.

        This method calculates the number of days in a year when the daily maximum temperature exceeds 25°C (298.15 K).

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the number of summer days for each year.
        :rtype: xarray.Dataset

        :Example:

        >>> su_index = climate_index.calculate_su()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.

        """
        summer_days = self.data[temp_var] > (25 + 273.15)  # 25°C in Kelvin
        summer_days_count = summer_days.groupby('time.year').sum(dim='time')
        summer_days_count.name = 'su'
        return xr.Dataset({'su': summer_days_count})

    def calculate_id(self, temp_var='t2m'):
        """
        Calculate the Number of Icing Days (ID) in a year.

        This method calculates the number of days in a year when the daily maximum temperature falls below 0°C (273.15 K).

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the number of icing days for each year.
        :rtype: xarray.Dataset

        :Example:

        >>> id_index = climate_index.calculate_id()

        :Warnings:
        An exception will be caught and a warning will be issued if the calculation cannot be performed.

        """
        icing_days = self.data[temp_var] < 273.15  # 0°C in Kelvin
        icing_days_count = icing_days.groupby('time.year').sum(dim='time')
        icing_days_count.name = 'id'
        return xr.Dataset({'id': icing_days_count})

    def calculate_tr(self, temp_var='t2m'):
        """
        Calculate the Number of Tropical Nights (TR) in a year.

        This method calculates the number of days in a year when the daily minimum temperature exceeds 20°C (293.15 K).

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the number of tropical nights for each year.
        :rtype: xarray.Dataset

        :Example:

        >>> tr_index = climate_index.calculate_tr()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.

        """
        tropical_nights = self.data[temp_var] > (20 + 273.15)  # 20°C in Kelvin
        tropical_nights_count = tropical_nights.groupby('time.year').sum(dim='time')
        tropical_nights_count.name = 'tr'
        return xr.Dataset({'tr': tropical_nights_count})

    def calculate_gsl(self, temp_var='t2m'):
        """
        Calculate the Length of the Growing Season (GSL) in days.

        This method calculates the length of the growing season based on daily mean temperature data.

        .. note:: 
        The input netCDF file should contain daily temperature data with a variable that can be accessed
        using the `temp_var` parameter.

        :param temp_var: The variable name in the netCDF file that contains daily temperature data. (default is 't2m')
        :type temp_var: str, optional

        :return: An xarray Dataset containing the length of the growing season in days.
        :rtype: xarray.Dataset

        :Example:

        >>> gsl_index = climate_index.calculate_gsl()

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
            Additional warnings will be issued if no growing season or end of growing season is found in the dataset.

        """
        try:
            # Resampling the daily mean temperature
            mean_temp = self.data[temp_var].resample(time='D').mean(dim='time')
            
            # Find the first run of days that start the growing season (mean temp > 5oC)
            growing_season_runs = first_run(mean_temp > (5 + 273.15), window=6)
            
            # Check if a growing season exists
            if growing_season_runs.size == 0:
                warnings.warn("No growing season found in the dataset.")
                return None
            
            # Take the first run as the start of the growing season
            growing_season_start = growing_season_runs[0]
            
            # Find the first run of days that end the growing season (mean temp < 5oC) after July 1st
            mean_temp_after_july = mean_temp.where(mean_temp['time.month'] > 6)
            end_season_runs = first_run(mean_temp_after_july < (5 + 273.15), window=6)
            
            # Check if an end of the growing season exists
            if end_season_runs.size == 0:
                warnings.warn("No end of the growing season found in the dataset.")
                return None
            
            # Take the first run as the end of the growing season
            end_season_start = end_season_runs[0]
            
            # Calculate the length of the growing season
            gsl = end_season_start - growing_season_start
            gsl.name = 'gsl'
            return xr.Dataset({'gsl': gsl})
            
        except Exception as e:
            warnings.warn(f"An error occurred while calculating GSL: {str(e)}")
            return None



    

    def calculate_percentile_p(self, precip_var='tp', percentile=95, wet_day_thresh_mm=1.0, reference_period=('1961-01-01', '1990-12-31'), resample_freq='Y'):
        """
        Calculate the Sum of Precipitation on Days Exceeding a Specified Percentile (Percentile Precipitation).

        This method calculates the sum of precipitation on days that exceed a given percentile for a specified period.

        .. note:: 
        The input netCDF file should contain daily precipitation data with a variable that can be accessed
        using the `precip_var` parameter.

        :param precip_var: The variable name in the netCDF file that contains daily precipitation data. (default is 'tp')
        :type precip_var: str, optional

        :param percentile: The percentile to use for identifying days with significant precipitation. (default is 95)
        :type percentile: int, optional

        :param wet_day_thresh_mm: The precipitation threshold in millimeters to consider a day as 'wet'. (default is 1.0)
        :type wet_day_thresh_mm: float, optional

        :param reference_period: The time period to use for calculating the percentile. (default is ('1961-01-01', '1990-12-31'))
        :type reference_period: tuple of str, optional

        :param resample_freq: The frequency for resampling the data. Common options include 'M' for monthly and 'Y' for yearly. (default is 'Y')
        :type resample_freq: str, optional

        :return: An xarray Dataset containing the sum of precipitation on days exceeding the percentile for each resampling period.
        :rtype: xarray.Dataset

        :Example:

        >>> rXXp_index_time_series = climate_index.calculate_percentile_p(precip_var='tp', percentile=95, wet_day_thresh_mm=1.0, reference_period=('1961-01-01', '1990-12-31'), resample_freq='Y')

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        try:
            pr = self.data[precip_var]  # Use the user-specified variable name

            # Convert the wet day threshold to meters
            wet_day_thresh_meters = wet_day_thresh_mm / 1000.0

            # Calculate the specified percentile for the reference period
            reference_period_data = pr.sel(time=slice(*reference_period))
            wet_days_reference = reference_period_data.where(reference_period_data >= wet_day_thresh_meters)
            p_reference = wet_days_reference.reduce(np.nanpercentile, q=percentile, dim='time')

            # Identify wet days
            wet_days = pr.where(pr >= wet_day_thresh_meters)

            # Compare each day with the corresponding percentile
            exceeding_days = wet_days > p_reference

            # Sum the precipitation for those days (RXXpTOT index)
            rXXp_index_time_series = (exceeding_days * pr).resample(time=resample_freq).sum()

            # Create a dataset to return
            rXXp_index_time_series.name = f'R{str(int(percentile))}pTOT'
            return xr.Dataset({f'R{str(int(percentile))}pTOT': rXXp_index_time_series})

        except Exception as e:
            warnings.warn(f"An error occurred while calculating R{str(int(percentile))}pTOT: {str(e)}")
            return None




    def calculate_precip_days(self, precip_var='tp', threshold_mm=10.0, resample_freq='Y'):
        """
        Count the Number of Days with Precipitation Above a Specified Threshold (Precipitation Days).

        This method counts the number of days in each resampling period where the precipitation amount exceeds a given threshold.

        .. note:: 
        The input netCDF file should contain daily precipitation data with a variable that can be accessed
        using the `precip_var` parameter.

        :param precip_var: The variable name in the netCDF file that contains daily precipitation data. (default is 'tp')
        :type precip_var: str, optional

        :param threshold_mm: The precipitation threshold in millimeters to consider a day as having significant precipitation. (default is 10.0)
        :type threshold_mm: float, optional

        :param resample_freq: The frequency for resampling the data. Common options include 'M' for monthly and 'Y' for yearly. (default is 'Y')
        :type resample_freq: str, optional

        :return: An xarray Dataset containing the count of days with precipitation above the threshold for each resampling period.
        :rtype: xarray.Dataset

        :Example:

        >>> precip_days_index = climate_index.calculate_precip_days(precip_var='tp', threshold_mm=10.0, resample_freq='Y')

        :Warnings:
        An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        try:
            # Convert the threshold to meters (assuming original data is in meters)
            threshold_meters = threshold_mm / 1000.0

            # Identify days with precipitation exceeding the threshold
            heavy_precip_days = self.data[precip_var].where(self.data[precip_var] >= threshold_meters)

            # Count the number of days exceeding the threshold per year (or any other resampling frequency)
            precip_days_count = heavy_precip_days.resample(time=resample_freq).count(dim='time')

            # Create a dataset to return
            precip_days_count.name = f'R{str(int(threshold_mm))}mm'
            return xr.Dataset({f'R{str(int(threshold_mm))}mm': precip_days_count})

        except Exception as e:
            warnings.warn(f"An error occurred while calculating R{str(int(threshold_mm))}mm: {str(e)}")
            return None

    
    def calculate_rxnday(self, precip_var='tp', n_days=5, resample_freq='M'):
        """
        Calculate Maximum Precipitation Sum for a Specified Number of Consecutive Days (RXnDay).
        
        This method finds the maximum sum of precipitation for a specified number of consecutive days,
        typically resampled over each month.

        .. note:: 
        The input netCDF file should contain daily precipitation data with a variable that can be accessed
        using the `precip_var` parameter.

        :param precip_var: The variable name in the netCDF file that contains daily precipitation data. (default is 'tp')
        :type precip_var: str, optional

        :param n_days: The number of consecutive days to consider for finding the maximum precipitation sum. (default is 5)
        :type n_days: int, optional

        :param resample_freq: The frequency for resampling the data. Common options include 'M' for monthly and 'Y' for yearly. (default is 'M')
        :type resample_freq: str, optional

        :return: An xarray Dataset containing the RXnDay values resampled according to the frequency.
        :rtype: xarray.Dataset

        :Example:

        >>> rxnday_index = climate_index.calculate_rxnday(precip_var='tp', n_days=5, resample_freq='M')

        :Warnings:
        An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        try:
            pr = self.data[precip_var]  # Use the user-specified variable name

            # Calculate the rolling n-day precipitation total
            rolling_nday_total = pr.rolling(time=n_days).sum()

            # Find the maximum n-day total for each time period (e.g., month)
            rxnday_index = rolling_nday_total.resample(time=resample_freq).max()

            # Convert to mm (assuming original data is in meters)
            rxnday_index_mm = rxnday_index * 1000

            # Create a dataset to return
            rxnday_index_mm.name = f'rx{str(n_days)}day'
            return xr.Dataset({f'rx{str(n_days)}day': rxnday_index_mm})

        except Exception as e:
            warnings.warn(f"An error occurred while calculating rx{str(n_days)}day: {str(e)}")
            return None
        


    # Maximum length of dry spell (CDD)
    def calculate_cdd(self, precip_var='tp'):
        """
        Calculate Consecutive Dry Days (CDD) based on daily precipitation data and a threshold.

        This method computes the CDD index, which represents the maximum length of dry spells,
        defined as periods with daily precipitation less than 1 mm, in a year.

        .. note:: 
        The input netCDF file should contain daily precipitation data with a variable that can be accessed
        using the `precip_var` parameter.

        :param precip_var: The variable name in the netCDF file that contains daily precipitation data. (default is 'tp')
        :type precip_var: str, optional
        
        :return: An xarray Dataset containing the CDD values for each year.
        :rtype: xarray.Dataset

        :Example:

        >>> cdd_index = climate_index.calculate_cdd(precip_var='tp')

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        try:
            pr = self.data[precip_var]  # Daily precipitation amount
            dry_days = pr < 1.0 / 1000  # Convert 1mm to meters
            cdd = dry_days.groupby('time.year').apply(longest_run, condition_func=lambda x: x)
            cdd.name = 'CDD'
            return xr.Dataset({'CDD': cdd})
        except Exception as e:
            warnings.warn(f"An error occurred while calculating CDD: {str(e)}")
            return None

    # Maximum length of wet spell (CWD)
    def calculate_cwd(self, precip_var='tp'):
        """
        Calculate Consecutive Wet Days (CWD) based on daily precipitation data and a threshold.

        This method computes the CWD index, which represents the maximum length of wet spells,
        defined as periods with daily precipitation equal to or greater than 1 mm, in a year.

        .. note:: 
        The input netCDF file should contain daily precipitation data with a variable that can be accessed
        using the `precip_var` parameter.

        :param precip_var: The variable name in the netCDF file that contains daily precipitation data. (default is 'tp')
        :type precip_var: str, optional
        
        :return: An xarray Dataset containing the CWD values for each year.
        :rtype: xarray.Dataset

        :Example:

        >>> cwd_index = climate_index.calculate_cwd(precip_var='tp')

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        try:
            pr = self.data[precip_var]  # Daily precipitation amount
            wet_days = pr >= 1.0 / 1000  # Convert 1mm to meters
            cwd = wet_days.groupby('time.year').apply(longest_run, condition_func=lambda x: x)
            cwd.name = 'CWD'
            return xr.Dataset({'CWD': cwd})
        except Exception as e:
            warnings.warn(f"An error occurred while calculating CWD: {str(e)}")
            return None

    # Simple Precipitation Intensity Index (SDII)
    def calculate_sdii(self, precip_var='tp'):
        """
        Calculate Simple Precipitation Intensity Index (SDII) based on daily precipitation data.

        This method computes the SDII index, which represents the mean precipitation amount on wet days
        (days with precipitation equal to or greater than 1 mm) in a year.

        .. note:: 
        The input netCDF file should contain daily precipitation data with a variable that can be accessed
        using the `precip_var` parameter.

        :param precip_var: The variable name in the netCDF file that contains daily precipitation data. (default is 'tp')
        :type precip_var: str, optional
        
        :return: An xarray Dataset containing the SDII values for each year.
        :rtype: xarray.Dataset

        :Example:

        >>> sdii_index = climate_index.calculate_sdii(precip_var='tp')

        :Warnings:
            An exception will be caught and a warning will be issued if the calculation cannot be performed.
        """
        try:
            pr = self.data[precip_var]  # Daily precipitation amount
            wet_days = pr.where(pr >= 1.0 / 1000)  # Select wet days and convert 1mm to meters
            wet_days_count = wet_days.groupby('time.year').count(dim='time')
            wet_days_sum = wet_days.groupby('time.year').sum(dim='time')
            sdii = wet_days_sum / wet_days_count
            sdii.name = 'SDII'
            return xr.Dataset({'SDII': sdii})
        except Exception as e:
            warnings.warn(f"An error occurred while calculating SDII: {str(e)}")
            return None        

    
    def test_indices(self, data_type='temperature'):
        # Methods to test for temperature data
        temp_methods_to_test = [
            self.calculate_tx10p,
            self.calculate_tx90p,
            self.calculate_tn10p,
            self.calculate_tn90p,
            self.calculate_frost_days,
            self.calculate_warm_spell,
            self.calculate_cold_spell,
            self.calculate_dtr,
            self.calculate_su,
            self.calculate_id,
            self.calculate_tr,
            self.calculate_gsl
        ]

        # Methods to test for precipitation data
        precip_methods_to_test = [
            self.calculate_percentile_p,
            self.calculate_precip_days,
            self.calculate_rxnday, 
            self.calculate_cdd,
            self.calculate_cwd,
            self.calculate_sdii
        ]

        # Choose the methods to test based on the data type
        if data_type == 'temperature':
            methods_to_test = temp_methods_to_test
        elif data_type == 'precipitation':
            methods_to_test = precip_methods_to_test
        else:
            raise ValueError("Invalid data type. Choose 'temperature' or 'precipitation'.")

        # Suppress warnings during testing
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Iterate through the methods and print success or error message
            for method in methods_to_test:
                try:
                    print(f"Testing {method.__name__}...")
                    _ = method()
                    print(f"Performed successful test for {method.__name__}!")
                except Exception as e:
                    print(f"An error occurred while testing {method.__name__}: {str(e)}")
                    continue