import xarray as xr
import numpy as np
from xclim.indices import tx10p, tx90p, tn10p, tn90p, warm_spell_duration_index, cold_spell_duration_index, frost_days
from xclim.core.calendar import percentile_doy
import warnings

def first_run(data, window, dim='time'):
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
    def __init__(self, datafile):
        self.datafile = datafile
        self.data = None

    def pre_process(self, time_range=None, fill_missing=None, resample_freq=None, months=None, season=None, resample_method='mean'):
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
        tasmax_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=10, window=5).sel(percentiles=10)
        result = self.data[temp_var].where(self.data[temp_var] < tasmax_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tx10p"
        return xr.Dataset({"tx10p": percentage_days})

    def calculate_tx90p(self, temp_var='t2m'):
        tasmax_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=90, window=5).sel(percentiles=90)
        result = self.data[temp_var].where(self.data[temp_var] > tasmax_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tx90p"
        return xr.Dataset({"tx90p": percentage_days})

    def calculate_tn10p(self, temp_var='t2m'):
        tasmin_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=10, window=5).sel(percentiles=10)
        result = self.data[temp_var].where(self.data[temp_var] < tasmin_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tn10p"
        return xr.Dataset({"tn10p": percentage_days})

    def calculate_tn90p(self, temp_var='t2m'):
        tasmin_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=90, window=5).sel(percentiles=90)
        result = self.data[temp_var].where(self.data[temp_var] > tasmin_per)
        percentage_days = (result.count(dim="time") / self.data[temp_var].count(dim="time")) * 100
        percentage_days.name = "tn90p"
        return xr.Dataset({"tn90p": percentage_days})


    def calculate_frost_days(self, temp_var='t2m'):
        frost_days_dat = frost_days(self.data[temp_var], thresh='273.15 K')
        frost_days_dat.name = 'frost_days_dat'
        return xr.Dataset({'frost_days_dat': frost_days_dat})

    def calculate_warm_spell(self, temp_var='t2m'):
        tasmax_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=90, window=5).sel(percentiles=90)
        warm_days = self.data[temp_var] > tasmax_per
        warm_spell_days = first_run(warm_days, window=6, dim='time')
        warm_spell_days.name = 'wsdi'
        return xr.Dataset({'wsdi': warm_spell_days})

    def calculate_cold_spell(self, temp_var='t2m'):
        tasmin_per = percentile_doy(self.data[temp_var].sel(time=slice('1961-01-01', '1990-12-31')), per=10, window=5).sel(percentiles=10)
        cold_days = self.data[temp_var] < tasmin_per
        cold_spell_days = first_run(cold_days, window=6, dim='time')
        cold_spell_days.name = 'csdi'
        return xr.Dataset({'csdi': cold_spell_days})

    
    def calculate_dtr(self, temp_var='t2m'):
        daily_max = self.data[temp_var].resample(time='D').max(dim='time')
        daily_min = self.data[temp_var].resample(time='D').min(dim='time')
        dtr = daily_max - daily_min
        dtr_monthly_mean = dtr.resample(time='M').mean(dim='time')
        dtr_monthly_mean.name = 'dtr'
        return xr.Dataset({'dtr': dtr_monthly_mean})

    def calculate_su(self, temp_var='t2m'):
        summer_days = self.data[temp_var] > (25 + 273.15)  # 25°C in Kelvin
        summer_days_count = summer_days.groupby('time.year').sum(dim='time')
        summer_days_count.name = 'su'
        return xr.Dataset({'su': summer_days_count})

    def calculate_id(self, temp_var='t2m'):
        icing_days = self.data[temp_var] < 273.15  # 0°C in Kelvin
        icing_days_count = icing_days.groupby('time.year').sum(dim='time')
        icing_days_count.name = 'id'
        return xr.Dataset({'id': icing_days_count})

    def calculate_tr(self, temp_var='t2m'):
        tropical_nights = self.data[temp_var] > (20 + 273.15)  # 20°C in Kelvin
        tropical_nights_count = tropical_nights.groupby('time.year').sum(dim='time')
        tropical_nights_count.name = 'tr'
        return xr.Dataset({'tr': tropical_nights_count})

    def calculate_gsl(self, temp_var='t2m'):
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