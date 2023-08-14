import xarray as xr
import numpy as np
from xclim.indices import tx10p, tx90p, tn10p, tn90p, warm_spell_duration_index, cold_spell_duration_index, frost_days
from xclim.core.calendar import percentile_doy
import warnings

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
        tasmax_per = percentile_doy(self.data[temp_var], per=10, window=5).sel(percentiles=10)
        return tx10p(self.data[temp_var], tasmax_per)

    def calculate_tx90p(self, temp_var='t2m'):
        tasmax_per = percentile_doy(self.data[temp_var], per=90, window=5).sel(percentiles=90)
        return tx90p(self.data[temp_var], tasmax_per)

    def calculate_tn10p(self, temp_var='t2m'):
        tasmin_per = percentile_doy(self.data[temp_var], per=10, window=5).sel(percentiles=10)
        return tn10p(self.data[temp_var], tasmin_per)

    def calculate_tn90p(self, temp_var='t2m'):
        tasmin_per = percentile_doy(self.data[temp_var], per=90, window=5).sel(percentiles=90)
        return tn90p(self.data[temp_var], tasmin_per)

    def calculate_frost_days(self, temp_var='t2m'):
        return frost_days(self.data[temp_var], thresh='273.15 K')

    def calculate_warm_spell(self, temp_var='t2m'):
        tasmax_per = percentile_doy(self.data[temp_var], per=90, window=5).sel(percentiles=90)
        return warm_spell_duration_index(self.data[temp_var], tasmax_per)

    def calculate_cold_spell(self, temp_var='t2m'):
        tasmin_per = percentile_doy(self.data[temp_var], per=10, window=5).sel(percentiles=10)
        return cold_spell_duration_index(self.data[temp_var], tasmin_per)
    

    def calculate_percentile_p(self, precip_var='tp', percentile=95, wet_day_thresh_mm=1.0, reference_period=('1961-01-01', '1990-12-31'), resample_freq='Y'):
        pr = self.data[precip_var] # Use the user-specified variable name
        # Exclude leap day (February 29th) if present
        pr = pr.where((pr['time.month'] != 2) | (pr['time.day'] != 29), drop=True)

        # Convert the wet day threshold to meters
        wet_day_thresh_meters = wet_day_thresh_mm / 1000.0

        # Calculate the specified percentile for the reference period
        reference_period_data = pr.sel(time=slice(*reference_period))
        wet_days_reference = reference_period_data.where(reference_period_data >= wet_day_thresh_meters)
        p_reference = wet_days_reference.groupby('time.dayofyear').reduce(np.nanpercentile, q=percentile, dim='time')

        # Identify wet days
        wet_days = pr.where(pr >= wet_day_thresh_meters)

        # Compare each day with the corresponding percentile
        exceeding_days = wet_days.groupby('time.dayofyear') > p_reference

        # Sum the precipitation for those days (rXXp index)
        rXXp_index_time_series = (exceeding_days * pr).resample(time=resample_freq).sum()

        # Count the number of days exceeding the specified percentile
        rXXp_days_count = exceeding_days.resample(time=resample_freq).sum(dim='time')

        return rXXp_index_time_series * 1000, rXXp_days_count



    def calculate_precip_days(self, precip_var='tp', threshold_mm=10.0):
        # Convert the threshold to meters
        threshold_meters = threshold_mm / 1000.0

        # Identify days with precipitation exceeding the threshold
        heavy_precip_days = self.data['tp'].where(self.data['tp'] >= threshold_meters)

        # Count the number of days exceeding the threshold per year
        precip_days_index = heavy_precip_days.groupby('time.year').count(dim='time')

        return precip_days_index
    
    def calculate_rxnday(self, precip_var='tp', n_days=5, threshold_mm=50.0, resample_freq='Y'):
        pr = self.data[precip_var]  # Use the user-specified variable name

        # Calculate the rolling n-day precipitation total
        rolling_nday_total = pr.rolling(time=n_days).sum()

        # Find the maximum n-day total for each time period (e.g., year)
        rxnday_index = rolling_nday_total.resample(time=resample_freq).max()

        # Count the number of n-day heavy precipitation periods
        heavy_precip_periods = rolling_nday_total.where(rolling_nday_total >= threshold_mm / 1000.0)
        heavy_precip_count = heavy_precip_periods.resample(time=resample_freq).count(dim='time')

        return rxnday_index * 1000, heavy_precip_count
    
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
        ]

        # Methods to test for precipitation data
        precip_methods_to_test = [
            self.calculate_percentile_p,
            self.calculate_precip_days,
            self.calculate_rxnday
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

