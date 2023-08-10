
# ClimateIndex calculator

The `ClimateIndex` class provides a systematic and convenient interface for preprocessing climate data and calculating various climate indices. It can handle different climate variables and allows for the selection, resampling, and filling of missing values according to user-defined parameters. It is currently designed to calculate indices such as `tx10p`, `tx90p`, `tn10p`, `tn90p`, `warm spell duration`, `cold spell duration`, and `frost days`.

# Class Definition

## Initialization
To initialize the ClimateIndex class, you must provide the path to the NetCDF data file containing the climate data you wish to process.

``` {.python language="Python"}
climate_index = ClimateIndex(datafile='/path/to/your/datafile.nc')
```

-   `datafile` (str): Path to the NetCDF file containing the climate
    data.

## Preprocessing Method

``` {.python language="Python"}
climate_index.pre_process(time_range=None, fill_missing=None, resample_freq=None, months=None, season=None, resample_method='mean')
```

-   `time_range` (tuple, optional): This parameter allows you to specify a start and end time for the data you want to analyze. By providing a tuple with two date strings, you can subset the dataset to include only the observations within that time frame.
    -   Example: `time_range=('1961-01-01', '1990-12-31')` will include only the data from January 1, 1961, to December 31, 1990.
    -   Use Case: This can be useful when analyzing climate patterns over specific periods or when comparing different decades.

-   `fill_missing` (int, float, or 'interpolate', optional): This parameter provides options to handle missing values (NaN) in the dataset.
    - If a numerical value is provided, all missing values will be replaced with that number.
    - If `interpolate` is specified, linear interpolation will be used to estimate missing values based on neighbouring observations.
    - Example: `fill_missing=0` will replace all NaN values with 0.
    - Use Case: Handling missing values is crucial in climate analysis to avoid bias and errors in subsequent calculations.

-   `resample_freq` (str, optional): This parameter defines the frequency at which the data should be resampled. It allows you to change the time-frequency of the dataset.
    - Example: `resample_freq='MS'` will resample the data to monthly frequency, starting at the beginning of each month.
    - Use Case: Resampling is useful when analysing data at a different temporal resolution, such as monthly or yearly.

-   `months` (list, optional): This parameter allows you to select specific months from the dataset.
    - Example: `months=[1, 2, 3]` will include only the data from January, February, and March.
    - Use Case: Useful for seasonal analysis or when studying climate behaviour during particular months.

-   `season` (str, optional): This parameter enables the selection of data for a specific meteorological season.
    - Example: season='JJA' will include only the data from the summer months (June, July, August).
    - Use Case: Helpful for analyzing seasonal patterns and variations in climate.

-   `resample_method` (str, optional, default='mean'): This parameter defines the method used for resampling when changing the time-frequency with `resample_freq`.
    - Default is `mean`, meaning that the mean value of the observations within each resampling period will be used.
    - Other methods, such as `sum`, `max`, etc., can be specified if needed.
    - Use Case: This provides flexibility in how the data is aggregated when resampling, allowing for different types of analysis and interpretation.

## Climate Indices Calculation Methods

### tx10p (Temperature Exceeding 10th Percentile)

``` {.python language="Python"}
tx10p_index = climate_index.calculate_tx10p()
```
- Description: The tx10p index calculates the number of days when the maximum temperature exceeds the 10th percentile of a reference period.
- Calculation: It's computed by first calculating the 10th percentile of daily maximum temperature over a reference period (e.g., 1961-1990) and then counting how many days in the target period exceed this threshold.


### tx90p (Temperature Exceeding 90th Percentile)

``` {.python language="Python"}
tx90p_index = climate_index.calculate_tx90p()
```
- Description: The tx90p index calculates the number of days when the maximum temperature exceeds the 90th percentile of a reference period.
- Calculation: Similar to tx10p, but using the 90th percentile as the threshold.

### tn10p (Temperature Below 10th Percentile)

``` {.python language="Python"}
tn10p_index = climate_index.calculate_tn10p()
```
- Description: The tn10p index calculates the number of days when the minimum temperature is below the 10th percentile of a reference period.
- Calculation: It's computed by first calculating the 10th percentile of daily minimum temperature over a reference period and then counting how many days in the target period fall below this threshold.

### tn90p (Temperature Exceeding 90th Percentile)

``` {.python language="Python"}
tn90p_index = climate_index.calculate_tn90p()
```
- Description:  The tn90p index calculates the number of days when the minimum temperature exceeds the 90th percentile of a reference period.
- Calculation: Similar to tn10p, but using the 90th percentile as the threshold.

### Frost Days

``` {.python language="Python"}
frost_days_index = climate_index.calculate_frost_days()
```
- Description:  The frost days index counts the number of days when the minimum temperature falls below 0°C (or another specified threshold).
- Calculation:  The number of days in the target period when the daily minimum temperature falls below the freezing point.

    
### Warm Spell Duration Index

``` {.python language="Python"}
warm_spell_index = climate_index.calculate_warm_spell()
```
- Description:  This index measures the duration of warm spells, defined as periods when the daily maximum temperature exceeds the 90th percentile for at least a specified number of consecutive days.
- Calculation: Length of consecutive periods when daily maximum temperature exceeds the 90th percentile of a reference period.

  ### Cold Spell Duration Index

``` {.python language="Python"}
cold_spell_index = climate_index.calculate_cold_spell()
```
- Description: This index measures the duration of cold spells, defined as periods when the daily minimum temperature is below the 10th percentile for at least a specified number of consecutive days.
- Calculation: Length of consecutive periods when daily minimum temperature falls below the 10th percentile of a reference period.

These indices provide various ways to quantify temperature extremes and their occurrence, helping in understanding and analyzing climate patterns, trends, and variability.

# Example Usage

``` {.python language="Python"}
climate_index = ClimateIndex(datafile='/path/to/your/datafile.nc')
climate_index.pre_process()
tx10p_index = climate_index.calculate_tx10p()

def plot_tx10p(data, title):
    # Take the mean along latitude and longitude
    data_mean = data.mean(dim=['latitude', 'longitude'])
    
    data_mean.plot(figsize=(10, 5))
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('tx10p Index')
    plt.grid(True)
    plt.show()

plot_tx10p(tx10p_index, "tx10p test")

```
This example demonstrates how to initialize the ClimateIndex class with a specified data file, preprocess the data, calculate the tx10p index, and then plot the resulting time series after taking the mean along latitude and longitude.

