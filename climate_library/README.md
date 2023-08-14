
# ClimateIndex Library

# ClimateIndex Class

The `ClimateIndex` class provides a systematic and convenient interface for preprocessing climate data and calculating a wide range of climate indices. It can handle different climate variables, including temperature and precipitation, and allows for the selection, resampling, and filling of missing values according to user-defined parameters.

## Temperature-related Indices

For temperature-related analyses, it can calculate indices such as:

- **tx10p**: The percentage of days when daily maximum temperature is below the 10th percentile
- **tx90p**: The percentage of days when daily maximum temperature is above the 90th percentile
- **tn10p**: The percentage of days when daily minimum temperature is below the 10th percentile
- **tn90p**: The percentage of days when daily minimum temperature is above the 90th percentile
- **frost_days**: Number of frost days (when daily minimum temperature is below 0°C)
- **warm_spell_duration_index**: Warm spell duration index
- **cold_spell_duration_index**: Cold spell duration index

## Precipitation-related Indices

For precipitation-related analyses, it can calculate indices such as:

- **percentile_p**: The sum of precipitation in wet days exceeding a given percentile during a time period
- **precip_days**: The number of days with precipitation exceeding a specified threshold
- **rxnday**: The highest n-day total precipitation amount for a given resampling frequency

## Testing

The class also includes a testing method to ensure all calculations are performed correctly. The modularity of the class allows for easy expansion to include additional climate indices in the future.


## Installation

### From PyPI (when available)

You can install the package using pip:

```bash
pip install climate_library
```

### From GitHub

Alternatively, you can install the library directly from the GitHub repository:

```bash
pip install git+https://github.com/shiv3679/climate-indices.git
```

### Dependencies

The following packages are required dependencies and will be installed automatically:

- xarray
- xclim
- numpy

### Compatibility

The library is compatible with Python 3.6 and higher.


## Usage

### Initialization

Create an instance of the `ClimateIndex` class by providing the path to your climate data file. This class is the main entry point for working with climate indices in this library.

#### Example

```python
from climate_library.climate_index import ClimateIndex

# Path to the NetCDF file containing climate data
datafile = 'path/to/your/datafile.nc'

# Create an instance of the ClimateIndex class
climate_index = ClimateIndex(datafile)
```

The `datafile` parameter should be the path to a NetCDF file that contains the climate data you want to analyze. The file must include dimensions for time, latitude, and longitude, and it must contain the relevant climate variables (e.g., temperature, precipitation) that you wish to use in calculating climate indices.

Once you have created an instance of the `ClimateIndex` class, you can use its methods to preprocess the data and calculate various climate indices, as described in the following sections.


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

#### Example Usage

```python
climate_index.pre_process(time_range=('2000-01-01', '2020-12-31'), resample_freq='M', fill_missing='interpolate')
```

This example demonstrates how to preprocess the data by selecting a time range, resampling to monthly frequency, and interpolating missing values. Adjust these parameters according to your specific needs and analysis goals.

### Climate Indices Calculation Methods

#### Temperature Indices

##### TX10P (Temperature Exceeding 10th Percentile)


```python
tx10p_index = climate_index.calculate_tx10p()
```
- **Description**: Calculates the number of days when the maximum temperature exceeds the 10th percentile of a reference period.
- **Calculation**: Computed by first calculating the 10th percentile of daily maximum temperature over a reference period and then counting how many days in the target period exceed this threshold.

##### TX90P (Temperature Exceeding 90th Percentile)


```python
tx90p_index = climate_index.calculate_tx90p()
```
- **Description**: Calculates the number of days when the maximum temperature exceeds the 90th percentile of a reference period.
- **Calculation**: Similar to TX10P, but using the 90th percentile as the threshold.

##### TN10P (Temperature Below 10th Percentile)

```python
tn10p_index = climate_index.calculate_tn10p()
```
- **Description**: Calculates the number of days when the minimum temperature is below the 10th percentile of a reference period.
- **Calculation**: Similar to TX10P, but for daily minimum temperature.

##### TN90P (Temperature Exceeding 90th Percentile)

```python
tn90p_index = climate_index.calculate_tn90p()
```
- **Description**: Calculates the number of days when the minimum temperature exceeds the 90th percentile of a reference period.
- **Calculation**: Similar to TN10P, but using the 90th percentile as the threshold.

##### Frost Days


```python
frost_days_index = climate_index.calculate_frost_days()
```
- **Description**: Counts the number of days when the minimum temperature falls below 0°C.
- **Calculation**: Number of days when daily minimum temperature falls below freezing.

##### Warm Spell Duration Index


```python
warm_spell_index = climate_index.calculate_warm_spell()
```
- **Description**: Measures the duration of warm spells when daily maximum temperature exceeds the 90th percentile for consecutive days.
- **Calculation**: Length of periods exceeding the 90th percentile of maximum temperature.

##### Cold Spell Duration Index


```python
cold_spell_index = climate_index.calculate_cold_spell()
```
- **Description**: Measures the duration of cold spells when daily minimum temperature falls below the 10th percentile for consecutive days.
- **Calculation**: Length of periods falling below the 10th percentile of minimum temperature.

#### Precipitation Indices

##### Percentile Precipitation


```python
rXXp_index_time_series, rXXp_days_count = climate_index.calculate_percentile_p(precip_var='tp', percentile=95, wet_day_thresh_mm=1.0, reference_period=('1961-01-01', '1990-12-31'), resample_freq='Y')
```
- **Description**: Calculates the sum of precipitation on days exceeding a specified percentile.
- **Calculation**: Sum of precipitation on days exceeding the threshold defined by the percentile.
- **Example**: To calculate the sum of precipitation on days exceeding the 95th percentile and the number of such days, with a wet day threshold of 1 mm, and reference period from 1961 to 1990, use percentile=95, wet_day_thresh_mm=1.0, reference_period=('1961-01-01', '1990-12-31').


##### Precipitation Days


```python
precip_days_index = climate_index.calculate_precip_days(precip_var='tp', threshold_mm=10.0)
```
- **Description**: Counts the number of days with precipitation above a specified threshold.
- **Calculation**: Number of days with precipitation exceeding the threshold.
- **Example**: To calculate the number of days with precipitation above 10 mm, use threshold_mm=10.0.

##### RXnDay

```python
rxnday_index, heavy_precip_count = climate_index.calculate_rxnday(precip_var='tp', n_days=5, threshold_mm=50.0, resample_freq='Y')
```
- **Description**: Finds the maximum precipitation sum for a specified number of consecutive days.
- **Calculation**: Maximum sum of precipitation for the defined number of consecutive days.
- **Example**: To calculate the maximum 5-day consecutive precipitation sum and count the number of periods with more than 50 mm, use n_days=5, threshold_mm=50.0.

These indices provide various ways to quantify temperature and precipitation extremes, trends, and variability, aiding in climate analysis, pattern recognition, and decision-making.

### Testing Indices

The `test_indices` method provides a way to test the implemented indices for a specific data type (either 'temperature' or 'precipitation').

```python
climate_index.test_indices(data_type='temperature')
```

#### Importance

- **Validation**: Tests the implemented indices to ensure that they are functioning correctly.
- **Debugging**: Helps in identifying and resolving any errors or inconsistencies.
- **Flexibility**: Allows testing for different data types (temperature or precipitation) based on the user's needs.

#### Functionality

- **Temperature Testing**: If `data_type='temperature'`, the following methods will be tested:
  - `calculate_tx10p`
  - `calculate_tx90p`
  - `calculate_tn10p`
  - `calculate_tn90p`
  - `calculate_frost_days`
  - `calculate_warm_spell`
  - `calculate_cold_spell`

- **Precipitation Testing**: If `data_type='precipitation'`, the following methods will be tested:
  - `calculate_percentile_p`
  - `calculate_precip_days`
  - `calculate_rxnday`

#### Example Usage

- To test all temperature-related indices:
  ```python
  climate_index.test_indices(data_type='temperature')
  ```

- To test all precipitation-related indices:
  ```python
  climate_index.test_indices(data_type='precipitation')
  ```

## Authors

### Shiv Shankar Singh
- **Role**: Lead Developer
- **Affiliation**: [IISER Mohali](https://www.iisermohali.ac.in/)
- **Contributions**: Developed the core algorithms and preprocessing methods, conducted testing, and contributed to the documentation.
- **Email**: shivshankarsingh.py@gmail.com
- **GitHub**: [shiv3679](https://github.com/shiv3679)

### Akshit Nanda
- **Role**: Co-Developer
- **Affiliation**: [IISER Mohali](https://www.iisermohali.ac.in/)
- **Contributions**: Contributed to the development of precipitation indices and data analysis
- **Email**: nandaakshit2000@gmail.com


The authors would like to acknowledge the support of Dr. Raju Attada and [Sree Hari](sreeharik@iisermohali.ac.in) from the WEATHER AND CLIMATE MODELLING RESEARCH GROUP (IISER MOHALI) in the development of this library.

### Contact

For any inquiries or collaboration, please feel free to reach out to the authors via the provided email addresses.


## License

This project is licensed under the MIT License.