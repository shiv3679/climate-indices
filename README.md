# Climate Library: A Python Library for Climate Indices

<!-- Badges section -->
[![PyPI Version](https://img.shields.io/pypi/v/climate-library.svg)](https://pypi.org/project/climate-library/)
[![Conda Forge Version](https://img.shields.io/conda/vn/conda-forge/climate-library.svg)](https://anaconda.org/conda-forge/climate-library)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/climate-library.svg)](https://pypi.org/project/climate-library/)

<!-- Documentation and Support -->
[![Documentation Status](https://readthedocs.org/projects/climate-library/badge/?version=latest)](https://climate-library.readthedocs.io/)
[![Gitter Chat](https://badges.gitter.im/climate-library/community.svg)](https://gitter.im/climate-library/community)

<!-- Open Source and Licensing -->
[![License](https://img.shields.io/github/license/yourusername/climate-library.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/xxx/xxx.svg)](https://doi.org/xxx/xxx)

<!-- Coding Standards -->
[![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

<!-- Development Status -->
[![Build Status](https://travis-ci.com/yourusername/climate-library.svg?branch=master)](https://travis-ci.com/yourusername/climate-library)
[![Coverage Status](https://coveralls.io/repos/github/yourusername/climate-library/badge.svg?branch=master)](https://coveralls.io/github/yourusername/climate-library?branch=master)

<!-- Project Description -->

# ClimateIndex Class

The `ClimateIndex` class provides a systematic and convenient interface for preprocessing climate data and calculating a wide range of climate indices. It can handle different climate variables, including temperature and precipitation, and allows for the selection, resampling, and filling of missing values according to user-defined parameters.

## Temperature-related Indices

For temperature-related analyses, it can calculate indices such as:

- **tx10p**: The percentage of days when daily maximum temperature is below the 10th percentile
- **tx90p**: The percentage of days when daily maximum temperature is above the 90th percentile
- **tn10p**: The percentage of days when daily minimum temperature is below the 10th percentile
- **tn90p**: The percentage of days when daily minimum temperature is above the 90th percentile
- **frost_days**: Number of frost days (when daily minimum temperature is below 0°C)
- **warm_spell_duration_index**: The maximum length of a run of days with daily maximum temperature above the 90th percentile
- **cold_spell_duration_index**: The maximum length of a run of days with daily minimum temperature below the 10th percentile.
- **Diurnal Temperature Range (DTR)**: The mean monthly diurnal temperature range.
- **Summer Days (SU)**: The annual count of days when daily maximum temperature is above 25°C.
- **Icing Days (ID)**: The annual count of days when daily maximum temperature is below 0°C.
- **Tropical Nights (TR)**: The annual count of days when daily minimum temperature is above 20°C.
- **Growing Season Length (GSL)**: The length of the growing season in day

## Precipitation-related Indices

For precipitation-related analyses, it can calculate indices such as:

- **percentile_p**: The sum of precipitation in wet days exceeding a given percentile during a time period
- **precip_days**: The number of days with precipitation exceeding a specified threshold
- **rxnday**: The highest n-day total precipitation amount for a given resampling frequency
- **Consecutive Dry Days (CDD)**: The maximum length of dry spell in days.
- **Consecutive Wet Days (CWD)**: The maximum length of wet spell in days.
- **Simple Precipitation Intensity Index (SDII)**: The annual mean precipitation intensity on wet days.

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

#### Additional Checks 
- **Units Check**: The function also performs a check on the unit of the data variables. It supports 'K' for temperature and 'm' for precipitation.
- **Special Rule**: If season is 'DJF', and `resample_freq` is not specified, it defaults to 'QS-DEC'.

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

##### Icing Days


```python
id_index = climate_index.calculate_id()
```
- **Description**: Counts the number of days when the maximum temperature falls below 0°C.
- **Calculation**: Number of days when daily maximum temperature is below freezing (273.15 K).


##### Summer Days


```python
su_index = climate_index.calculate_su()
```
- **Description**: Counts the number of days when the maximum temperature exceeds 25°C.
- **Calculation**: Number of days when daily maximum temperature is above 25°C (298.15 K).

##### Tropical Nights


```python
tr_index = climate_index.calculate_tr()
```
- **Description**: Counts the number of days when the minimum temperature exceeds 20°C.
- **Calculation**: Number of days when daily minimum temperature is above 20°C (293.15 K).

##### Growing Season Length (GSL)


```python
gsl_index = climate_index.calculate_gsl()
```
- **Description**: Calculates the length of the growing season in days.
- **Calculation**: The growing season starts with the first run of at least 6 days where the daily mean temperature exceeds 5°C and ends with the first run of at least 6 days where the daily mean temperature falls below 5°C after July 1st.

##### Diurnal Temperature Range (DTR)


```python
dtr_index = climate_index.calculate_dtr()
```
- **Description**: Calculates the mean monthly diurnal temperature range.
- **Calculation**: The difference between the daily maximum and minimum temperature is averaged over each month.

#### Precipitation Indices

##### Percentile Precipitation


```python
rXXp_index_time_series = climate_index.calculate_percentile_p(precip_var='tp', percentile=95, wet_day_thresh_mm=1.0, reference_period=('1961-01-01', '1990-12-31'), resample_freq='Y')
```
- **Description**: Calculates the sum of precipitation on days exceeding a specified percentile.
- **Calculation**: Sum of precipitation on days exceeding the threshold defined by the percentile.
- **Example**: To calculate the sum of precipitation on days exceeding the 95th percentile, with a wet day threshold of 1 mm, and reference period from 1961 to 1990, use percentile=95, wet_day_thresh_mm=1.0, and reference_period=('1961-01-01', '1990-12-31').


##### Precipitation Days


```python
precip_days_index = climate_index.calculate_precip_days(precip_var='tp', threshold_mm=10.0, resample_freq='Y')
```
- **Description**: Counts the number of days with precipitation above a specified threshold.
- **Calculation**: Number of days with precipitation exceeding the threshold.
- **Example**: To calculate the number of days with precipitation above 10 mm, use threshold_mm=10.0.

##### RXnDay

```python
rxnday_index = climate_index.calculate_rxnday(precip_var='tp', n_days=5, resample_freq='M')
```
- **Description**: Finds the maximum precipitation sum for a specified number of consecutive days.
- **Calculation**: Maximum sum of precipitation for the defined number of consecutive days.
- **Example**: To calculate the maximum 5-day consecutive precipitation sum, use n_days=5.

##### Consecutive Wet Days (CWD)

```python
cwd_index = climate_index.calculate_cwd(precip_var='tp')
```
- **Description**: Calculates the maximum length of consecutive wet days in a year.
- **Calculation**:Length of the longest sequence of wet days (precipitation equal to or greater than 1 mm).

##### Consecutive Dry Days (CDD)

```python
cdd_index = climate_index.calculate_cdd(precip_var='tp')
```
- **Description**: Calculates the maximum length of consecutive dry days in a year.
- **Calculation**:Length of the longest sequence of dry days (precipitation less than 1 mm).

##### Simple Precipitation Intensity Index (SDII)

```python
sdii_index = climate_index.calculate_sdii(precip_var='tp')
```
- **Description**: Calculates the average precipitation intensity on wet days.
- **Calculation**:Sum of precipitation on wet days (≥ 1 mm) divided by the total number of wet days in each year.


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
  - `calculate_dtr`
  - `calculate_su`
  - `calculate_id`
  - `calculate_tr`
  - `calculate_gsl`

- **Precipitation Testing**: If `data_type='precipitation'`, the following methods will be tested:
  - `calculate_percentile_p`
  - `calculate_precip_days`
  - `calculate_rxnday`
  - `calculate_cdd`
  - `calculate_cwd`
  - `calculate_sdii`

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
