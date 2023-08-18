import os
import xarray as xr
import warnings
from climate_library.climate_index import ClimateIndex

def test_pre_process():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Path to the test data file
        datafile = os.path.join(os.path.dirname(__file__), 'data', 'cir.nc')
        
        # Create an instance of the ClimateIndex class
        ci = ClimateIndex(datafile)

        # Test the pre_process method
        ci.pre_process(resample_freq='M')
        assert ci.data is not None, "Data loading failed."
        assert 'time' in ci.data.dims, "Time dimension not found in the dataset."
        assert len(ci.data['time']) == 132, "Time selection or resampling failed."

def test_selected_temperature_indices():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Path to the test data file
        datafile = os.path.join(os.path.dirname(__file__), 'data', 'cir.nc')
        
        # Create an instance of the ClimateIndex class
        ci = ClimateIndex(datafile)
        ci.pre_process()

        # Test the selected temperature indices
        assert isinstance(ci.calculate_tx10p(), xr.DataArray), "Failed to calculate tx10p."
        assert isinstance(ci.calculate_frost_days(), xr.DataArray), "Failed to calculate frost days."
        assert isinstance(ci.calculate_warm_spell(), xr.DataArray), "Failed to calculate warm spell."
