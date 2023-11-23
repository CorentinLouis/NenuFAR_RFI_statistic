from scipy.io import readsav

class LazyLoader:
    def __init__(self, data_file, loader_function):
        self.data_file = data_file
        self._data = None  # Placeholder for the loaded data
        self.loader_function = loader_function

    @property
    def data(self):
        if self._data is None:
            # Load the data using the provided loader function when it's first requested
            self._data = self.loader_function(self.data_file)
        return self._data
    
def load_rfistat_data(data_file):
    # Load data from a .sav file
    data = readsav(data_file)
    return data
