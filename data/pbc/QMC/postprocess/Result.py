from dataclasses import dataclass
from Measurement import Measurement
from typing import Dict, Any
from pathlib import Path
import zipfile_deflate64 as z

@dataclass
class Parameters:
    D:int 
    L:int
    N:int 
    M:int 
    l:int 
    U:float
    t:float = 1.0 
    bin_size:int = 10
    subgeometry:str = "square"
    trial_state:str = "constant"
    
    @classmethod
    def from_filename(cls, filename:str):
        import re
        
        pattern = r"(.+?)_(\d+)D_L(\d+)_N(\d+)_M(\d+)_l(\d+)_U([\d.]+)_t([\d.]+)_beta([\d.]+)_binsize(\d+)_seed(\d+)_(\w+)_trialstate_([\w-]+)"

        match = re.match(pattern, filename)

        if match:
            params = { 
                'D': int(match.group(2)),
                'L': int(match.group(3)),
                'N': int(match.group(4)),
                'M': int(match.group(5)),
                'l': int(match.group(6)),
                'U': float(match.group(7)),
                't': float(match.group(8)), 
                'bin_size': int(match.group(10)), 
                'subgeometry': match.group(12),
                'trial_state': match.group(13) 
            }

            print("Extracted Parameters:")
            print(params)
            return cls(**params)
        else:
            raise ValueError("Filename does not match the expected pattern.")
            
            

class Result:
    measurements: Dict[str, Measurement]
    calculation_label: str
    path: str
    
    def __init__(self, parameters:Parameters, path:Path):
        self.measurements = {}
        self.calculation_label = self.get_calulation_label(parameters)
        self.path = path
         
    def get_filename(self, measurement_name:str, beta:float, seed:int):
        return self.measurements[measurement_name].get_filename(seed=seed, beta=beta)
    
    def add_measurement(self, measurement_factory, path:Path=None, additional_kwargs:Dict[str,Any]=None):
        if additional_kwargs is None:
            additional_kwargs = {}
        if path is None:
            path = self.path
        measurement = measurement_factory(calculation_label=self.calculation_label, path=path, **additional_kwargs)
        self.measurements[measurement.name] = measurement
        
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        for measurement in self.measurements.values():
            measurement.add_seed(seed=seed, beta=beta, n_bins_threshold=n_bins_threshold, zip=zip)
    
    def evaluate(self):
        for measurement in self.measurements.values():
            measurement.evaluate()
            
    def get_result(self, measurement_name:str):
        return self.measurements[measurement_name].results
        
    def get_calulation_label(self, parameters:Parameters):
        p = parameters
        return f"{p.D:d}D_L{p.L:d}_N{p.N:d}_M{p.M:d}_l{p.l:d}_U{p.U:.4f}_t{p.t:.4f}_beta{{beta:.4f}}_binsize{p.bin_size:d}_seed{{seed:d}}_{p.subgeometry}_trialstate_{p.trial_state}"

    def clear(self):
        for measurement in self.measurements.values():
            measurement.clear()