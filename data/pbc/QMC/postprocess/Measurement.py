from typing import Protocol, Dict 
from pathlib import Path
import numpy as np
import zipfile_deflate64 as z

class Measurement(Protocol):
    _calculation_label: str 
    zip:z.ZipFile
    _data: Dict[float, np.ndarray]
    results: Dict[float, np.ndarray]
    name: str
    
    def __init__(self, calculation_label:str, path:Path):
        """Add filename template (caluclation_label) and data folder (path) to load data for different seeds and beta from."""
    def get_filename(self, seed:int, beta:float):
        """Return the full filepath for a given seed and beta."""
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        """Add results for a seed and beta to the data dictionary (key:beta, value:concatenated array for all seeds)."""
    def evaluate(self):
        """Average the seeds, finite size scale beta, and calculate the error."""
        
class KineticEnergy:
    def __init__(self, calculation_label:str, path:Path):
        self._calculation_label = calculation_label
        self.path = path
        self._data = {}
        self._means = {}
        self.name = "K" 
        self.results = {}
    def get_filename(self, seed:int, beta:float):
        calculation_label = self._calculation_label.format(seed=seed, beta=beta)
        fn = f"K_{calculation_label}.dat"
        return fn
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        if zip is None:
            zip = z.ZipFile(self.path(beta=beta,seed=seed))
        data = np.loadtxt(str(zip.read(self.get_filename(seed, beta))).split("\\n")[1:-1]) 
        if len(data)>n_bins_threshold:
            self._data[beta] = np.concatenate((self._data[beta], data)) if beta in self._data else data 
            if beta in self._means:
                self._means[beta].append(np.mean(data))  
            else:
                self._means[beta] = [np.mean(data)]
        else: 
            print(f"Skipping seed {seed} for beta {beta} because it has only {len(data)} bins.")
    def clear(self):
        self._data = {}
        self._means = {}
        self.results = {}
    def evaluate(self):
        for beta in self._data:
            self.results[beta] = np.mean(self._data[beta]), np.std(self._means[beta])/np.sqrt(len(self._means[beta]))
    
class PotentialEnergy:
    def __init__(self, calculation_label:str, path:Path):
        self._calculation_label = calculation_label
        self.path = path
        self._data = {}
        self._means = {}
        self.name = "V" 
        self.results = {}
    def get_filename(self, seed:int, beta:float):
        calculation_label = self._calculation_label.format(seed=seed, beta=beta)
        fn = f"V_{calculation_label}.dat"
        return fn
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        if zip is None:
            zip = self.zip
        data = np.loadtxt(str(zip.read(self.get_filename(seed, beta))).split("\\n")[1:-1])  
        if len(data)>n_bins_threshold:
            self._data[beta] = np.concatenate((self._data[beta], data)) if beta in self._data else data 
            if beta in self._means:
                self._means[beta].append(np.mean(data))  
            else: 
                self._means[beta] = [np.mean(data)] 
        else: 
            print(f"Skipping seed {seed} for beta {beta} because it has only {len(data)} bins.")
    def clear(self):
        self._data = {}
        self._means = {}
        self.results = {}
    def evaluate(self):
        for beta in self._data:
            self.results[beta] = np.mean(self._data[beta]), np.std(self._means[beta])/np.sqrt(len(self._means[beta]))

class Energy:
    def __init__(self, calculation_label:str, path:Path):
        self._calculation_label = calculation_label
        self.path = path
        self._data = {}
        self._means = {}
        self.name = "E" 
        self.results = {}
    def get_filename(self, seed:int, beta:float):
        calculation_label = self._calculation_label.format(seed=seed, beta=beta)
        fn_V = f"V_{calculation_label}.dat"
        fn_K = f"K_{calculation_label}.dat"
        return fn_K, fn_V
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        if zip is None:
            zip = self.zip
        fn_K, fn_V = self.get_filename(seed, beta)
        data_K = np.loadtxt(str(zip.read(fn_K)).split("\\n")[1:-1])  
        data_V = np.loadtxt(str(zip.read(fn_V)).split("\\n")[1:-1])  
        if len(data_K)>n_bins_threshold:
            Ks = np.concatenate((self._data[beta][0], data_K)) if beta in self._data else data_K 
            Vs = np.concatenate((self._data[beta][1], data_V)) if beta in self._data else data_V
            self._data[beta] = (Ks, Vs)
            if beta in self._means:
                self._means[beta][0].append(np.mean(data_K))
                self._means[beta][1].append(np.mean(data_V)) 
            else: 
                self._means[beta] = ([np.mean(data_K)], [np.mean(data_V)]) 
        else: 
            print(f"Skipping seed {seed} for beta {beta} because it has only {len(data_K)} bins.")
    def clear(self):
        self._data = {}
        self._means = {}
        self.results = {}
    def evaluate(self):
        for beta in self._data:
            # means
            ex_K = np.mean(self._data[beta][0])
            ex_V = np.mean(self._data[beta][1])
            # errors
            dK = np.std(self._means[beta][0])/np.sqrt(np.size(self._means[beta][0]))
            dV = np.std(self._means[beta][1])/np.sqrt(np.size(self._means[beta][1]))
            # result = (mean, error)
            self.results[beta] =  (ex_K + ex_V, dK+dV)
    
class ParticleNumberFluctuationsFromCorrelationsPBC:
    def __init__(self, calculation_label:str, path:Path):
        self._calculation_label = calculation_label
        self.path = path
        self._data = {}
        self._means = {}
        self.name = "F_corr_pbc" 
        self.results = {}
    def get_filename(self, seed:int, beta:float):
        calculation_label = self._calculation_label.format(seed=seed, beta=beta) 
        fn_F = f"sigma2_{calculation_label}.dat"
        return fn_F
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        if zip is None:
            zip = self.zip
        fn_F = self.get_filename(seed, beta) 
        data_F = np.loadtxt(str(zip.read(fn_F)).split("\\n")[1:-1])   
        if len(data_F)>n_bins_threshold:
            if len(data_F):
                F = np.concatenate((self._data[beta], data_F)) if beta in self._data else data_F 
                self._data[beta] = F 
                if beta in self._means:
                    self._means[beta].append(np.mean(data_F,axis=0)) 
                else: 
                    self._means[beta] = [np.mean(data_F,axis=0)] 
        else: 
            print(f"Skipping seed {seed} for beta {beta} because it has only {len(data_F)} bins.")
    def clear(self):
        self._data = {}
        self._means = {}
        self.results = {}
    def evaluate(self):
        for beta in self._data:
            # means
            F = np.mean(self._data[beta],axis=0)
            # errors 
            dF = np.std(self._means[beta],axis=0)/np.sqrt(np.size(self._means[beta],axis=0)) 
            # naive
            self.results[beta] = (F, dF)
        
     
class ParticleNumberFluctuations:
    def __init__(self, calculation_label:str, path:Path):
        self._calculation_label = calculation_label
        self.path = path
        self._data = {}
        self._means = {}
        self.name = "F" 
        self.results = {}
    def get_filename(self, seed:int, beta:float):
        calculation_label = self._calculation_label.format(seed=seed, beta=beta)
        fn_n = f"n_{calculation_label}.dat"
        fn_n2 = f"n_squared_{calculation_label}.dat"
        return fn_n, fn_n2
    def add_seed(self, seed:int, beta:float, n_bins_threshold:int=0, zip:z.ZipFile=None):
        if zip is None:
            zip = self.zip
        fn_n, fn_n2 = self.get_filename(seed, beta) 
        data_n = np.loadtxt(str(zip.read(fn_n)).split("\\n")[1:-1]) 
        data_n2 = np.loadtxt(str(zip.read(fn_n2)).split("\\n")[1:-1]) 
        if len(data_n)>n_bins_threshold:
            if len(data_n) and len(data_n2):
                ns = np.concatenate((self._data[beta][0], data_n)) if beta in self._data else data_n
                n2s = np.concatenate((self._data[beta][1], data_n2)) if beta in self._data else data_n2
                self._data[beta] = (ns, n2s)
                if beta in self._means:
                    self._means[beta][0].append(np.mean(data_n,axis=0))
                    self._means[beta][1].append(np.mean(data_n2,axis=0)) 
                else: 
                    self._means[beta] = ([np.mean(data_n,axis=0)], [np.mean(data_n2,axis=0)]) 
        else: 
            print(f"Skipping seed {seed} for beta {beta} because it has only {len(data_n)} bins.")
    def clear(self):
        self._data = {}
        self._means = {}
        self.results = {}
    def evaluate(self):
        for beta in self._data:
            # means
            ex_n = np.mean(self._data[beta][0],axis=0)
            ex_n2 = np.mean(self._data[beta][1],axis=0)
            # errors 
            dn = np.std(self._means[beta][0],axis=0)/np.sqrt(np.size(self._means[beta][0],axis=0))
            dn2 = np.std(self._means[beta][1],axis=0)/np.sqrt(np.size(self._means[beta][1],axis=0))
            # naive
            self.results[beta] = (ex_n2 - ex_n**2, dn2 + np.abs(2*ex_n)*dn)
    