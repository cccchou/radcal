import os
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp1d
import time
from .abs import abc_rad
from tqdm import tqdm

class colyear(abc_rad):
    def __init__(self,name,cc=1,delta_14R=False,
            dcp=False,
            f_14R=False,
            yrsteps=1,
            ):
        super().__init__(name=name,cc=cc,
            delta_14R=delta_14R,
            dcp=dcp,
            f_14R=f_14R,
            yrsteps=yrsteps,
            )
        match self.cc:
            case 1:
                self.cc = "IntCal20.14C"
            case 2:
                self.cc = "SHCal13.14C"
            case 3:
                self.cc = "IntCal13.14C"
            case _:
                raise ValueError(
                    "Invalid value for cc. Please check the manual.")
        self._assert
        self.dets = pd.read_csv(f"InputData/{self.name}/{self.name}{self.ext}", sep=self.sep)
  
        self.calcurve = pd.read_table(self.cc, header=None)
        
    def _assert(self):
        if self.dets.shape[1] != 4:
            raise ValueError("Your input file should contain 4 columns. Please, check the template designed for the colyear function.")
        
        if np.min(self.dets.iloc[:,1]) < 0:
            raise ValueError("At least one of your radiocarbon data is negative.\n Colyear is not designed to work with post-bomb samples. Remove such data from the input file.")
        
        if(max(self.dets.iloc[:,3]) > 1950) or (min(self.dets.iloc[:,3]) < 1500):
            raise ValueError("Colyear function is only designed for sample with calendar age between 1950 and 1500 AD. Remove such data from the intput file")
        
        if self.dets.iloc[:,3] <= 450:
            raise ValueError("Did you already convert sample collection year from AD/BC timescale to the BP timescale?\n Please let collection year data in the AD/BC timescale.")
            
        if self.yrsteps != 1:
            raise ValueError("The collection year (colyear) function requires the yrsteps argument is set to 1.\n In calling function, set yrsteps=1.")
        
        
    def colyear(self):
        print(f"Reservoir age offsets are calculated using the {self.cc} calibration curve.\n")
  
        # interpolate calibration curve between 0 and 500 BP, i.e. between 1950 and 1500 AD
        x = np.arange(0, 451, self.yrsteps)
        y = interp1d(self.calcurve.iloc[:,0], self.calcurve.iloc[:,1])(x)
        sd = interp1d(self.calcurve.iloc[:,0], self.calcurve.iloc[:,2])(x)
        interp = np.column_stack((x, y, sd))
  
        # now work in the BP timescale
        colyr_AD = self.dets["Collection_year_AD"]
        Collection_year_BP = 1950 - colyr_AD
  
        # prepare for reservoir age calculation 
        rad_res = self.dets["Radiocarbon_age"]
        rad_res_sd = self.dets["Radiocarbon_age_error"]
        
        # calculate sample reservoir ages
        Atmospheric_14C_age = np.round(interp[Collection_year_BP,1]) # find the 14C age of the atmosphere during collection year... 
        Atmospheric_14C_error = np.round(interp[Collection_year_BP,2]) # ... and associated uncertainty
        Reservoir_age = round(rad_res - Atmospheric_14C_age) # calculate reservoir age
        Reservoir_age_error = round(np.sqrt(rad_res_sd*rad_res_sd + Atmospheric_14C_error*Atmospheric_14C_error)) # propagate errors through quadratic sum
        
        fm_res = np.exp(rad_res/-8033)
        fm_res_sd =np.exp(-(rad_res-rad_res_sd)/8033) - fm_res
        fm_atm = np.exp(Atmospheric_14C_age/-8033)
        fm_atm_sd = np.exp(-(Atmospheric_14C_age-Atmospheric_14C_error)/8033) - fm_atm
        
        F14R  = round(fm_res/fm_atm,4) # calculate F14R rounded to the fourth digit
        F14R_error = round(F14R*np.sqrt(fm_res_sd*fm_res_sd/(fm_res*fm_res)+fm_atm_sd*fm_atm_sd/(fm_atm*fm_atm)),4) # propagates errors
        delta14R = 1000*(F14R-1) # calculate delta14R
        delta14R_error = 1000*F14R_error # propagate errors
        DCP = 100*(1-F14R) # calculate dcp
        DCP_error = 100*F14R_error
        
        output = np.column_stack((self.dets, Collection_year_BP, Atmospheric_14C_age, Atmospheric_14C_error, Reservoir_age, Reservoir_age_error)) 
        header=[id for id in self.dets.columns]
        header.extend(["Collection_year_BP", "Atmospheric_14C_age", "Atmospheric_14C_error", "Reservoir_age", "Reservoir_age_error"])
        if self.f_14R:
            output = np.column_stack((output, F14R, F14R_error))
            header.extend(["F14R","F14R_error"])
        if self.delta_14R: 
            output = np.column_stack((output, delta14R, delta14R_error))
            header.extend(["delta14R","delta14R_error"])
        if self.dcp: 
            output = np.column_stack((output, DCP, DCP_error))
            header.extend(["DCP","DCP_error"])
        
        
        
        output=pd.DataFrame(output,columns=header)
        output.to_csv(f"InputData/{self.name}/{self.name}_out{self.ext}",index=False)
        # np.savetxt( f"InputData/{self.name}/{self.name}_out{self.ext}", output, delimiter="," )
        
        # Done message
        print("Work has been done successfully, check output in your folder.\n\n")
    #Since the abstract function requires _hpd it needs to be written as a blank function
    def _hpd(self):
        return
        
        



