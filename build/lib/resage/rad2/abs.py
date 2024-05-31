from abc import ABCMeta, abstractmethod
import pandas as pd
class abc_rad(metaclass=ABCMeta):
    def __init__(self,
                 name="example_radcal_Southon2012",
                 cc=1,
                 yrsteps=1,
                 prob=0.95,
                 threshold=1e-6,
                 hpdsteps=1,
                 delta_14R=False,
                 dcp=False,
                 f_14R=False,
                 ext =".csv",
                 sep = ","
                 ):
        self.name = name
        self.cc = cc
        self.yrsteps = yrsteps
        self.prob = prob
        self.threshold = threshold
        self.hpdsteps = hpdsteps
        self.delta_14R = delta_14R
        self.dcp = dcp
        self.f_14R = f_14R
        self.ext = ext
        self.sep = sep
        
        
    
    @abstractmethod
    def _hpd(self):
        pass

    @property
    def _assert(self):
        dets = pd.read_csv(f"InputData/{self.name}/{self.name}{self.ext}",
                           sep=self.sep)

        # check for consistency
        if dets.shape[1] != 5:
            raise ValueError(
                "Your input file should contain 5 columns.\nPlease, check the template designed for the radcal function."
            )

        if min(dets.iloc[:, 1]) < 0 or min(dets.iloc[:, 3]) < 0:
            raise ValueError(
                "At least one of your radiocarbon data or one of your calendar age is negative."
                "\nResAge is not designed to work with post-bomb samples."
                "\nRemove such data from the input file.")

        if self.delta_14R and (self.dcp or self.f_14R or self.dcp or self.f_14R
                               or self.dcp and self.f_14R and self.delta_14R):
            raise ValueError(
                "Please, choose only one metric: reservoir age offset (default), f.14R, delta.14R, or dcp"
            )