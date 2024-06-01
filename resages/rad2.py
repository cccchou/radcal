import os
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp1d
import time
from .abs import abc_rad
from tqdm import tqdm
class  Rad2(abc_rad):
    def __init__(self,name,cc=1,threshold=1e-6,delta_14R=False,
            dcp=False,
            f_14R=False,
            yrsteps=1,
            prob=0.95,
            hpdsteps=1,calibrate=True,export_cal_pdf=False):
        super().__init__(name=name,cc=cc,threshold=threshold,
            delta_14R=delta_14R,
            dcp=dcp,
            f_14R=f_14R,
            yrsteps=yrsteps,
            prob=prob,
            hpdsteps=hpdsteps)
        self.calibrate = calibrate
        self.export_cal_pdf =export_cal_pdf
 
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

        self.calcurve = pd.read_table(self.cc, header=None)
        self.theta = self.calcurve.iloc[:, 0]
        self.f_mu = np.exp(-self.calcurve.iloc[:, 1] / 8033)
        self.f_sigma = np.exp(-(self.calcurve.iloc[:, 1] - self.calcurve.iloc[:, 2]) / 8033) - self.f_mu

        
		
		

    def rad2(self):
        dets = pd.read_csv(f"InputData/{self.name}/{self.name}{self.ext}", sep=self.sep)

        




        self.id_col = dets["id"]
        self.rad_res, self.rad_res_sd, self.rad_atm, self.rad_atm_sd = dets["Reservoir_derived_14C_age"], dets[
            "Reservoir_derived_14C_error"], dets["Atmosphere_derived_14C_age"], dets["Atmosphere_derived_14C_error"]
        self.fm_res, self.fm_res_sd, self.fm_atm, self.fm_atm_sd = np.exp(self.rad_res / -8033), np.exp(-(self.rad_res - self.rad_res_sd) / 8033) - np.exp(
            self.rad_res / -8033), np.exp(self.rad_atm / -8033), np.exp(-(self.rad_atm - self.rad_atm_sd) / 8033) - np.exp(self.rad_atm / -8033)
        self.f_cage = np.exp(-self.rad_atm / 8033)
        self.f_error = np.exp(-(self.rad_atm - self.rad_atm_sd) / 8033) - self.f_cage
        self.dat = pd.DataFrame({"f_cage": self.f_cage, "f_error": self.f_error})
        # calculate sample reservoir ages
        print("Calculating reservoir age offsets...\n")
        self.Reservoir_age_offset = np.round(self.rad_res - self.rad_atm)
        self.Reservoir_age_offset_error = np.round(np.sqrt(self.rad_res_sd ** 2 + self.rad_atm_sd ** 2))

        self.F14R = np.round(self.fm_res / self.fm_atm, 4)
        self.F14R_error = np.round(self.F14R * np.sqrt(self.fm_res_sd ** 2 / (self.fm_res ** 2) + self.fm_atm_sd ** 2 / (self.fm_atm ** 2)), 4)
        self.Delta_14R = 1000 * (self.F14R - 1)
        self.Delta_14R_error = 1000 * self.F14R_error
        self.DCP = 100 * (1 - self.F14R)
        self.DCP_error = 100 * self.F14R_error

        output = pd.concat([dets, pd.DataFrame(
            {"Reservoir_age_offset": self.Reservoir_age_offset, "Reservoir_age_offset_error": self.Reservoir_age_offset_error})],
                           axis=1)
        if self.f_14R:
            output = pd.concat([output, pd.DataFrame({"F14R": self.F14R, "F14R_error": self.F14R_error})], axis=1)
        if self.delta_14R:
            output = pd.concat([output, pd.DataFrame({"delta14R": self.Delta_14R, "delta14R_error": self.Delta_14R_error})],
                               axis=1)
        if self.dcp:
            output = pd.concat([output, pd.DataFrame({"DCP": self.DCP, "DCP_error": self.DCP_error})], axis=1)
        output.to_csv(f"InputData/{self.name}/{self.name}_out{self.ext}", index=False)
        if self.calibrate:
            self._calibrate
        if self.export_cal_pdf:
            self._export_cal_pdf
        print("\nWork has been done successfully, check outputs in your folder.\n")
        return

    def _caldist(self,f_cage,f_error,f_sigma):
        cal = pd.DataFrame({"theta": self.theta, "pdf": norm.pdf(self.f_mu, f_cage, np.sqrt(f_error**2 + f_sigma**2))})
        cal=cal[cal["pdf"]>0]
        f = interp1d(cal.iloc[:, 0], cal.iloc[:, 1])
        newx=np.arange(cal.iloc[:, 0].min(), cal.iloc[:, 0].max()+1, self.yrsteps)
        cal=f(newx)
        cal = pd.DataFrame(data=np.column_stack((newx, cal)), columns=["theta","pdf"])
        cal["pdf"]/= cal["pdf"].sum()
        return cal[cal["pdf"] > self.threshold]

    def _hpd(self,dat):
        f=interp1d(dat["theta"],dat["pdf"])
        newx=np.arange(min(dat["theta"]),max(dat["theta"])+1,self.yrsteps)

        dat = f(newx)

        dat=pd.DataFrame(data=np.column_stack((newx,dat)),columns=["theta","pdf"])
        dat = dat.sort_values(by="pdf", ascending=False)
        dat["pdf"]/=dat["pdf"].sum()

        dat = dat[dat["pdf"].cumsum() <= self.prob]
        dat=dat.sort_values(by="theta",ascending=True)


        dif = np.where(np.diff(dat["theta"]) > self.hpdsteps)[0]

        if len(dif) == 0:
            hpds = pd.DataFrame({"min": [dat["theta"].min()], "max": [dat["theta"].max()], "prob": [100 * self.prob]})
        else:

            dif = np.concatenate((np.array([dat.iloc[0,0]]),np.sort(np.append(dat.iloc[dif,0],dat.iloc[dif+1,0])),dat.iloc[len(dat)-1,0].reshape(-1)))

            dif = dif.reshape(-1, 2)
            probs = [round(100 * dat[(dat["theta"] >= i[0]) & (dat["theta"] <= i[1])]["pdf"].sum(), 1) for i in dif]
            hpds = pd.DataFrame({"min": dif[:, 0], "max": dif[:, 1], "prob": probs})

        return hpds

    @property
    def _calibrate(self):




        print(f"Calibrating atmospheric 14C ages using the {self.cc} calibration curve...\n")

        border = 0
        if any((self.rad_atm - self.rad_atm_sd) < min(self.calcurve.iloc[:, 1] + self.calcurve.iloc[:, 2])):
            if any((self.rad_atm + self.rad_atm_sd) > min(self.calcurve.iloc[:, 1] - self.calcurve.iloc[:, 2])):
                border = 1
            else:
                border = 2
        if any((self.rad_atm + self.rad_atm_sd) > max(self.calcurve.iloc[:, 1] - self.calcurve.iloc[:, 2])):
            if any((self.rad_atm - self.rad_atm_sd) < max(self.calcurve.iloc[:, 1] + self.calcurve.iloc[:, 2])):
                border = 1
            else:
                border = 2
        if border == 1:
            print("\nAt least one 14C age falls partly beyond the calibration curve and will be truncated")
        if border == 2:
            raise ValueError("\nAt least one 14C age cannot be calibrated because it falls beyond the calibration curve! Please remove those.")


        temp=dict()
        print("Loading", end="")
        for i in tqdm(range(len(self.f_cage)),
                      desc="processing",
                      total=len(self.f_cage),
                      colour="blue"):
            calib = self._caldist(f_cage=self.f_cage[i],f_error=self.f_error[i],f_sigma=self.f_sigma[i])
            temp[f"calib_{i}"] = calib["pdf"]
            temp[f"hpd_{i}"] = self._hpd(calib)
            temp[f"mid1_{i}"] = round((temp[f"hpd_{i}"]["min"].values[0] + temp[f"hpd_{i}"]["max"].values[0]) / 2)
            yrs = calib["theta"]
            temp[f"mid2_{i}"] = round(np.mean([max(yrs), min(yrs)]))
            temp[f"wmn_{i}"] = round(np.average(calib["theta"], weights=1/calib["pdf"]))
            temp[f"med_{i}"] =calib[calib["pdf"].cumsum()<=0.5]["theta"].iloc[-1]
            temp[f"mode_{i}"]=calib[calib["pdf"]==calib["pdf"].max()]["theta"].values
        self.temp=temp
        print(f"Calibrating atmospheric 14C ages using the {self.cc} calibration curve has done")
        return




    @property
    def _export_cal_pdf(self):
        print("Exporting your calibrated ages...")
        mini=[]
        maxi=[]
        for i in range(len(self.f_cage)):

            mini.append(min(self.temp[f"calib_{i}"].iloc[:,0]))
            maxi.append(max(self.temp[f"calib_{i}"].iloc[:,0]))
        min_val = min(mini)
        max_val = max(maxi)
        header=[]
        f=interp1d(self.temp[f"calib_{i}"].iloc[:,0],self.temp[f"calib_{i}"].iloc[:,1])
        
        cal_age = np.zeros(
            (int(np.floor(max_val - min_val + 1)), 1 + self.dat.shape[0]))
        cal_age[:, 0] = np.arange(min_val, max_val + 1, self.radsteps)
        for i in range(len(self.dat)):

            d = max(np.diff(pd.DataFrame(self.temp[f"calib_{i}"]).iloc[:,0]))
            if d > 1:
                raise ValueError("Please, decrease the threshold parameter: e.g. threshold = 1e-14 in the command line.")
            cal_age[int(
                np.floor((
                        np.min(self.temp[f'uncalib_{i}'].iloc[:, 0]) -
                        min_val))):int(
                np.floor(
                    np.max(self.temp[f'uncalib_{i}'].iloc[:, 0]) -
                    min_val + 1)),
            i + 1] = self.temp[f'uncalib_{i}'].iloc[:, 1]
            header.append(self.id_col[i])
        cal_age.columns=np.append("cal_year_BP",header)


        cal_age.to_csv(f"InputData/{self.name}/{self.name}{self.cc}_Calibrated_ages_pdfs{self.ext}", index=False)

        for i in range(len(self.f_cage)):
            hpds = self.temp[f"hpd_{i}"]
            test = 1 if len(hpds) == 1 else 2
            if test == 2:
                break

        if test == 1:
            for i in range(len(self.dat)):
                export_calib = pd.DataFrame({"CalibratedMinRange_at_{}%".format(100 * self.prob): hpds["min"].values,
                                                "CalibratedMaxRange_at_{}%".format(100 * self.prob): hpds["max"].values,
                                                "Mode": self.temp[f"mode_{i}"],
                                                "MidRange": self.temp[f"mid1_{i}"],
                                                "Median": self.temp[f"med_{i}"]})

            colnames = [f"{self.id_col[i]}_{col}" for col in export_calib.columns]
            export_calib.columns = colnames

            cal_output = pd.concat([pd.DataFrame(self.id_col),export_calib], axis=1)
            cal_output.to_csv(f"InputData/{self.name}/{self.name}{self.cc}_Calibrated_age_ranges_at_{100 * self.prob}%{self.ext}", index=False)
        else:
            ref = "Calibrated_"
            hpd_file = open(f"InputData/{self.name}/{self.name}{self.cc}_{self.ref}ranges.txt", "w")
            hpd_file.write(f"{ref}ranges at {100 * self.prob}% confidence intervals\n")

            for i in range(len(self.f_cage)):
                hpd_file.write(
                    f"\n\nIdentification_number: {self.id_col[i]}\nCal_year_MinRange\tCal_year_MaxRange\tprobability\n")
                hpds = self.temp[f"hpd_{i}"]
                pdf_info = pd.DataFrame(
                    {"Mode": [self.temp[f"mode_{i}"]], "MidRange": [self.temp[i,f"mid1_{i}"]], "Median": [self.temp[f"med_{i}"]]})

                for j in range(len(hpds)):
                    for k in range(3):
                        hpd_file.write(f"{hpds.iloc[j, k]}\t")
                    hpd_file.write("\n")

                hpd_file.write("Mode\tMidRange\tMedian\n")
                for l in range(len(pdf_info.columns)):
                    hpd_file.write(f"{pdf_info.iloc[0, l]}\t")

                hpd_file.write("\n")

            hpd_file.close()
        print("Exporting your calibrated ages has done")

        

