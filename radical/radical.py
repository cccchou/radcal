import numpy as np
from scipy.stats import norm
import pandas as pd
from scipy.interpolate import interp1d
import os
import time
from tqdm import tqdm
from colorama import Fore


class Radical(object):

    def __init__(self,
                 name="example_radcal_Southon2012",
                 cc=1,
                 yrsteps=1,
                 uncalsteps=5,
                 radsteps=1,
                 convolsteps=1,
                 prob=0.95,
                 threshold=1e-6,
                 hpdsteps=1,
                 storedat=False,
                 export_uncal_pdf=False,
                 mixture_pdf=False,
                 delta_14R=False,
                 dcp=False,
                 f_14R=False,
                 export_resage_pdf=True):
        """
        :param name: file name
        :param cc: calibrate name
        :param cc1: IntCal20.14C
        :param cc2: SHCal13.14C
        :param cc3: IntCal13.14C
        :param ext: file format
        :param sep: read file separate way
        :param yrsteps: resolution by calibrate age
        :param uncalsteps: resolution by uncalibrate
        :param radsteps:resolution by radiocarbon
        :param convolsteps:resolution by convolution
        :param prob: the maximum possibility
        :param threshold: threshold
        :param hpdsteps: resolution by hpd
        :param storedat: whether store data
        :param export_uncal_pdf: whether export pdf
        :param mixture_pdf: whether mixture pdf
        :param delta_14R: whether export 14C
        :param dcp:whether export dcp
        :param f_14R:whether export f_14C

        Soulet G, 2015. Methods and codes for reservoir-atmosphere 14C age offset calculations.Quaternary Geochronology 29:97-103, doi: 10.1016/j.quageo.2015.05.023
        """
        self.name = name
        match cc:
            case 1:
                self.cc = "IntCal20.14C"
            case 2:
                self.cc = "SHCal13.14C"
            case 3:
                self.cc = "IntCal13.14C"
            case _:
                raise ValueError(
                    "Invalid value for cc. Please check the manual.")
        self.calcurve = pd.read_table(self.cc, header=None)
        self.ext = ".csv"
        self.sep = ","
        self.yrsteps = yrsteps
        self.uncalsteps = uncalsteps
        self.radsteps = radsteps
        self.convolsteps = convolsteps
        self.prob = prob
        self.threshold = threshold
        self.hpdsteps = hpdsteps
        self.storedat = storedat
        self.export_uncal_pdf = export_uncal_pdf
        self.export_resage_pdf = export_resage_pdf
        self.mixture_pdf = mixture_pdf
        self.delta_14R = delta_14R
        self.dcp = dcp
        self.f_14R = f_14R

        self.theta = self.calcurve.iloc[:, 0]
        self.f_mu = np.exp(-self.calcurve.iloc[:, 1] / 8033)
        self.f_sigma = np.exp(
            -(self.calcurve.iloc[:, 1] - self.calcurve.iloc[:, 2]) /
            8033) - self.f_mu
        self._assert

    def radcal(self):
        #1.read the data
        #2.canculate pdf(uncalidist/res)
        #3.canculate hpd

        #updated data by interpolate the step year

        calcurve = self.calcurve
        self.theta = np.arange(0, self.theta.max() + 1, self.yrsteps)
        self.rad = np.round(
            interp1d(calcurve.iloc[:, 0], calcurve.iloc[:, 1])(self.theta))
        self.rad_sigma = np.round(
            interp1d(calcurve.iloc[:, 0], calcurve.iloc[:, 2])(self.theta))
        self.f_mu = interp1d(calcurve.iloc[:, 0], self.f_mu,
                             kind='linear')(self.theta)
        self.f_sigma = interp1d(calcurve.iloc[:, 0],
                                self.f_sigma,
                                kind='linear')(self.theta)

        dets = pd.read_csv(f"InputData/{self.name}/{self.name}{self.ext}",
                           sep=self.sep)
        self.dets = dets
        id_col = dets["id"]
        self.id_col = id_col
        rad_res = dets["Reservoir_derived_14C_age"]
        self.rad_res = rad_res
        rad_res_sd = dets["Reservoir_derived_14C_error"]
        self.rad_res_sd = rad_res_sd
        cal_age = dets["Calendar_age"]
        self.cal_age = cal_age
        cal_sigma = dets["Calendar_age_error"]
        self.cal_sigma = cal_sigma

        self.dat = pd.DataFrame({
            'id': self.id_col,
            'cal_age': cal_age,
            'cal_sigma': cal_sigma,
            'rad_res': rad_res,
            'rad_res_sd': rad_res_sd
        })

        temp = dict()
        print("\n.running")
        now = time.time()
        for i in tqdm(range(len(self.dat['cal_age'])),
                      desc="processing",
                      total=len(self.dat['cal_age']),
                      colour="blue"):
            uncalib = self._uncaldist(self.dat['cal_age'][i],
                                      self.dat['cal_sigma'][i])  #prior
            temp[
                f"uncalib_{i}"] = uncalib  # resulting uncalibrated (atm radiocarbon distribution) distribution
            temp[f"hpd_uncal_{i}"] = self._hpd(
                data=uncalib
            )  # uncalibrated highest probability density corresponding to prob posterior
            
            self.dat.loc[i, "mid1"] = round(
                (temp[f"hpd_uncal_{i}"]["min"].values[0] +
                 temp[f"hpd_uncal_{i}"]["max"].values[0]) /
                2)  # midpoints of uncalibrated ranges
            radyrs = uncalib["theta"]
            self.dat.loc[i, "mid2"] = round(
                np.mean([max(radyrs), min(radyrs)])
            )  # midpoints of entire uncalibrated distributions (with probabilities beyond threshold)
            self.dat.loc[i, "wmn"] = round(
                np.average(
                    uncalib["theta"], weights=1 /
                    uncalib["pdf"]))  # weighted means of uncalibrated ranges
            self.dat.loc[i, "med"] = uncalib[uncalib["pdf"].cumsum() <= 0.5][
                "theta"].iloc[-1]  # medians of uncalibrated distributions
            self.dat.loc[i, "mode"] = uncalib[
                uncalib["pdf"] == uncalib["pdf"].max(
                )]["theta"].values  # maximum densities of uncalibrated distributions

            # Implement or replace convolution function
        
            resage = self._convolution(uncalib, self.dat['rad_res'][i],
                                       self.dat['rad_res_sd'][i])
            temp[f'resage_{i}'] = resage  # resulting reservoir age

            temp[f"hpd_resage_{i}"] = self._hpd(
                data=resage
            )  # reservoir age highest probability density corresponding to prob

            self.dat.loc[i, 'mid1_resage'] = round(
                (temp[f"hpd_resage_{i}"]["min"].values[0] +
                 temp[f"hpd_resage_{i}"]["max"].values[0]) /
                2)  # midpoints of reservoir ranges
            resyrs = resage.iloc[:, 0]
            self.dat.loc[i, 'mid2_resage'] = round(
                np.mean([max(resyrs), min(resyrs)])
            )  # midpoints of entire reservoir age distributions (with probabilities beyond threshold)
            self.dat.loc[i, 'wmn_resage'] = round(
                np.average(
                    resage.iloc[:, 0], weights=1 /
                    resage.iloc[:, 1]))  # weighted means of reservoir ranges
            self.dat.loc[i, 'med_resage'] = resage[resage["pdf"].cumsum(
            ) <= 0.5]["theta"].iloc[-1]  # medians of reservoir distributions
            self.dat.loc[
                i, 'mode_resage'] = resage[resage["pdf"] == resage["pdf"].max(
                )]["theta"].values  # maximum densities of reservoir distributions

            

        print(
            f"the total spending time is {(time.time()-now)//60}min{round(time.time()-now-60*((time.time()-now)//60),1)}s "
        )
        self.temp = temp

        if self.export_uncal_pdf:

            self._export_uncal_pdf
        if self.export_resage_pdf:

            self._export_resage_pdf

        if self.mixture_pdf:

            self._mixture_pdf

        self.dat.to_csv(
            f"InputData/{self.name}/{self.name}_using_{self.cc}_output.csv",
            index=False)

        print(
            "Work has been done successfully, check outputs in your folder.\n")

    def _uncaldist(self, cal_age, cal_sigma):
        """
        
        :param cal_age: the calculate age, theta for N
        :param cal_sigma: the calculate sigma for N
        :param theta: theta值；t for N
        :param f_mu: 平均值
        :param f_sigma: 方差
        :param yrsteps: 年分辨率
        :param threshold: 阈值
        :param uncalsteps: 未矫正的年龄分辨率
        :param radsteps: 储库年分辨率
        :return: dataframe{"theta","pdf"}

        This function calculates the uncalibrated distribution of radiocarbon age offsets
        given a set of calibrated age and uncertainty.

        P(Y) and P(t) == const
        P(Y|t) = N(t, theta, sigma)
        P(r|t) = N(r, p(t), sigma(t))
        P(Y|t,r) & P(Y|t) * p(r|t) * p(t)
        P(Y|r) & ∫ P(Y|t) * P(r|t) * dt

        To reduce calculation time, we only work over 6 time sigma of the sample calendar age
        """
        # 1) work over 6 time sigma of the sample calendar age
        # Reduce calculation time
        #cal_age is observed calendar age
        cal_span = self.theta[(self.theta >= cal_age - 6 * cal_sigma)
                              & (self.theta <= cal_age + 6 * cal_sigma)]
        # 2) select to suitable range of the radiocarbon time scale to calibrate and convert it into F14C
        rad_span = np.round(
            np.arange(
                min(self.rad[cal_span]) - 50 * max(self.rad_sigma[cal_span]),
                max(self.rad[cal_span]) + 50 * max(self.rad_sigma[cal_span]) +
                1, self.uncalsteps))
        
        f_span = np.exp(-rad_span / 8033)

        # Find how far each calendar time of cal.span is from cal.age
        LKH_cal_age = np.column_stack(
            (cal_span, norm.pdf(cal_span, cal_age, cal_sigma)))

        # Uncalibrate
        uncal = []
        
        for i in range(len(f_span)):
            # How far is a single F14C from the calibration curve f.mu
            LKH_cc = norm.pdf(self.f_mu[cal_span], f_span[i],
                              self.f_sigma[cal_span])

            LKH_ri = self.yrsteps * np.sum(LKH_cal_age[:, 1] * LKH_cc)

            uncal.append(LKH_ri)

        uncal = np.column_stack((rad_span, uncal))

        # Interpolate and normalise uncalibrated distribution to 1
        uncal = uncal[uncal[:, 1] > 0, :]  # remove unnecessary data

        f = interp1d(uncal[:, 0], uncal[:, 1])

        uncal = np.column_stack(
            (np.arange(min(uncal[:, 0]),
                       max(uncal[:, 0]) + 1, self.radsteps),
             f(np.arange(min(uncal[:, 0]),
                         max(uncal[:, 0]) + 1, self.radsteps)) /
             np.sum(
                 f(
                     np.arange(min(uncal[:, 0]),
                               max(uncal[:, 0]) + 1, self.radsteps)))))

        # Only report those normalised uncalibrated probabilities beyond a threshold
        uncal = pd.DataFrame(data=uncal, columns=["theta", "pdf"])
        return uncal[uncal["pdf"] > self.threshold]


    def _convolution(self, uncal_dist: pd.DataFrame, rad_res: float,
                     rad_res_sd: float) -> pd.DataFrame:
        """
        Convolution between uncalibrated distribution of ages and
        reservoir 14C ages to obtain the posterior distribution of ages.

        :param uncal_dist: Uncalibrated distribution of ages
        :param rad_res: Mean of reservoir 14C ages
        :param rad_res_sd: Standard deviation of reservoir 14C ages
        :param threshold: Minimum probability to be used in the convolution
        :param convolsteps: Stepsize for the convolution
        :param radsteps: Stepsize for the creation of the spl_pdf array
        :return: Posterior distribution of ages

        The convolution is done according to the formula:

            p(d14R|Y) = p(r|Y)res*(-1*Ratm*p(r|Y)atm

        where Y is the uncalibrated age.
        """
        spl_span = np.arange(rad_res - 10 * rad_res_sd,
                             rad_res + 10 * rad_res_sd + 1, self.radsteps)
        spl_pdf = np.column_stack(
            (spl_span, norm.pdf(spl_span, rad_res, rad_res_sd)))
        spl_pdf[:, 1] /= np.sum(spl_pdf[:, 1])
        spl_pdf = spl_pdf[spl_pdf[:, 1] > self.threshold, :]

        # Move the end of matrix dat (uncalibrated age) to the start of
        # reservoir 14C age spl.pdf
        new_origin = np.min(spl_pdf[:, 0]) - np.max(uncal_dist.iloc[:, 0])
        uncal_conv = np.column_stack(
            (uncal_dist.iloc[:, 0] + new_origin - 1, uncal_dist.iloc[:, 1]))

        # Make matrices same dimension so that they overlap
        m, n = spl_pdf.shape[0], uncal_conv.shape[0]
        s = np.column_stack((uncal_conv[:, 0], np.zeros(n)))

        # Print uncal_conv[:, 1]
        spl_pdf_new = np.row_stack((s, spl_pdf))
        o = np.column_stack((spl_pdf[:, 0], np.zeros(m)))
        uncal_conv_new = np.row_stack((uncal_conv, o))
        spl_density, uncal_density = spl_pdf_new[:, 1], uncal_conv_new[:, 1]

        # Convolution !
        LKH_res = []

        for i in range(0, n + m + n, self.convolsteps):

            spl_density = np.append(spl_density, np.zeros(self.convolsteps))

            uncal_density = np.append(np.zeros(self.convolsteps),
                                      uncal_density)

            LKH_res.append(np.sum(spl_density * uncal_density))

        res_age = np.column_stack(
            (new_origin - 1 + np.arange(1, n + m + n + 1, self.convolsteps),
             LKH_res))

        res_age = res_age[res_age[:, 1] > 0, :]
        res_age = pd.DataFrame(data=res_age, columns=["theta", "pdf"])

        return res_age[res_age["pdf"] > self.threshold]

    def _hpd(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate the Highest Posterior Density (HPD) interval for a given
        uncalibrated distribution of ages.

        The HPD interval is the smallest interval that contains a specified
        probability of the posterior distribution. The interval is defined by
        the minimum and maximum values of the distribution, and the probability
        of the interval being contained in the distribution.

        Parameters
        ----------
        data : pd.DataFrame
            Uncalibrated distribution of ages. The dataframe should have two
            columns: "theta" and "pdf". The "theta" column should contain the
            age values, and the "pdf" column should contain the corresponding
            probability density values.

        Returns
        ------- pd.DataFrame
            A dataframe with three columns: "min", "max", and "prob". The
            "min" and "max" columns contain the lower and upper bounds of the
            HPD interval, respectively. The "prob" column contains the
            probability of the HPD interval being contained in the
            distribution.
        """
        # sort the data by age
        dat = pd.DataFrame(data=data, columns=["theta", "pdf"])
        f = interp1d(dat["theta"], dat["pdf"])
        newx = np.arange(min(dat["theta"]),
                         max(dat["theta"]) + 1, self.yrsteps)

        dat = f(newx)

        dat = pd.DataFrame(data=np.column_stack((newx, dat)),
                           columns=["theta", "pdf"])

        # normalize the pdf column to get the cumulative probability
        dat["pdf"] /= dat["pdf"].sum()

        # select the data points that are below the specified probability
        # threshold
        dat = dat[dat["pdf"].cumsum() <= self.prob]
        # sort the data points by age again
        dat = dat.sort_values(by="theta", ascending=True)

        # calculate the HPD interval
        if len(dat) == 1:
            # if there is only one data point, the HPD interval is the entire
            # age range
            hpds = pd.DataFrame({
                "min": [dat["theta"].min()],
                "max": [dat["theta"].max()],
                "prob": [100 * self.prob]
            })
        else:
            # find the indices of the data points where the age increases
            # by more than the specified step size
            dif = np.where(np.diff(dat["theta"]) > self.hpdsteps)[0]

            # calculate the HPD interval from the indices of the data points
            if len(dif) == 0:
                # if the age increases monotonically, the HPD interval is
                # the entire age range
                hpds = pd.DataFrame({
                    "min": [dat["theta"].min()],
                    "max": [dat["theta"].max()],
                    "prob": [100 * self.prob]
                })
            else:
                # construct the HPD interval from the indices of the data points
                dif = np.concatenate(
                    (np.array([dat.iloc[0, 0]]),
                    np.sort(np.append(dat.iloc[dif, 0],
                                      dat.iloc[dif + 1,
                                              0])), dat.iloc[len(dat) - 1,
                                                              0].reshape(-1, 1))
                dif = dif.reshape(-1, 2)
                # calculate the probability of each HPD interval
                probs = [
                    round(
                        100 * dat[(dat["theta"] >= i[0])
                                  & (dat["theta"] <= i[1])]["pdf"].sum(), 1)
                    for i in dif
                ]
                hpds = pd.DataFrame({
                    "min": dif[:, 0],
                    "max": dif[:, 1],
                    "prob": probs
                })

        return hpds


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

        # prepare for reservoir age calculation

    @property
    def _export_uncal_pdf(self):
        """Export the uncalibrated ages.

        Returns:
            A pandas DataFrame containing the minimum, maximum, and median (mode and midrange) uncalibrated ages for each sample.
            The output is saved as a CSV file with the name "InputData/{name}/{name}_{cc}_Uncalibrated_age_ranges_at_{100 * prob}%{ext}".
            Alternatively, the output can be saved as a text file "InputData/{name}/{name}_{cc}_Uncalibrated_ranges.txt" with the confidence intervals for each sample.

        """
        print("Exporting your uncalibrated ages")
        mini = []
        maxi = []
        for i in range(len(self.dat)):
            mini.append(min(self.temp[f"uncalib_{i}"]["theta"]))
            maxi.append(max(self.temp[f"uncalib_{i}"]["theta"]))
        min_val = min(mini)
        max_val = max(maxi)
        header = []
        cal_age = np.zeros(
            (int(np.floor(max_val - min_val + 1)), 1 + self.dat.shape[0]))
        cal_age[:, 0] = np.arange(min_val, max_val + 1, self.radsteps)
        # cal_age = pd.DataFrame({"uncal_year_BP": np.arange(min_val, max_val + 1, self.yrsteps)})
        for i in range(len(self.dat)):
            d = max(np.diff(self.temp[f"uncalib_{i}"]["theta"]))
            if d > 1:
                raise ValueError(
                    "Please, decrease the threshold parameter: e.g. threshold = 1e-14 in the command line."
                )
            cal_age[int(
                np.floor((
                    np.min(self.temp[f'uncalib_{i}'].iloc[:, 0]) -
                    min_val))):int(
                        np.floor(
                            np.max(self.temp[f'uncalib_{i}'].iloc[:, 0]) -
                            min_val + 1)),
                    i + 1] = self.temp[f'uncalib_{i}'].iloc[:, 1]

            header.append(self.id_col[i])


        for i in range(len(self.dat)):
            hpds = self.temp[f"hpd_uncal_{i}"]
            test = 1 if len(hpds) == 1 else 2
            if test == 2:
                break

        if test == 1:
            # Export the uncalibrated age ranges as a single CSV file
            export_calib = pd.DataFrame({
                "UncalibratedMinRange_at_{}%".format(100 * self.prob):
                    hpds["min"].values,
                "UncalibratedMaxRange_at_{}%".format(100 * self.prob):
                    hpds["max"].values,
                "Mode":
                    self.dat[f"mode"][i],
                "MidRange":

                    self.dat[f"mid1"][i],
                "Median":

                    self.dat[f"med"][i]
            })

            colnames = [
                f"{self.id_col[i]}_{col}" for col in export_calib.columns
            ]
            export_calib.columns = colnames

            cal_output = pd.concat([pd.DataFrame(self.id_col), export_calib], axis=1)
            cal_output.to_csv(
                f"InputData/{self.name}/{self.name}{self.cc}_Uncalibrated_age_ranges_at_{100 * self.prob}%{self.ext}",
                index=False
            )
        else:
            # Export the uncalibrated age ranges as a text file
            ref = "Uncalibrated_"
            hpd_file = open(
                f"InputData/{self.name}/{self.name}{self.cc}_{ref}ranges.txt",
                "w"
            )
            hpd_file.write(
                f"{ref}ranges at {100 * self.prob}% confidence intervals\n"
            )

            for i in range(len(self.dat)):
                hpd_file.write(
                    f"\n\nIdentification_number: {self.id_col[i]}\nCal_year_MinRange\tCal_year_MaxRange\tprobability\n"
                )
                hpds = self.temp[f"hpd_uncal_{i}"]
                pdf_info = pd.DataFrame({
                    "Mode": [self.dat.loc[i, f"mode"]],
                    "MidRange": [self.dat.loc[i, f"mid1"]],
                    "Median": [self.dat.loc[i, f"med"]]
                })

                for j in range(len(hpds)):
                    for k in range(3):
                        hpd_file.write(f"{hpds.iloc[j, k]}\t")
                    hpd_file.write("\n")

                hpd_file.write("Mode\tMidRange\tMedian\n")
                for l in range(len(pdf_info.columns)):
                    hpd_file.write(f"{pdf_info.iloc[0, l]}\t")

                hpd_file.write("\n")

            hpd_file.close()

            print("Export completed successfully.")

            print("Export completed successfully.")

    @property
    def _mixture_pdf(self):
        """
        Calculate a mixture of the reservoir age density probabilities and
        return the mode, midrange, median and HPDs.

        Returns:
            A pandas DataFrame with the mixture pdf information.
        """
        mix = np.column_stack([np.sum(self.res_age[:, 1:], axis=1)])
        # Normalize the mixture pdf
        mix /= np.sum(mix)
        # Stack the reservoir age array with the mixture pdf
        mix = np.column_stack([self.res_age[:, 0], mix])

        hpd_mix = self._hpd(mix)

        # Calculate the midrange of the HPDs
        mid1_mix = np.round(
            (hpd_mix["min"] + hpd_mix.iloc[:, 2 * len(hpd_mix) - 1]) / 2)
        # Calculate the midrange of the mixture pdf
        mid2_mix = np.round(np.mean([np.max(mix[:, 0]), np.min(mix[:, 0])]))
        # Calculate the weighted mean of the reservoir ages
        wmn_mix = np.round(np.average(mix[:, 0], weights=1 / mix[:, 1]))
        # Calculate the median of the mixture pdf
        med_mix = mix[np.max(np.where(np.cumsum(mix[:, 1]) <= 0.5)), 0]
        # Calculate the mode of the mixture pdf
        mode_mix = mix[np.where(mix[:, 1] == np.max(mix[:, 1])), 0][0]

        # Stack the mixture pdf information into a numpy array
        mix_pdf_info = np.column_stack([mode_mix, mid1_mix, med_mix])

        # Append the mixture pdf to the reservoir age array
        self.res_age = np.column_stack([self.res_age, mix[:, 1]])
        # Update the header
        header = [id_i for id_i in self.id_col]
        header.append("mixture_pdf")
        self.header = header  # resage和mixture公用的header，但是不同的数据集可能不一样

        ref = "Mixture_reservoir_age_offset"
        if self.f_14R:
            ref += 'F14R_percent_'
            # Convert the mixture pdf to delta14R in per mil
            mix_pdf_info = np.round(mix_pdf_info * 1000 / -8033, 4)
            # Convert the HPDs to delta14R in per mil
            hpd_mix.iloc[:, 0] = np.round(hpd_mix.iloc[:, 1] * 1000 / -8033, 4)
            hpd_mix.iloc[:, 1] = np.round(hpd_mix.iloc[:, 0] * 1000 / -8033, 4)
        elif self.delta_14R:
            ref += f"delta14R_permil"
            # Convert the mixture pdf to delta14R in per mil
            mix_pdf_info = np.round((np.exp(mix_pdf_info / -8033) - 1) * 1000,
                                    1)
            # Convert the HPDs to delta14R in per mil
            hpd_mix.iloc[:, 0] = np.round(
                (np.exp(hpd_mix.iloc[:, 1] / -8033) - 1) * 1000, 1)
            hpd_mix.iloc[:, 1] = np.round(
                (np.exp(hpd_mix.iloc[:, 0] / -8033) - 1) * 1000, 1)
        elif self.dcp:
            ref += f"DCP_percent"
            # Convert the mixture pdf to DCP in percent
            mix_pdf_info = np.round((1 - np.exp(mix_pdf_info / -8033)) * 100,
                                    2)
            # Convert the HPDs to DCP in percent
            hpd_mix.iloc[:, 0:2] = np.round(
                (1 - np.exp(hpd_mix.iloc[:, 0:2] / -8033)) * 100, 2)

        # Export file with information about the mixture pdf
        with open(
                f"InputData/{self.name}/{self.name}{self.cc}_{ref} Ranges.txt",
                "w") as hpd_file:
            # Write the header
            hpd_file.write(
                f"{ref}ranges at {100 * self.prob}% confidence intervals\n")
            hpd_file.write(f"Mode\tMidRange\tMedian\n")
            # Write the mixture pdf information
            for j in range(len(hpd_mix)):
                for k in range(3):
                    hpd_file.write(f"{hpd_mix.iloc[j, k]}\t")
                hpd_file.write("\n")

    @property
    def _export_resage_pdf(self):
        """Export the reservoir age offsets in the form of a probability distribution
        function (PDF) and/or histogram of posterior distribution of the reservoir age offsets.

        This function will create a CSV file with the reservoir age offsets and a text file with the
        posterior distribution of the reservoir age offsets.

        Parameters
        ----------
        self : Radical object
            Object containing all the data and methods of the Radical class.

        Returns
        -------
        None
        """
        # Find the min and max reservoir ages among all the reservoir age density probabilities
        mini = []
        maxi = []

        for i in range(len(self.dat)):
            #print(temp[f"resage_{i}"])
            mini.append(self.temp[f"resage_{i}"]["theta"].min())
            maxi.append(self.temp[f"resage_{i}"]["theta"].max())

        min_res_age, max_res_age = min(mini), max(maxi)

        # Create a matrix of zeros and paste inside all the reservoir age density
        
        res_age = np.zeros((int(np.floor(max_res_age - min_res_age + 1)),
                            1 + self.dat.shape[0]))
        res_age[:, 0] = np.arange(min_res_age, max_res_age + 1, self.radsteps)

        for i in range(len(self.dat)):
            d = np.max(np.diff(self.temp[f'uncalib_{i}']["theta"]))
            if d > 1:
                raise ValueError(
                    "Please decrease threshold parameter, e.g., threshold = 1e-14 in the command line"
                )
            res_age[int(np.floor((np.min(self.temp[f'resage_{i}'].iloc[:, 0]) - min_res_age))): int(np.floor(np.max(self.temp[f'resage_{i}'].iloc[:, 0]) - min_res_age+1)), i+1] = \
            self.temp[f'resage_{i}'].iloc[:,1]

        header = [id_i for id_i in self.id_col]
        self.header = header
        self.res_age = res_age
        ref, ref2 = "Reservoir_14C_Years", "Reservoir_age_offset_14C_yrs_"

        if self.delta_14R:
            self.res_age[:,
                         0] = 1000 * (np.exp(self.res_age[:, 0] / -8033) - 1)
            ref, ref2 = "delta14R_permil", "delta14R_"

        if self.f_14R:
            self.res_age[:, 0] = np.exp(self.res_age[:, 0] / -8033)
            ref, ref2 = "F14R", "F14R_"

        if self.dcp:
            self.res_age[:, 0] = 100 * (1 - np.exp(self.res_age[:, 0] / -8033))
            ref, ref2 = "dcp_percent", "dcp_"

        header = [ref] + self.header
        res_age = pd.DataFrame(self.res_age, columns=header)
        res_age.to_csv(
            f"InputData/{self.name}/{self.name}{self.cc}_{ref2}pdfs{self.ext}",
            index=False)
        print("\n")

        #export res_age
        for i in range(self.dat.shape[0]):

            hpds = self.temp[f"hpd_resage_{i}"]
            if hpds.shape[0] == 1:
                test = 1
            else:
                test = 2
                break

        if test == 1:
            export_resage = np.zeros((self.dat.shape[0], 5))
            ref, ref2, ref3 = "Reservoir_Min_Range", "Reservoir_Max_Range", "Reservoir_age_offset_14C_yrs_"

            for i in range(self.dat.shape[0]):
                hpds = self.temp[f"hpd_resage_{i}"]
                export_resage[i, 0] = hpds.iloc[:, 0]
                export_resage[i, 1] = hpds.iloc[:, 1]
                export_resage[i, 2] = self.dat["mode_resage"][i]
                export_resage[i, 3] = self.dat["mid1_resage"][i]
                export_resage[i, 4] = self.dat["med_resage"][i]

                if self.f_14R:
                    export_resage[i, 0] = np.round(np.exp(hpds[:, 1] / -8033),
                                                   4)
                    export_resage[i, 1] = np.round(np.exp(hpds[:, 0] / -8033),
                                                   4)
                    export_resage[i, 2] = np.round(
                        np.exp(self.dat["mode_resage"][i] / -8033), 4)
                    export_resage[i, 3] = np.round(
                        np.exp(self.dat["mid1_resage"][i] / -8033), 4)
                    export_resage[i, 4] = np.round(
                        np.exp(self.dat["med_resage"][i] / -8033), 4)
                    ref, ref2, ref3 = "F14R_Min_Range", "F14R_Max_Range", "F14R_"

                if self.delta_14R:
                    export_resage[i, 0] = np.round(
                        1000 * (np.exp(hpds[:, 1] / -8033) - 1), 1)
                    export_resage[i, 1] = np.round(
                        1000 * (np.exp(hpds[:, 0] / -8033) - 1), 1)
                    export_resage[i, 2] = np.round(
                        1000 *
                        (np.exp(self.dat["mode_resage"][i] / -8033) - 1), 1)
                    export_resage[i, 3] = np.round(
                        1000 *
                        (np.exp(self.dat["mid1_resage"][i] / -8033) - 1), 1)
                    export_resage[i, 4] = np.round(
                        1000 * (np.exp(self.dat["med_resage"][i] / -8033) - 1),
                        1)
                    ref, ref2, ref3 = "delta14R_Min_Range", "delta14R_Max_Range", "delta14R_permil_"

                if self.dcp:
                    export_resage[i, 0] = np.round(
                        100 * (1 - np.exp(hpds[:, 0] / -8033)), 2)
                    export_resage[i, 1] = np.round(
                        100 * (1 - np.exp(hpds[:, 1] / -8033)), 2)
                    export_resage[i, 2] = np.round(
                        100 * (1 - np.exp(self.dat["mode_resage"][i] / -8033)),
                        2)
                    export_resage[i, 3] = np.round(
                        100 * (1 - np.exp(self.dat["mid1_resage"][i] / -8033)),
                        2)
                    export_resage[i, 4] = np.round(
                        100 * (1 - np.exp(self.dat["med_resage"][i] / -8033)),
                        2)
                    ref, ref2, ref3 = "dcp_Min_Range", "dcp_Max_Range", "dcp_percent_"

            header = [
                f"{ref}_at_{100 * self.prob}%",
                f"{ref2}_at_{100 * self.prob}%", "Mode", "MidRange", "Median"
            ]
            export_resage_df = pd.DataFrame(export_resage, columns=header)
            output = pd.concat(
                [self.id_col, self.cal_age, self.cal_sigma, export_resage_df],
                axis=1)
            output.columns = ["id", "Calendar_age", "Calendar_age_error"
                              ] + header
            output.to_csv(
                f"InputData/{self.name}/{self.name}{self.cc}_{ref3}ranges_at_{100 * self.prob}%.csv",
                index=False)

        else:
            ref = "Reservoir_age_offset_14C_yrs_"

            if self.f_14R:
                ref = "F14R_"

            if self.delta_14R:
                ref = "delta14R_permil_"

            if self.dcp:
                ref = "dcp_percent_"

            
            hpd_file = open(
                f"InputData/{self.name}/{self.name}{self.cc}_{ref}ranges.txt",
                "w")
            hpd_file.write(
                f"{ref}ranges at {100 * self.prob}% confidence intervals\n")
            

            for i in range(self.dat.shape[0]):
                hpd_file.write(
                    f"\n\nIdentification_number: {self.dets['id'][i]}\nMinRange\tMaxRange\tprobability\n"
                )
                hpds = self.temp[f'hpd_resage_{i}']
                pdf_info = np.column_stack(
                    (self.dat['mode_resage'][

    def __repr__(self):
        s="\nHi, ResAge is here to help you with reservoir age offset calculation!\n"+\
        "The radcal function is designed to work with pairs of reservoir-derived 14C age and corresponding calendar age.\n"+\
        f"\nUncalibrating calendar ages using the {self.cc[:-4]} calibration curve and calculating reservoir age offsets..."

        return s


if __name__ == "__main__":
    A = Radical("Skiner2023", threshold=1e-14)
    A.radcal()
