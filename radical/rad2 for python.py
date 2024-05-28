import os
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp1d
import time
'''
这个版本是rad2的python版本，经过验证结果完全正确！！！
radcal等有时间再转换，源代码无法跑通
editor by 周郑鹏
联系方式：alex_chou@sjtu.edu.cn
'''


def rad2(name="example_rad2", cc=1, cc1="IntCal20.14C", cc2="SHCal13.14C", cc3="IntCal09.14C", calibrate=False, ext=".csv",
         sep=",", yrsteps=1, prob=0.95, threshold=1e-6, hpdsteps=1,
         storedat=False, export_cal_pdf=False, delta_14R = False, f_14R=False, dcp=False):
    
    dets = pd.read_csv(f"InputData/{name}/{name}{ext}", sep=sep)

    if dets.shape[1] != 5:
        raise ValueError("Your input file should contain 5 columns. Please check the template designed for the radrad function.")

    if np.min(dets.iloc[:, 1]) < 0 or np.min(dets.iloc[:, 3]) < 0:
        raise ValueError("At least one of your radiocarbon data is negative.ResAge is not designed to work with post-bomb samples. Remove such data from the input file.")

    id_col = dets["id"]
    rad_res, rad_res_sd, rad_atm, rad_atm_sd = dets["Reservoir_derived_14C_age"], dets["Reservoir_derived_14C_error"], dets["Atmosphere_derived_14C_age"], dets["Atmosphere_derived_14C_error"]
    fm_res, fm_res_sd, fm_atm, fm_atm_sd = np.exp(rad_res / -8033), np.exp(-(rad_res - rad_res_sd) / 8033) - np.exp(rad_res / -8033), np.exp(rad_atm / -8033), np.exp(-(rad_atm - rad_atm_sd) / 8033) - np.exp(rad_atm / -8033)

    # calculate sample reservoir ages
    print("Calculating reservoir age offsets...\n")
    Reservoir_age_offset = np.round(rad_res - rad_atm)
    Reservoir_age_offset_error = np.round(np.sqrt(rad_res_sd**2 + rad_atm_sd**2))

    F14R = np.round(fm_res / fm_atm, 4)
    F14R_error = np.round(F14R * np.sqrt(fm_res_sd**2 / (fm_res**2) + fm_atm_sd**2 / (fm_atm**2)), 4)
    Delta_14R = 1000 * (F14R - 1)
    Delta_14R_error = 1000 * F14R_error
    DCP = 100 * (1 - F14R)
    DCP_error = 100 * F14R_error

    output = pd.concat([dets, pd.DataFrame({"Reservoir_age_offset": Reservoir_age_offset, "Reservoir_age_offset_error": Reservoir_age_offset_error})], axis=1)

    if f_14R:
        output = pd.concat([output, pd.DataFrame({"F14R": F14R, "F14R_error": F14R_error})], axis=1)
    if delta_14R:
        output = pd.concat([output, pd.DataFrame({"delta14R": Delta_14R, "delta14R_error": Delta_14R_error})], axis=1)
    if dcp:
        output = pd.concat([output, pd.DataFrame({"DCP": DCP, "DCP_error": DCP_error})], axis=1)

    def caldist(f_cage, f_error, theta, f_mu, f_sigma, yrsteps, threshold):
        cal = pd.DataFrame({"theta": theta, "pdf": norm.pdf(f_mu, f_cage, np.sqrt(f_error**2 + f_sigma**2))})
        cal=cal[cal["pdf"]>0]
        f = interp1d(cal.iloc[:, 0], cal.iloc[:, 1])
        newx=np.arange(cal.iloc[:, 0].min(), cal.iloc[:, 0].max()+1, yrsteps)
        cal=f(newx)
        cal = pd.DataFrame(data=np.column_stack((newx, cal)), columns=["theta","pdf"])
        cal["pdf"]/= cal["pdf"].sum()
        return cal[cal["pdf"] > threshold]

    def hpd(dat, prob, hpdsteps, yrsteps):
        f=interp1d(dat["theta"],dat["pdf"])
        newx=np.arange(min(dat["theta"]),max(dat["theta"])+1,yrsteps)

        dat = f(newx)

        dat=pd.DataFrame(data=np.column_stack((newx,dat)),columns=["theta","pdf"])
        #print(dat)
        dat = dat.sort_values(by="pdf", ascending=False)
        dat["pdf"]/=dat["pdf"].sum()

        dat = dat[dat["pdf"].cumsum() <= prob]
        dat=dat.sort_values(by="theta",ascending=True)
        #s=np.where(np.diff(dat["theta"]) > hpdsteps)

        dif = np.where(np.diff(dat["theta"]) > hpdsteps)[0]

        if len(dif) == 0:
            hpds = pd.DataFrame({"min": [dat["theta"].min()], "max": [dat["theta"].max()], "prob": [100 * prob]})
        else:

            dif = np.concatenate((np.array([dat.iloc[0,0]]),np.sort(np.append(dat.iloc[dif,0],dat.iloc[dif+1,0])),dat.iloc[len(dat)-1,0].reshape(-1)))

            dif = dif.reshape(-1, 2)
            probs = [round(100 * dat[(dat["theta"] >= i[0]) & (dat["theta"] <= i[1])]["pdf"].sum(), 1) for i in dif]
            hpds = pd.DataFrame({"min": dif[:, 0], "max": dif[:, 1], "prob": probs})

        return hpds

    if calibrate:
        calcurve = pd.read_table(cc1,header=None) if cc == 1 else pd.read_table(cc2,header=None) if cc == 2 else pd.read_table(cc3,header=None) if cc == 3 else None
        ccname = "IntCal20" if cc == 1 else "SHCal13" if cc == 2 else "IntCal09" if cc == 3 else None
        cc = f"_using_{ccname}" if ccname else None

        if calcurve is not None and ccname is not None:
            print(f"Calibrating atmospheric 14C ages using the {ccname} calibration curve...\n")

            border = 0
            if any((rad_atm - rad_atm_sd) < min(calcurve.iloc[:, 1] + calcurve.iloc[:, 2])):
                if any((rad_atm + rad_atm_sd) > min(calcurve.iloc[:, 1] - calcurve.iloc[:, 2])):
                    border = 1
                else:
                    border = 2
            if any((rad_atm + rad_atm_sd) > max(calcurve.iloc[:, 1] - calcurve.iloc[:, 2])):
                if any((rad_atm - rad_atm_sd) < max(calcurve.iloc[:, 1] + calcurve.iloc[:, 2])):
                    border = 1
                else:
                    border = 2
            if border == 1:
                print("\nAt least one 14C age falls partly beyond the calibration curve and will be truncated")
            if border == 2:
                raise ValueError("\nAt least one 14C age cannot be calibrated because it falls beyond the calibration curve! Please remove those.")

            theta = calcurve.iloc[:, 0]
            f_mu = np.exp(-calcurve.iloc[:, 1] / 8033)
            f_sigma = np.exp(-(calcurve.iloc[:, 1] - calcurve.iloc[:, 2]) / 8033) - f_mu

            f_cage = np.exp(-rad_atm / 8033)
            f_error = np.exp(-(rad_atm - rad_atm_sd) / 8033) - f_cage

            dat = pd.DataFrame({"f_cage": f_cage, "f_error": f_error})
            temp=dict()
            print("Loading", end="")
            for i in range(len(f_cage)):
                calib = caldist(dat["f_cage"][i], dat["f_error"][i], theta, f_mu, f_sigma, yrsteps, threshold)
                #print(calib.shape)
                temp[f"calib_{i}"] = calib["pdf"]
                temp[f"hpd_{i}"] = hpd(calib, prob=prob, hpdsteps=hpdsteps, yrsteps=yrsteps)
                dat.loc[i,"mid1"] = round((temp[f"hpd_{i}"]["min"].values[0] + temp[f"hpd_{i}"]["max"].values[0]) / 2)
                yrs = calib["theta"]
                dat.loc[i,"mid2"] = round(np.mean([max(yrs), min(yrs)]))
                dat.loc[i,"wmn"] = round(np.average(calib["theta"], weights=1/calib["pdf"]))
                dat.loc[i,f"med"] =calib[calib["pdf"].cumsum()<=0.5]["theta"].iloc[-1]
                dat.loc[i,"mode"]=calib[calib["pdf"]==calib["pdf"].max()]["theta"].values
                print(".", end='', flush=True)
                #time.sleep(0.1)
                #dat[f"mode_{i}"] = calib[calib["pdf"].idxmax()]["theta"]


            if storedat:
                globals()["dat"] = dat

            if export_cal_pdf:
                print("Exporting your calibrated ages...")
                mini=[]
                maxi=[]
                for i in range(len(f_cage)):
                    a=temp[f"calib_{i}"]
                    mini.append(min(temp[f"calib_{i}"].index))
                    maxi.append(max(temp[f"calib_{i}"].index))
                min_val = min(mini)
                max_val = max(maxi)
                header=[]

                cal_age = pd.DataFrame({"cal_year_BP": np.arange(min_val, max_val + 1, yrsteps)})
                for i in range(len(dat)):
                    #d = max(np.diff(temp[f"calib_{i}"].dropna().index))
                    d = max(np.diff(pd.DataFrame(temp[f"calib_{i}"]).index))
                    if d > 1:
                        raise ValueError("Please, decrease the threshold parameter: e.g. threshold = 1e-14 in the command line.")
                    cal_age[dat.columns[i + 2]] = temp[f"calib_{i}"].reindex(cal_age.index).interpolate(method='linear', axis=0)
                    header.append(id_col[i])
                cal_age.columns=np.append("cal_year_BP",header)


                cal_age.to_csv(f"InputData/{name}/{name}{cc}_Calibrated_ages_pdfs{ext}", index=False)

            for i in range(len(f_cage)):
                hpds = temp[f"hpd_{i}"]
                test = 1 if len(hpds) == 1 else 2
                if test == 2:
                    break

            if test == 1:
                for i in range(len(dat)):
                    export_calib = pd.DataFrame({"CalibratedMinRange_at_{}%".format(100 * prob): hpds["min"].values,
                                                    "CalibratedMaxRange_at_{}%".format(100 * prob): hpds["max"].values,
                                                    "Mode": dat[f"mode"][i],
                                                    "MidRange": dat[f"mid1"][i],
                                                    "Median": dat[f"med"][i]})

                colnames = [f"{id_col[i]}_{col}" for col in export_calib.columns]
                export_calib.columns = colnames

                cal_output = pd.concat([pd.DataFrame(id_col),export_calib], axis=1)
                cal_output.to_csv(f"InputData/{name}/{name}{cc}_Calibrated_age_ranges_at_{100 * prob}%{ext}", index=False)
            else:
                ref = "Calibrated_"
                hpd_file = open(f"InputData/{name}/{name}{cc}_{ref}ranges.txt", "w")
                hpd_file.write(f"{ref}ranges at {100 * prob}% confidence intervals\n")

                for i in range(len(f_cage)):
                    hpd_file.write(
                        f"\n\nIdentification_number: {id_col[i]}\nCal_year_MinRange\tCal_year_MaxRange\tprobability\n")
                    hpds = temp[f"hpd_{i}"]
                    pdf_info = pd.DataFrame(
                        {"Mode": [dat.loc[i,f"mode"]], "MidRange": [dat.loc[i,f"mid1"]], "Median": [dat.loc[i,f"med"]]})

                    for j in range(len(hpds)):
                        for k in range(3):
                            hpd_file.write(f"{hpds.iloc[j, k]}\t")
                        hpd_file.write("\n")

                    hpd_file.write("Mode\tMidRange\tMedian\n")
                    for l in range(len(pdf_info.columns)):
                        hpd_file.write(f"{pdf_info.iloc[0, l]}\t")

                    hpd_file.write("\n")

                hpd_file.close()

        # export output as .csv file
        output.to_csv(f"InputData/{name}/{name}_out{ext}", index=False)

    # Done message
    print("\n\nWork has been done successfully, check outputs in your folder.\n\n")

# list the available data
data = os.listdir("InputData/")
print(data)

# Welcome
print("Hi, the ResAge package is here to help you with reservoir age calculation!")
print("The rad2 function is designed to work with pairs of reservoir-derived and atmospheric 14C ages.\n\n")

if __name__ == "__main__":
    rad2(calibrate=True)
