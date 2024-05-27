import numpy as np
from scipy.stats import norm
import pandas as pd
from scipy.interpolate import interp1d
import os
import time
from tqdm import tqdm
from colorama import Fore
'''
editor by Alex_chou
mail：alex_chou@sjtu.edu.cn
'''


def radcal(name="example_radcal_Southon2012", cc=1, cc1="IntCal20.14C", cc2="SHCal13.14C", cc3="IntCal13.14C",
           ext=".csv", sep=",", yrsteps=1, uncalsteps=5, radsteps=1, convolsteps=1, prob=0.95, threshold=1e-6,
           hpdsteps=1, storedat=False, export_uncal_pdf=False, mixture_pdf=False, delta_14R=False, dcp=False, f_14R=False):
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
    """
    _radcal(name, cc, cc1, cc2, cc3, ext, sep, yrsteps, uncalsteps, radsteps, convolsteps, prob, threshold,
            hpdsteps, storedat, export_uncal_pdf, mixture_pdf, delta_14R, dcp, f_14R)

def _radcal(name, cc, cc1, cc2, cc3, ext, sep, yrsteps, uncalsteps, radsteps, convolsteps, prob, threshold,
            hpdsteps, storedat, export_uncal_pdf, mixture_pdf, delta_14R, dcp, f_14R):
    dets = pd.read_csv(f"InputData/{name}/{name}{ext}", sep=sep)



    # check for consistency
    if dets.shape[1] != 5:
        raise ValueError("Your input file should contain 5 columns.\nPlease, check the template designed for the radcal function.")

    if min(dets.iloc[:, 1]) < 0 or min(dets.iloc[:, 3]) < 0:
        raise ValueError("At least one of your radiocarbon data or one of your calendar age is negative."
                         "\nResAge is not designed to work with post-bomb samples."
                         "\nRemove such data from the input file.")

    if delta_14R and (dcp or f_14R or dcp or f_14R or dcp and f_14R and delta_14R):
        raise ValueError("Please, choose only one metric: reservoir age offset (default), f.14R, delta.14R, or dcp")

    # prepare for reservoir age calculation
    id_col = dets["id"]
    rad_res = dets["Reservoir_derived_14C_age"]
    rad_res_sd = dets["Reservoir_derived_14C_error"]
    cal_age = dets["Calendar_age"]
    cal_sigma = dets["Calendar_age_error"]




    def uncaldist(cal_age, cal_sigma, theta, f_mu, f_sigma, yrsteps, threshold, uncalsteps, radsteps):
        """

        :param cal_age: the calculate age
        :param cal_sigma: the calculate 方差
        :param theta: θ值
        :param f_mu: 平均值
        :param f_sigma: 方差
        :param yrsteps: 年分辨率
        :param threshold: 阈值
        :param uncalsteps: 未矫正的年龄分辨率
        :param radsteps: 储库年分辨率
        :return:
        """
        # Let's reduce calculation time

        # 1) work over 6 time sigma of the sample calendar age
        # Reduce calculation time
        cal_span = theta[(theta >= cal_age - 6 * cal_sigma) & (theta <= cal_age + 6 * cal_sigma)]
        # 2) select to suitable range of the radiocarbon time scale to calibrate and convert it into F14C
        rad_span = np.round(np.arange(min(rad[cal_span]) - 50 * max(rad_sigma[cal_span]),
                                      max(rad[cal_span]) + 50 * max(rad_sigma[cal_span])+1, uncalsteps))
        f_span = np.exp(-rad_span / 8033)

        # Find how far each calendar time of cal.span is from cal.age
        LKH_cal_age = np.column_stack((cal_span, norm.pdf(cal_span, cal_age, cal_sigma)))

        # Uncalibrate
        uncal = []

        for i in range(len(f_span)):
            # How far is a single F14C from the calibration curve f.mu
            LKH_cc = norm.pdf(f_mu[cal_span], f_span[i], f_sigma[cal_span])
            LKH_ri = yrsteps * np.sum(LKH_cal_age[:, 1] * LKH_cc)
            uncal.append(LKH_ri)

        uncal = np.column_stack((rad_span, uncal))

        # Interpolate and normalise uncalibrated distribution to 1
        uncal = uncal[uncal[:, 1] > 0, :]  # remove unnecessary data
        f = interp1d(uncal[:, 0], uncal[:, 1])


        uncal = np.column_stack((np.arange(min(uncal[:, 0]), max(uncal[:, 0])+1, radsteps),
                                 f(np.arange(min(uncal[:, 0]), max(uncal[:, 0])+1, radsteps)) / np.sum(
                                     f(np.arange(min(uncal[:, 0]), max(uncal[:, 0])+1, radsteps)))))
        #print(uncal[uncal[:, 1] > threshold, :].shape)

        # Only report those normalised uncalibrated probabilities beyond a threshold
        uncal=pd.DataFrame(data=uncal,columns=["theta","pdf"])
        return uncal[uncal["pdf"] > threshold]

    # Find reservoir age: subtract uncalibrated age to reservoir-derived 14C age
    def convolution(dat, rad_res, rad_res_sd, threshold, convolsteps, radsteps):
        """

        :param dat:
        :param rad_res:
        :param rad_res_sd:
        :param threshold:
        :param convolsteps:
        :param radsteps:
        :return:
        """
        # prepare for convolution
        # create the normal distribution of the reservoir-derived 14C age
        spl_span = np.arange(rad_res - 10 * rad_res_sd, rad_res + 10 * rad_res_sd+1, radsteps)
        spl_pdf = np.column_stack((spl_span,norm.pdf(spl_span, rad_res, rad_res_sd)))
        spl_pdf[:,1] /= np.sum(spl_pdf[:,1])
        spl_pdf = spl_pdf[spl_pdf[:,1] > threshold,:]

        # move the end of matrix dat (uncalibrated age) to the start of reservoir 14C age spl.pdf
        new_origin = np.min(spl_pdf[:,0]) - np.max(dat.iloc[:, 0])
        uncal_conv = np.column_stack((dat.iloc[:, 0] + new_origin - 1, dat.iloc[:, 1]))

        # make matrices same dimension so that they overlap
        m, n = spl_pdf.shape[0], uncal_conv.shape[0]
        s=np.column_stack((uncal_conv[:, 0], np.zeros(n)))

        #print(uncal_conv[:, 1])
        spl_pdf_new = np.row_stack((s, spl_pdf))
        o=np.column_stack((spl_pdf[:,0], np.zeros(m)))
        uncal_conv_new = np.row_stack((uncal_conv, o))
        spl_density, uncal_density = spl_pdf_new[:, 1], uncal_conv_new[:, 1]

        # convolution !
        LKH_res = []

        for i in range(0, n + m + n, convolsteps):

            spl_density = np.append(spl_density, np.zeros(convolsteps))

            uncal_density = np.append(np.zeros(convolsteps), uncal_density)

            LKH_res.append(np.sum(spl_density * uncal_density))

        res_age = np.column_stack((new_origin - 1 + np.arange(1, n + m + n + 1, convolsteps), LKH_res))

        res_age = res_age[res_age[:, 1] > 0, :]
        res_age=pd.DataFrame(data=res_age,columns=["theta","pdf"])

        return res_age[res_age["pdf"] > threshold]




    def hpd(dat, prob, hpdsteps, yrsteps):
        dat=pd.DataFrame(data=dat,columns=["theta","pdf"])
        f = interp1d(dat["theta"], dat["pdf"])
        newx = np.arange(min(dat["theta"]), max(dat["theta"]) + 1, yrsteps)

        dat = f(newx)

        dat = pd.DataFrame(data=np.column_stack((newx, dat)), columns=["theta", "pdf"])
        # print(dat)
        dat = dat.sort_values(by="pdf", ascending=False)
        dat["pdf"] /= dat["pdf"].sum()

        dat = dat[dat["pdf"].cumsum() <= prob]
        dat = dat.sort_values(by="theta", ascending=True)
        # s=np.where(np.diff(dat["theta"]) > hpdsteps)

        dif = np.where(np.diff(dat["theta"]) > hpdsteps)[0]

        if len(dif) == 0:
            hpds = pd.DataFrame({"min": [dat["theta"].min()], "max": [dat["theta"].max()], "prob": [100 * prob]})
        else:

            dif = np.concatenate((
                                 np.array([dat.iloc[0, 0]]), np.sort(np.append(dat.iloc[dif, 0], dat.iloc[dif + 1, 0])),
                                 dat.iloc[len(dat) - 1, 0].reshape(-1)))

            dif = dif.reshape(-1, 2)
            probs = [round(100 * dat[(dat["theta"] >= i[0]) & (dat["theta"] <= i[1])]["pdf"].sum(), 1) for i in dif]
            hpds = pd.DataFrame({"min": dif[:, 0], "max": dif[:, 1], "prob": probs})

        return hpds


############################################################################################################
#######################################################################################################
    # Now, let's do the job now

        # Read calibration curve based on cc value
    now=time.time()
    if cc == 1:
        calcurve = pd.read_table(cc1, header=None)
    elif cc == 2:
        calcurve = pd.read_table(cc2, header=None)
    elif cc == 3:
        calcurve = pd.read_table(cc3, header=None)
    else:
        raise ValueError("Invalid value for cc. Please check the manual.")

# Set ccname based on cc value
    if cc == 1:
        ccname = "IntCal20"
    elif cc == 2:
        ccname = "SHCal13"
    elif cc == 3:
        ccname = "IntCal13"
    else:
        raise ValueError("Invalid value for cc.")

# Set cc based on cc value
    if cc == 1:
        cc = "_using_IntCal20"
    elif cc == 2:
        cc = "_using_SHCal13"
    elif cc == 3:
        cc = "_using_IntCal13"
    else:
        raise ValueError("Invalid value for cc.")

    print(f"\nUncalibrating calendar ages using the {ccname} calibration curve and calculating reservoir age offsets...")

# Prepare calcurve for calculation in F14C and interpolate
    theta = calcurve.iloc[:, 0]
    f_mu = np.exp(-calcurve.iloc[:, 1] / 8033)
    f_sigma = np.exp(-(calcurve.iloc[:, 1] - calcurve.iloc[:, 2]) / 8033) - f_mu

# Interpolate calibration each yrsteps (required for uncalibration function)
    theta_interp = np.arange(0, theta.max()+1, yrsteps)
    rad = np.round(interp1d(calcurve.iloc[:, 0],calcurve.iloc[:, 1])(theta_interp))
    rad_sigma = np.round(interp1d(calcurve.iloc[:, 0],calcurve.iloc[:, 2])(theta_interp))
    f_mu = interp1d(calcurve.iloc[:, 0], f_mu, kind='linear')(theta_interp)
    f_sigma = interp1d(calcurve.iloc[:, 0], f_sigma, kind='linear')(theta_interp)

    dat = pd.DataFrame({
    'id':id_col,
    'cal_age': cal_age,
    'cal_sigma': cal_sigma,
    'rad_res': rad_res,
    'rad_res_sd': rad_res_sd
    })

    temp=dict()
    #print("\n.running")
    iter_now = time.time()
    for i in tqdm(range(len(dat['cal_age'])),desc="processing",total=len(dat['cal_age']),colour="blue"):
        #print(".")
        # Implement or replace uncal_dist function

        uncalib = uncaldist(dat['cal_age'][i], dat['cal_sigma'][i], theta_interp, f_mu, f_sigma, yrsteps, threshold, uncalsteps, radsteps)



        temp[f"uncalib_{i}"] = uncalib  # resulting uncalibrated (atm radiocarbon distribution) distribution
        temp[f"hpd_uncal_{i}"] = hpd(uncalib, prob, hpdsteps, radsteps)  # uncalibrated highest probability density corresponding to prob
        #print(temp['hpd_uncal_{i}'][0])
        dat.loc[i,"mid1"] = round((temp[f"hpd_uncal_{i}"]["min"].values[0] + temp[f"hpd_uncal_{i}"]["max"].values[0]) / 2)  # midpoints of uncalibrated ranges
        radyrs = uncalib["theta"]
        dat.loc[i,"mid2"] = round(np.mean([max(radyrs), min(radyrs)]))  # midpoints of entire uncalibrated distributions (with probabilities beyond threshold)
        dat.loc[i,"wmn"] = round(np.average(uncalib["theta"], weights=1/uncalib["pdf"]))  # weighted means of uncalibrated ranges
        dat.loc[i,"med"] =uncalib[uncalib["pdf"].cumsum()<=0.5]["theta"].iloc[-1]  # medians of uncalibrated distributions
        dat.loc[i,"mode"]=uncalib[uncalib["pdf"]==uncalib["pdf"].max()]["theta"].values  # maximum densities of uncalibrated distributions

        # Implement or replace convolution function
        #print(dat['rad_res'][i])
        resage = convolution(uncalib, dat['rad_res'][i], dat['rad_res_sd'][i], threshold, convolsteps, radsteps)
        temp[f'resage_{i}'] = resage  # resulting reservoir age

        temp[f"hpd_resage_{i}"] = hpd(resage, prob, hpdsteps, radsteps)  # reservoir age highest probability density corresponding to prob

        dat.loc[i,'mid1_resage'] = round((temp[f"hpd_resage_{i}"]["min"].values[0] + temp[f"hpd_resage_{i}"]["max"].values[0]) / 2)  # midpoints of reservoir ranges
        resyrs = resage.iloc[:, 0]
        dat.loc[i,'mid2_resage'] = round(np.mean([max(resyrs), min(resyrs)]))  # midpoints of entire reservoir age distributions (with probabilities beyond threshold)
        dat.loc[i,'wmn_resage'] = round(np.average(resage.iloc[:,0], weights=1/resage.iloc[:,1]))  # weighted means of reservoir ranges
        dat.loc[i,'med_resage'] = resage[resage["pdf"].cumsum()<=0.5]["theta"].iloc[-1]  # medians of reservoir distributions
        dat.loc[i,'mode_resage'] = resage[resage["pdf"]==resage["pdf"].max()]["theta"].values  # maximum densities of reservoir distributions

        if (i+1)%100 == 0:

            #print(".", end='')
            print(f"\t已经完成了{i+1}/{len(dat)}",end="")
            print(f"\tthe iter spending time is {round((time.time()-iter_now)//60)}min{round(time.time()-iter_now-60*((time.time()-iter_now)//60),1)}s ")
            iter_now=time.time()

    print(f"the total spending time is {(time.time()-now)//60}min{round(time.time()-now-60*((time.time()-now)//60),1)}s ")




    dat.to_csv(f"InputData/{name}/{name}{cc}_output.csv", index=False)

    print("Work has been done successfully, check outputs in your folder.\n")
    if storedat:
        globals()["dat"]=dat

    ###################################################################################################################3
    

# Export the uncalibrated probability density if required

    if export_uncal_pdf:
        print("Exporting your uncalibrated ages")
        mini = []
        maxi = []
        for i in range(len(dat)):
            mini.append(min(temp[f"uncalib_{i}"].index))
            maxi.append(max(temp[f"uncalib_{i}"].index))
        min_val = min(mini)
        max_val = max(maxi)
        header = []

        cal_age = pd.DataFrame({"uncal_year_BP": np.arange(min_val, max_val + 1, yrsteps)})
        for i in range(len(dat)):
            # d = max(np.diff(temp[f"calib_{i}"].dropna().index))
            d = max(np.diff(pd.DataFrame(temp[f"uncalib_{i}"]).index))
            if d > 1:
                raise ValueError(
                    "Please, decrease the threshold parameter: e.g. threshold = 1e-14 in the command line.")
            cal_age[dat.columns[i + 2]] = temp[f"uncalib_{i}"].reindex(cal_age.index).interpolate(method='linear', axis=0)
            header.append(id_col[i])
        cal_age.columns = np.append("uncal_year_BP", header)

        cal_age.to_csv(f"InputData/{name}/{name}{cc}_Uncalibrated_ages_pdfs{ext}", index=False)

        for i in range(len(dat)):
            hpds = temp[f"hpd_{i}"]
            test = 1 if len(hpds) == 1 else 2
            if test == 2:
                break

        if test == 1:
            for i in range(len(dat)):
                export_calib = pd.DataFrame({"UncalibratedMinRange_at_{}%".format(100 * prob): hpds["min"].values,
                                             "UncalibratedMaxRange_at_{}%".format(100 * prob): hpds["max"].values,
                                             "Mode": dat[f"mode"][i],
                                             "MidRange": dat[f"mid1"][i],
                                             "Median": dat[f"med"][i]})

            colnames = [f"{id_col[i]}_{col}" for col in export_calib.columns]
            export_calib.columns = colnames

            cal_output = pd.concat([pd.DataFrame(id_col), export_calib], axis=1)
            cal_output.to_csv(f"InputData/{name}/{name}{cc}_Uncalibrated_age_ranges_at_{100 * prob}%{ext}", index=False)
        else:
            ref = "Uncalibrated_"
            hpd_file = open(f"InputData/{name}/{name}{cc}_{ref}ranges.txt", "w")
            hpd_file.write(f"{ref}ranges at {100 * prob}% confidence intervals\n")

            for i in range(len(dat)):
                hpd_file.write(
                    f"\n\nIdentification_number: {id_col[i]}\nCal_year_MinRange\tCal_year_MaxRange\tprobability\n")
                hpds = temp[f"hpd_{i}"]
                pdf_info = pd.DataFrame(
                    {"Mode": [dat.loc[i, f"mode"]], "MidRange": [dat.loc[i, f"mid1"]], "Median": [dat.loc[i, f"med"]]})

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

    ########################################################################################################################################
    # Export reservoir age distributions
    print("Exporting reservoir age offsets...")

    # Find the min and max reservoir ages among all the reservoir age density probabilities

    mini=[]
    maxi=[]

    for i in range(len(dat)):
        #print(temp[f"resage_{i}"])
        mini.append(temp[f"resage_{i}"]["theta"].min())
        maxi.append(temp[f"resage_{i}"]["theta"].max())


    min_res_age, max_res_age = min(mini), max(maxi)

    # Create a matrix of zeros and paste inside all the reservoir age densities
    #s=np.floor(max_res_age - min_res_age + 1)
    res_age = np.zeros((int(np.floor(max_res_age - min_res_age + 1)), 1 + dat.shape[0]))
    res_age[:, 0] = np.arange(min_res_age, max_res_age + 1, radsteps)


    for i in range(len(dat)):
        d = np.max(np.diff(temp[f'uncalib_{i}']["theta"]))
        if d > 1:
            raise ValueError("Please decrease threshold parameter, e.g., threshold = 1e-14 in the command line")
        res_age[int(np.floor((np.min(temp[f'resage_{i}'].iloc[:, 0]) - min_res_age))): int(np.floor(np.max(temp[f'resage_{i}'].iloc[:, 0]) - min_res_age+1)), i+1] = temp[f'resage_{i}'].iloc[:,1]


    header = [id_i for id_i in id_col]
    ###########################################################################################################################
        



# Mixture PDF
    if mixture_pdf:
        mix = np.column_stack([np.sum(res_age[:, 1:], axis=1)])
        mix /= np.sum(mix)
        mix = np.column_stack([res_age[:, 0], mix])

        hpd_mix = hpd(mix, prob=prob, hpdsteps=hpdsteps, yrsteps=radsteps)

        mid1_mix = np.round((hpd_mix["min"] + hpd_mix.iloc[:,2 * len(hpd_mix)-1]) / 2)
        mid2_mix = np.round(np.mean([np.max(mix[:, 0]), np.min(mix[:, 0])]))
        wmn_mix = np.round(np.average(mix[:, 0], weights=1/mix[:, 1]))
        med_mix = mix[np.max(np.where(np.cumsum(mix[:, 1]) <= 0.5)), 0]
        mode_mix = mix[np.where(mix[:, 1] == np.max(mix[:, 1])), 0][0]

        mix_pdf_info = np.column_stack([mode_mix, mid1_mix, med_mix])

        res_age = np.column_stack([res_age, mix[:, 1]])
        header.append("mixture_pdf")

        ref = f"Mixture_reservoir_age_offset_{'F14R_percent_' if f_14R else 'delta14R_permil_' if delta_14R else 'dcp_percent_'}"

        if any([f_14R, delta_14R, dcp]):
            mix_pdf_info = np.round(mix_pdf_info * 1000 / -8033, 4) if f_14R else np.round((np.exp(mix_pdf_info / -8033) - 1) * 1000, 1) if delta_14R else np.round((1 - np.exp(mix_pdf_info / -8033)) * 100, 2)

        hpd_mix.iloc[:, 0:2] = np.round(hpd_mix.iloc[:, 0:2] * 1000 / -8033, 4) if f_14R else np.round((np.exp(hpd_mix.iloc[:, 0:2] / -8033) - 1) * 1000, 1) if delta_14R else np.round((1 - np.exp(hpd_mix.iloc[:, 0:2] / -8033)) * 100, 2)

        # Export file with information about the mixture pdf
        hpd_file = f"InputData/{name}/{name}{cc}_{ref}Ranges.txt"
        np.savetxt(hpd_file, np.column_stack([hpd_mix]), fmt='%s', delimiter='\t')

        mix_pdf_info_file = f"InputData/{name}/{name}{cc}_{ref}Info.txt"
        np.savetxt(mix_pdf_info_file, mix_pdf_info, fmt='%s', delimiter='\t')

    # Reservoir age offset, F14R, delta14R, or dcp metrics
    ref, ref2 = "Reservoir_14C_years", "Reservoir_age_offset_"

    if delta_14R:
        res_age[:, 0] = 1000 * (np.exp(res_age[:, 0] / -8033) - 1)
        ref, ref2 = "delta14R_permil", "delta14R_"

    if f_14R:
        res_age[:, 0] = np.exp(res_age[:, 0] / -8033)
        ref, ref2 = "F14R", "F14R_"

    if dcp:
        res_age[:, 0] = 100 * (1 - np.exp(res_age[:, 0] / -8033))
        ref, ref2 = "dcp_percent", "dcp_"

        header = [ref] + header
        res_age = pd.DataFrame(res_age, columns=header)
        res_age.to_csv(f"InputData/{name}/{name}{cc}_{ref2}pdfs{ext}", index=False)
    print("\n")
    ##########################################################################################


# Export Reservoir age ranges of all dates


    for i in range(dat.shape[0]):

        hpds = temp[f"hpd_resage_{i}"]
        if hpds.shape[0] == 1:
            test = 1
        else:
            test = 2
            break

    if test == 1:
        export_resage = np.zeros((dat.shape[0], 5))
        ref, ref2, ref3 = "Reservoir_Min_Range", "Reservoir_Max_Range", "Reservoir_age_offset_14C_yrs_"



        for i in range(dat.shape[0]):
            hpds = temp[f"hpd_resage_{i}"]
            export_resage[i, 0] = hpds.iloc[:, 0]
            export_resage[i, 1] = hpds.iloc[:, 1]
            export_resage[i, 2] = dat["mode_resage"][i]
            export_resage[i, 3] = dat["mid1_resage"][i]
            export_resage[i, 4] = dat["med_resage"][i]

            if f_14R:
                export_resage[i, 0] = np.round(np.exp(hpds[:, 2] / -8033), 4)
                export_resage[i, 1] = np.round(np.exp(hpds[:, 1] / -8033), 4)
                export_resage[i, 2] = np.round(np.exp(dat["mode_resage"][i] / -8033), 4)
                export_resage[i, 3] = np.round(np.exp(dat["mid1_resage"][i] / -8033), 4)
                export_resage[i, 4] = np.round(np.exp(dat["med_resage"][i] / -8033), 4)
                ref, ref2, ref3 = "F14R_Min_Range", "F14R_Max_Range", "F14R_"

            if delta_14R:
                export_resage[i, 0] = np.round(1000 * (np.exp(hpds[:, 2] / -8033) - 1), 1)
                export_resage[i, 1] = np.round(1000 * (np.exp(hpds[:, 1] / -8033) - 1), 1)
                export_resage[i, 2] = np.round(1000 * (np.exp(dat["mode_resage"][i] / -8033) - 1), 1)
                export_resage[i, 3] = np.round(1000 * (np.exp(dat["mid1_resage"][i] / -8033) - 1), 1)
                export_resage[i, 4] = np.round(1000 * (np.exp(dat["med_resage"][i] / -8033) - 1), 1)
                ref, ref2, ref3 = "delta14R_Min_Range", "delta14R_Max_Range", "delta14R_permil_"

            if dcp:
                export_resage[i, 0] = np.round(100 * (1 - np.exp(hpds[:, 1] / -8033)), 2)
                export_resage[i, 1] = np.round(100 * (1 - np.exp(hpds[:, 2] / -8033)), 2)
                export_resage[i, 2] = np.round(100 * (1 - np.exp(dat["mode_resage"][i] / -8033)), 2)
                export_resage[i, 3] = np.round(100 * (1 - np.exp(dat["mid1_resage"][i] / -8033)), 2)
                export_resage[i, 4] = np.round(100 * (1 - np.exp(dat["med_resage"][i] / -8033)), 2)
                ref, ref2, ref3 = "dcp_Min_Range", "dcp_Max_Range", "dcp_percent_"

        header = [f"{ref}_at_{100 * prob}%", f"{ref2}_at_{100 * prob}%", "Mode", "MidRange", "Median"]
        export_resage_df = pd.DataFrame(export_resage, columns=header)
        output = pd.concat([id_col, cal_age, cal_sigma, export_resage_df], axis=1)
        output.columns = ["id", "Calendar_age", "Calendar_age_error"] + header
        output.to_csv(f"InputData/{name}/{name}{cc}_{ref3}ranges_at_{100 * prob}%.csv", index=False)

    else:
        ref = "Reservoir_age_offset_14C_yrs_"

        if f_14R:
            ref = "F14R_"

        if delta_14R:
            ref = "delta14R_permil_"

        if dcp:
            ref = "dcp_percent_"

        #hpd_file=pd.DataFrame(data=np.zeros([len(dets),7]))
        hpd_file = open(f"InputData/{name}/{name}{cc}_{ref}ranges.txt", "w")
        hpd_file.write(f"{ref}ranges at {100 * prob}% confidence intervals\n")
        #hpd_file.columns = ["id", "MinRange", "MaxRange", "probability", "Mode", "MidRange", "Median"]
        #hpd_file["id"]=dets["id"]

        for i in range(dat.shape[0]):

            hpd_file.write(f"\n\nIdentification_number: {dets['id'][i]}\nMinRange\tMaxRange\tprobability\n")
            hpds = temp[f'hpd_resage_{i}']
            pdf_info = np.column_stack((dat['mode_resage'][i], dat['mid1_resage'][i],dat['med_resage'][i]))
            pdf_info=pd.DataFrame(data=pdf_info)

            if f_14R:
                hpds_save = hpds.iloc[:, 0].copy()
                hpds.iloc[:, 0] = np.round(np.exp(hpds.iloc[:, 1] / -8033), 4)
                hpds.iloc[:, 1] = np.round(np.exp(hpds_save / -8033), 4)
                pdf_info = np.round(np.exp(pdf_info / -8033), 4)

            if delta_14R:
                hpds_save = hpds.iloc[:, 0].copy()
                hpds.iloc[:, 0] = np.round(1000 * (np.exp(hpds.iloc[:, 1] / -8033) - 1), 1)
                hpds.iloc[:, 1] = np.round(1000 * (np.exp(hpds_save / -8033) - 1), 1)
                pdf_info = np.round(1000 * (np.exp(pdf_info / -8033) - 1), 1)

            if dcp:
                hpds.iloc[:, 0] = np.round(100 * (1 - np.exp(hpds.iloc[:, 0] / -8033)), 2)
                hpds.iloc[:, 1] = np.round(100 * (1 - np.exp(hpds.iloc[:, 1] / -8033)), 2)
                pdf_info = np.round(100 * (1 - np.exp(pdf_info / -8033)), 2)

            for j in range(hpds.shape[0]):
                #hpd_file.iloc[i*hpds.shape[0]+j, 1]=hpds.iloc[j, 0]
                #hpd_file.iloc[i*hpds.shape[0]+j, 2] = hpds.iloc[j, 1]
                #hpd_file.iloc[i*hpds.shape[0]+j, 3] = hpds.iloc[j, 2]
                hpd_file.write(f"{hpds.iloc[j, 0]}\t{hpds.iloc[j, 1]}\t{hpds.iloc[j, 2]}\n")

            hpd_file.write("Mode\tMidRange\tMedian\n")
            for l in range(pdf_info.shape[0]):
                #hpd_file.iloc[i*pdf_info.shape[0]+l, 4] = pdf_info.iloc[l, 0]
                #hpd_file.iloc[i*pdf_info.shape[0]+l, 5] = pdf_info.iloc[l, 1]
                #hpd_file.iloc[i*pdf_info.shape[0]+l, 6] = pdf_info.iloc[l, 2]

                hpd_file.write(f"{pdf_info.iloc[l, 0]}\t{pdf_info.iloc[l, 1]}\t{pdf_info.iloc[l, 2]}\n")

        hpd_file.close()
        #hpd_file.(f"InputData/{name}/{name}{cc}_{ref}ranges_at_{100 * prob}%.csv", index=False)

    
    # Done message
    print("Work has been done successfully, check outputs in your folder.\n\n")

# List the available data
Data = os.listdir("InputData/")

# Welcome
print("\nHi, ResAge is here to help you with reservoir age offset calculation!\n")
print("The radcal function is designed to work with pairs of reservoir-derived 14C age and corresponding calendar age.\n")

if __name__ == "__main__":
    radcal("Skiner2023",threshold=1e-14)





    