import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plotspec",
]


#
# Function to display the information stored in the .spec file.
# For a given galaxy, a .spec file is created with a lot of information concerning
# the input data as well as the result of the fit.
# filename in input: is the name of one .spec file.
#
# The first plot shows the observed fluxes as a function of wavelength,
# as well as the best-fit templates on the top of it
# 6 classes could be displayed
#  GAL-1     First galaxy solution
#  GAL-2     Second galaxy solution,
#  GAL-FIR   Far-infrared library like Dale, etc
#  GAL-STOCH The stochastic library (not used anymore)
#  QSO       QSO/AGN library
#  STAR      Star library
#
# The inset panel shows the PDF(z) with two different definitions of the PDF(z)
# - Bayesian posterior
# - Profile likelihood
#
def plotspec(filename):
    ##############  READING INFORMATION  ###############

    ### Open .spec file[s] and read the parameters
    with open(filename, "r") as fsp:
        # print("File:", filename)

        # Information stored in lines 1-2 (Id, spec-z, photo-z)
        bid = fsp.readline()  # header1
        line = fsp.readline()
        line = line.split()
        id = line[0]
        zspec = float(line[1])
        zphot = float(line[2])

        # Information stored in lines 3-4 (number of filters)
        bid = fsp.readline()  # header2 (seems a wrong header create by LePHARE)
        line = fsp.readline()
        line = line.split()
        nfilt = int(line[1])

        # Information stored in lines 5-6 (number of PDF steps)
        bid = fsp.readline()  # header3
        line = fsp.readline()
        line = line.split()
        npdf = int(line[1])

        # Information stored in lines 7-13 (Information on the best-fit templates)
        bid = fsp.readline()
        # header4:  Type Nline Model Library Nband    Zphot Zinf Zsup Chi2  PDF
        #     Extlaw EB-V Lir Age  Mass SFR SSFR
        models_info = []
        # Loop over GAL-1, GAL-2, GAL-FIR, GAL-STOCH, QSO, STAR
        for _ in range(6):
            line = fsp.readline()
            model = line.split()
            # information on the best-fit template, as the number of lines used after
            models_info.append(model)

        # Read observed mag, err_mag, lambda_eff and width, modeled magnitude
        mag = np.zeros(nfilt)
        em = np.zeros(nfilt)
        lf = np.zeros(nfilt)
        dlf = np.zeros(nfilt)
        mmod = np.zeros(nfilt)
        mfir = np.zeros(nfilt)
        mphys = np.zeros(nfilt)
        for i in range(nfilt):  # Loop over the filters
            line = fsp.readline()
            line = line.split()
            mag[i] = float(line[0])  # Mag
            em[i] = float(line[1])  # Error mag
            lf[i] = float(line[2])  # Effective wavelength of the filter in A
            dlf[i] = float(line[3])  # Width of the filter in A
            mmod[i] = float(line[4])
            mfir[i] = float(line[5])
            mphys[i] = float(line[6])
            mmod[i] = float(line[8])  # Predicted magnitude from the best-fit galaxy template

        # convert mag(AB syst.) in log(flux)
        ibad = np.where((mmod <= 0) | (mmod > 45))
        mmod = -0.4 * (mmod - 23.91)  # uJy
        mmod[ibad] = -10.0
        ibad = np.where((mphys <= 0) | (mphys > 45))
        mphys = -0.4 * (mphys - 23.91)  # uJy
        mphys[ibad] = -10.0
        ibad = np.where((mfir <= 0) | (mfir > 45))
        mfir = -0.4 * (mfir - 23.91)  # uJy
        mfir[ibad] = -10.0

        # Read the P(z) as an array
        # First column is the redshift, then the Bayesian posterior, then the profile likelihood
        zpdf = np.zeros([3, npdf])
        for i in range(npdf):
            line = fsp.readline()
            zpdf[:, i] = np.array(line.split())

        # Read the best-fit template [lambda(A), Mag(AB)]
        # convert in log(F_nu) = -0.4*(mab-23.91) [uJy]
        # Loop over the 6 models (GAL-1, GAL-2, GAL-FIR, GAL-STOCH, QSO, STAR)
        lg = []
        mg = []
        for m in range(6):
            nline = int(models_info[m][1])  # Get the number of lines of the templates given above
            bid = np.zeros([2, nline])
            if nline > 0:
                for i in range(nline):
                    line = fsp.readline()
                    bid[:, i] = np.array(line.split())
                    if bid[1, i] > 35:
                        bid[1, i] = -10.0
                    else:
                        bid[1, i] = -0.4 * (bid[1, i] - 23.91)  # log(flux)
            lg.append(bid[0, :] / 10000.0)  # From anstrom to micron
            mg.append(bid[1, :])  # log(flux)

    # fsp.close()

    ##############  PLOT THE MAG AND SED VERSUS LAMBDA ###############

    ### Initialise the figure
    fig = plt.figure()

    ### Main panel to plot log(flux) versus lambda, and the associated best-fit templates
    ax1 = fig.add_axes(
        [0.1, 0.1, 0.78, 0.78],
        xscale="log",
        xlabel=r"$\lambda$ [$\mu$m]",
        ylabel="log(F$_{\\nu}$) [$\\mu$Jy]",
    )

    # only the reliable obs mag (0<mag<35) will be plotted:
    em = em * 2.0
    dlf = dlf / 2.0
    mag1 = mag[(mag > 0.0) & (mag < 35) & (em > -3) & (dlf > 50)]
    em1 = em[(mag > 0.0) & (mag < 35) & (em > -3) & (dlf > 50)]
    lf1 = lf[(mag > 0.0) & (mag < 35) & (em > -3) & (dlf > 50)] / 10000.0  # in micron
    dlf1 = dlf[(mag > 0.0) & (mag < 35) & (em > -3) & (dlf > 50)] / 10000.0  # in micron

    # Define the limits of the axis, for magnitude having error<1 (if none, try 10)
    if len(mag1[em1 < 1] > 0):
        ymin = max(mag1[em1 < 1] + 3.0)  # 3 mag added to the faintest magnitude
        ymax = min(mag1[em1 < 1] - 2.0)  # 2 mag subsracted to the brightest magnitude
    else:
        # if no magnitude having error<1, try 10
        if len(mag1[em1 < 10] > 0):
            ymin = max(mag1[em1 < 10] + 3.0)
            ymax = min(mag1[em1 < 10] - 2.0)
        else:
            ymin = 15
            ymax = 35
    if ymin > 60:
        ymin = 30

    # high S/N bands:
    ic = (em1 >= 0.0) & (em1 < 2.0)
    lf2 = lf1[ic]  # effective wavelength in micron of each filter
    dlf2 = dlf1[ic]  # filter width
    mag2 = -0.4 * (mag1[ic] - 23.91)  # mag -> log(flux)
    em2 = 0.4 * em1[ic]  # err_mag -> err of log(flux)
    # low S/N bands:
    ic2 = (em1 >= 2.0) & (em1 < 8.0)
    lf2b = lf1[ic2]
    mag2b = -0.4 * (mag1[ic2] - 23.91)
    em2b = 0.4 * em1[ic2]
    dlf2b = dlf1[ic2]

    # set the main plot aspect
    if len(lf1) > 0:
        ax1.axis([min(lf1) * 0.85, max(lf1) * 2, -0.4 * (ymin - 23.91), -0.4 * (ymax - 23.91)])
    else:
        ax1.axis([0, 100000, -0.4 * (ymin - 23.91), -0.4 * (ymax - 23.91)])
    ### plot SED and print info of best-fit models
    # The different colors correspond to the 6 classes:
    # GAL-1, GAL-2, GAL-FIR, GAL-STOCH, QSO, STAR
    col_lst = ["r", "g", "b", "m", "y", "c"]  # each one with a different color
    w_lst = [2, 0.7, 0.7, 1, 1, 0.7]  # different line weight
    a_lst = [1, 0.7, 0.7, 0.7, 0.7, 0.7]  # different alpha
    # Display the information for each of the 6 classes
    plt.figtext(
        0.10,
        0.98,
        " Type: (Model Library Nband) z_phot  Chi^2,  Extlaw  EB-V  Lir  Age  logM*  logSFR",
        size="small",
    )
    plt.figtext(0.13, 0.8, "ID=" + id)
    iml = 0
    # Loop over the 6 classes  GAL-1, GAL-2, GAL-FIR, GAL-STOCH, QSO, STAR
    for im in range(5, -1, -1):
        if int(models_info[im][2]) == -1:
            continue  # print only models really used
        iml = iml + 1  # counter of models used
        # plot the SED: lambda in micron versus log(flux)
        ax1.plot(lg[im], mg[im], color=col_lst[im], alpha=a_lst[im], linewidth=w_lst[im])
        del models_info[im][6:8]  # do not print z_inf and z_sup
        del models_info[im][-1]  # nor sSFR
        info1 = ("  ".join(["%.3f"] * len(models_info[im][5:7]))) % tuple(
            [float(j) for j in models_info[im][5:7]]
        )
        if float(models_info[im][8]) >= 0.0:  # additional information
            info2 = ("   ".join(["%.2f"] * len(models_info[im][8:]))) % tuple(
                [float(j) for j in models_info[im][8:]]
            )
            info2 = ",  " + info2 + "."
        else:
            info2 = "."
        infol = models_info[im][0] + ": (" + " ".join(models_info[im][2:5]) + ")  " + info1 + info2
        plt.figtext(0.15, 0.98 - 0.022 * iml, infol, color=col_lst[im], size="x-small")  # print the rest

    # plot the observed log(flux) ...
    ax1.scatter(lf / 10000.0, mmod, color="b", marker="s", s=30)  # modeled log(flux)
    ax1.errorbar(lf2b, mag2b, yerr=em2b, xerr=dlf2b, fmt="o", color="0.6")  # observed log(flux) / low S/N
    ax1.errorbar(lf2, mag2, yerr=em2, xerr=dlf2, fmt="o", color="0.")  # observed log(flux) / high S/N
    # ... and upper limits (define with a negative error in flux or mag)
    iu = np.where(em1 < 0)
    if len(iu[0]) > 0:
        lf3 = lf1[iu]
        mag3 = -0.4 * (mag1[iu] - 23.91)
        # Make an arrow when upper-limit
        ax1.quiver(lf3, mag3, 0, -1, units="height", width=0.004, headwidth=5, color="k", pivot="tip")

    ##############  PLOT THE PDF(z) ###############

    ### 2nd panel (inset) showing PDF(z)
    ax2 = fig.add_axes([0.55, 0.17, 0.3, 0.25])
    poscond = np.where((zpdf[1, :] / max(zpdf[1, :]) > 0.002) | (zpdf[2, :] / max(zpdf[2, :]) > 0.002))
    zaxe = zpdf[0, poscond]
    ax2.axis([min(zaxe[0]) - 0.1, max(zaxe[0] + 0.1), 0, 1.1])
    # Plot the PDF as the sum of all the probability at a given redshift (Bayesian)
    ax2.plot(zpdf[0, :], zpdf[1, :] / max(zpdf[1, :]), color="k", linestyle="solid", label="Bayesian")
    # Plot the PDF as the minimum chi2 converted into probability at every redshift
    ax2.plot(zpdf[0, :], zpdf[2, :] / max(zpdf[2, :]), color="k", linestyle="dashed", label="Profile")
    # Plot a vertical line with the spec-z value if it exists
    if zspec > 0:
        ax2.plot([zspec, zspec], [0, 1.1], color="r", linestyle="solid", label="spec z")
    ax2.plot([zphot, zphot], [0, 1.1], color="b", linestyle="solid", label="photo z")
    ax2.legend(fontsize="10")
    ax2.set_title("PDF(z)")
