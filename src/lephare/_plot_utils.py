# from math import *

from math import ceil

import matplotlib.pyplot as plt
import numpy as np


#
#
# Class to generate many diagnostics on the output table created by LePHARE
#
#
# Example of run
#
#   t = Table.read("outputphotoz.fits")
#   utils = lp.PlotUtils(t, sel_filt=3,pos_filt=[0,1,2,4,5,5],
#           range_z=[0,0.5,1,1.5,3],range_mag=[19,20.5,21.5,22.5,25])
#   utils.zml_zs()
#   utils.zp_zs()
#   utils.zml_zp()
#   utils.distz()
#   utils.chi2dist()
#   utils.dist_filt()
#   utils.dist_model()
#   utils.dist_ebv()
#   utils.secondpeak()
#   utils.bzk()
#   utils.absmag_z()
#   utils.rf_color()
#   utils.william()
#   utils.cumulative68()
#   utils.check_error()
#   utils.errormag()
#   utils.errorz()
#
class PlotUtils:
    # Initialisation of the variables
    # Reading the input table
    # Define the filter in which we do the plots
    #
    # Input
    #   t         : table with the output from the LePHARE run
    #   sel_filt   : index of the filter used to select in magnitude (starting at 0)
    #   asignFilt : array with the index of the filters u, g, r, z, J, Ks
    #               used to create some plot in colors (starting at 0)
    #   range_z   : range in redshift to be considered (by default, compute 4 bins using quantiles)
    #   range_mag : range in mag to be considered (by default, compute 4 bins using quantiles)
    #
    def __init__(self, t, sel_filt=0, pos_filt=None, range_z=None, range_mag=None):
        if pos_filt is None:
            pos_filt = [0, 0, 0, 0, 0, 0]
        if range_z is None:
            range_z = [-1]
        if range_mag is None:
            range_mag = [-1]
        # Number of the filter start at 0
        self.sel_filt = sel_filt  # filter for the selection in mag
        self.uFilt = pos_filt[0]
        self.bFilt = pos_filt[1]
        self.rFilt = pos_filt[2]
        self.zFilt = pos_filt[3]
        self.jFilt = pos_filt[4]
        self.kFilt = pos_filt[5]

        # Read the input file
        self.Id = t["IDENT"]
        self.zp = t["Z_BEST"]
        self.zl68 = t["Z_BEST68_LOW"]
        self.zu68 = t["Z_BEST68_HIGH"]
        self.zml = t["Z_MED"]
        self.zmll68 = t["Z_MED68_LOW"]
        self.zmlu68 = t["Z_MED68_HIGH"]
        self.chi = t["CHI_BEST"]
        self.mod = t["MOD_BEST"]
        self.law = t["EXTLAW_BEST"]
        self.ebv = t["EBV_BEST"]
        self.zp2 = t["Z_SEC"]
        self.chi2 = t["CHI_SEC"]
        self.mod2 = t["MOD_SEC"]
        self.ebv2 = t["EBV_SEC"]
        self.zq = t["ZQ_BEST"]
        self.chiq = t["CHI_QSO"]
        self.modq = t["MOD_QSO"]
        self.mods = t["MOD_STAR"]
        self.chis = t["CHI_STAR"]
        self.mag = t["MAG_OBS()"][:, self.sel_filt]
        self.magu = t["MAG_OBS()"][:, self.uFilt]
        self.magb = t["MAG_OBS()"][:, self.bFilt]
        self.magr = t["MAG_OBS()"][:, self.rFilt]
        self.magz = t["MAG_OBS()"][:, self.zFilt]
        self.magj = t["MAG_OBS()"][:, self.jFilt]
        self.magk = t["MAG_OBS()"][:, self.kFilt]
        self.mabsu = t["MAG_ABS()"][:, self.uFilt]
        self.mabsb = t["MAG_ABS()"][:, self.bFilt]
        self.mabsr = t["MAG_ABS()"][:, self.rFilt]
        self.mabsz = t["MAG_ABS()"][:, self.zFilt]
        self.mabsj = t["MAG_ABS()"][:, self.jFilt]
        self.mabsk = t["MAG_ABS()"][:, self.kFilt]
        self.scale = t["SCALE_BEST"]
        self.nbFilt = t["NBAND_USED"]
        self.context = t["CONTEXT"]
        self.zs = t["ZSPEC"]
        self.ageb = t["AGE_BEST"]
        self.agel = t["AGE_INF"]
        self.agem = t["AGE_MED"]
        self.ages = t["AGE_SUP"]
        self.ldustb = t["LDUST_BEST"]
        self.ldustl = t["LDUST_INF"]
        self.ldustm = t["LDUST_MED"]
        self.ldusts = t["LDUST_SUP"]
        self.ltirb = t["LUM_TIR_BEST"]
        self.ltirl = t["LUM_TIR_INF"]
        self.ltirm = t["LUM_TIR_MED"]
        self.ltirs = t["LUM_TIR_SUP"]
        self.massb = t["MASS_BEST"]
        self.massl = t["MASS_INF"]
        self.massm = t["MASS_MED"]
        self.masss = t["MASS_SUP"]
        self.sfrb = t["SFR_BEST"]
        self.sfrl = t["SFR_INF"]
        self.sfrm = t["SFR_MED"]
        self.sfrs = t["SFR_SUP"]
        self.ssfrb = t["SSFR_BEST"]
        self.ssfrl = t["SSFR_INF"]
        self.ssfrm = t["SSFR_MED"]
        self.ssfrs = t["SSFR_SUP"]
        self.Lnuv = t["LUM_NUV_BEST"]
        self.Lr = t["LUM_R_BEST"]
        self.Lk = t["LUM_K_BEST"]

        # Define the panels with the binning in redshift an magnitude
        if len(range_z) == 1:
            self.range_z = np.quantile(self.zs[(self.zs > -1) & (self.zs < 9)], [0, 0.25, 0.5, 0.75, 1])
        else:
            self.range_z = range_z
        if len(range_mag) == 1:
            self.range_mag = np.quantile(self.mag[(self.mag > 10) & (self.mag < 40)], [0, 0.25, 0.5, 0.75, 1])
        else:
            self.range_mag = range_mag
        # Define the maximums in redshift and magnitude
        self.z_min, self.z_max = np.amin(self.range_z), np.amax(self.range_z)
        self.mag_min, self.mag_max = np.amin(self.range_mag), np.amax(self.range_mag)

        # Define the number panels
        # case with several bin of mag
        if len(self.range_mag) > 1:
            # Number of row with 2 columns
            self.nbRowM = int(ceil(float(len(self.range_mag) - 1) / 2.0))
            self.nbColM = 2
        else:
            self.nbColM = 1
            self.nbRowM = 1
        # case with several bin of z
        if len(self.range_z) > 1:
            # Number of row with 2 columns
            self.nbRowZ = int(ceil(float(len(self.range_z) - 1) / 2.0))
            self.nbColZ = 2
        else:
            self.nbColZ = 1
            self.nbRowZ = 1

        # Condition for selection
        # general condition to select the galaxies in the expected z/mag range
        self.cond = (
            (self.zp > self.z_min)
            & (self.zp < self.z_max)
            & (self.zml > self.z_min)
            & (self.zml < self.z_max)
            & (self.mag > self.mag_min)
            & (self.mag < self.mag_max)
        )
        # condition to select stars
        self.condstar = self.chis < self.chi
        # condition to select galaxies
        self.condgal = ~self.condstar
        # condition to select spectroscopic redshifts
        self.condspec = (self.zs > 0) & (self.zs < 9)

        return

    # Plot photo-z (median PDF) versus spec-z
    def zml_zs(self):
        # Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "$z_{spec}$", ha="center", fontsize=15)
        f.text(0.04, 0.5, r"$z_{phot}\; median\; PDF(z)$", va="center", rotation="vertical", fontsize=15)

        # Loop over the magnitude bins
        for rm in range(len(self.range_mag)):
            # Define the subplots and pass the panel
            ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

            # mag limits
            magmin = self.range_mag[rm - 1]
            magmax = self.range_mag[rm]

            ax.axis([self.z_min, self.z_max, self.z_min, self.z_max])
            conda = self.cond & self.condspec & (self.mag > magmin) & (self.mag < magmax)
            ax.scatter(self.zs[conda], self.zml[conda], s=1, color="b", alpha=0.5, marker="s")

            # Check that we have some sources before performing statistics
            ngal = len(self.zml[conda])
            if ngal > 0:
                # statistics
                arg_bias = self.zml[conda] - self.zs[conda]
                arg_std = arg_bias / (1.0 + self.zs[conda])
                nmad = 1.4821 * np.median(abs(arg_std))
                cond_outl = abs(arg_std) > 0.15
                outl_rate = len(arg_std[cond_outl]) / float(ngal)
                ax.annotate(
                    r"$"
                    + str(round(magmin, 2))
                    + " < mag < "
                    + str(round(magmax, 2))
                    + "$ \n"
                    + "$N_{gal}  = "
                    + str(ngal)
                    + "$ \n"
                    + r"$\eta  ="
                    + str(round(100 * outl_rate, 2))
                    + "  \\%$\n"
                    + r"$ \sigma_{\Delta z /(1+z)}  = "
                    + str(round(nmad, 5))
                    + "$",
                    xy=(0.05 * self.z_max, 0.65 * self.z_max),
                    color="black",
                    fontsize=12,
                )

            # Trace the limits 0.15(1+z)
            x_zs = np.array([0, 6])
            ax.plot(x_zs, x_zs * 1.15 + 0.15, "c--")
            ax.plot(x_zs, x_zs, "r-")
            ax.plot(x_zs, x_zs * 0.85 - 0.15, "c--")

        return

    # Plot photo-z (minimum chi2) versus spec-z
    def zp_zs(self):
        # Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "$z_{spec}$", ha="center", fontsize=15)
        f.text(0.04, 0.5, r"$z_{phot}\; median\; PDF(z)$", va="center", rotation="vertical", fontsize=15)

        # Loop over the magnitude bins
        for rm in range(len(self.range_mag)):
            # Define the subplots and pass the panel
            ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

            # mag limits
            magmin = self.range_mag[rm - 1]
            magmax = self.range_mag[rm]

            ax.axis([self.z_min, self.z_max, self.z_min, self.z_max])
            conda = self.cond & self.condspec & (self.mag > magmin) & (self.mag < magmax)
            ax.scatter(self.zs[conda], self.zp[conda], s=1, color="b", alpha=0.5, marker="s")

            # Check that we have some sources before performing statistics
            ngal = len(self.zp[conda])
            if ngal > 0:
                # statistics
                arg_bias = self.zp[conda] - self.zs[conda]
                arg_std = arg_bias / (1.0 + self.zs[conda])
                nmad = 1.4821 * np.median(abs(arg_std))
                cond_outl = abs(arg_std) > 0.15
                outl_rate = len(arg_std[cond_outl]) / float(ngal)
                ax.annotate(
                    r"$"
                    + str(round(magmin, 2))
                    + " < mag < "
                    + str(round(magmax, 2))
                    + "$ \n"
                    + "$N_{gal}  = "
                    + str(ngal)
                    + "$ \n"
                    + r"$\eta  ="
                    + str(round(100 * outl_rate, 2))
                    + "  \\%$\n"
                    + r"$ \sigma_{\Delta z /(1+z)}  = "
                    + str(round(nmad, 5))
                    + "$",
                    xy=(0.05 * self.z_max, 0.65 * self.z_max),
                    color="black",
                    fontsize=12,
                )

            # Trace the limits 0.15(1+z)
            x_zs = np.array([0, 6])
            ax.plot(x_zs, x_zs * 1.15 + 0.15, "c--")
            ax.plot(x_zs, x_zs, "r-")
            ax.plot(x_zs, x_zs * 0.85 - 0.15, "c--")

        return

    # Plot photo-z (median of the PDF) versus photo-z (minimum chi2)
    def zml_zp(self):
        # Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, r"$z_{phot}\; minimum\; \chi^2$", ha="center", fontsize=15)
        f.text(0.04, 0.5, r"$z_{phot}\; median\; PDF(z)$", va="center", rotation="vertical", fontsize=15)

        # Loop over the magnitude bins
        for rm in range(len(self.range_mag)):
            if rm > 0:
                # Define the subplots and pass the panel
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # mag limits
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]

                # set the axis
                ax.axis([self.z_min, self.z_max, self.z_min, self.z_max])

                # new condition with the magnitude range
                conda = self.cond & self.condgal & (self.mag > magmin) & (self.mag < magmax)

                # Plot photo-z versus spec-z
                ax.scatter(self.zp[conda], self.zml[conda], s=1, color="b", alpha=0.5, marker="s")

                # Trace the limits 0.15(1+z)
                x_zs = np.array([0, 6])
                ax.plot(x_zs, x_zs * 1.15 + 0.15, "c--")
                ax.plot(x_zs, x_zs, "r-")
                ax.plot(x_zs, x_zs * 0.85 - 0.15, "c--")

                # labels
                ax.annotate(
                    "$" + str(round(magmin, 2)) + " < mag < " + str(round(magmax, 2)) + "$",
                    xy=(0.1 * self.z_max, 0.8 * self.z_max),
                    color="black",
                    fontsize=15,
                )

        return

    # Redshift distribution of the spec-z sample, zml and zbest
    def distz(self):
        # redshift step for the histogram
        nstep = 20

        ## Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "$Redshift$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$N_{normalized}$", va="center", rotation="vertical", fontsize=15)

        # Loop over the magnitude bins
        for rm in range(len(self.range_mag)):
            # No legend
            leg = 0

            if rm > 0:
                if rm == 1:
                    leg = 1  # legend in the first panel

                # Define the subplots and pass the panel
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # mag limits
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]

                # set the axis
                ax.axis([0, self.z_max, 0, 2.9])

                # new condition with the magnitude range
                conda = self.cond & self.condgal & (self.mag > magmin) & (self.mag < magmax) & (self.zml > 0)
                # Check that some object are present
                if len(self.zp[conda]) > 0:
                    # Histogram with the photometric redshifts median PDF
                    ax.hist(
                        self.zml[conda],
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="b",
                        label=r"$z_{phot}\; median\; PDF(z)$",
                    )

                    # Histogram with the photometric redshifts minimum chi2
                    ax.hist(
                        self.zp[conda],
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="r",
                        alpha=0.5,
                        label=r"$z_{phot}\; minimum\; \chi^2$",
                    )

                # new condition with the magnitude range and spectro-z
                conda = self.cond & self.condspec & (self.mag > magmin) & (self.mag < magmax)
                # Check that some object are present
                if len(self.zp[conda]) > 0:
                    # Histogram with the photometric redshifts zbest
                    ax.hist(
                        self.zs[conda],
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="g",
                        alpha=0.2,
                        label="$z_{spectro}$",
                    )
                    # labels
                    ax.annotate(
                        "$" + str(round(magmin, 2)) + " < mag < " + str(round(magmax, 2)) + "$",
                        xy=(0.5 * self.z_max, 1.5),
                        color="black",
                        fontsize=14,
                    )

                # print the legend
                if leg == 1:
                    ax.legend(prop={"size": 8})

        return

    # Chi2 distribution per magnitude bin
    def chi2dist(self):
        # nb ebv step for the histogram
        nstep = 20

        ## Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, r"log($\chi^2$)", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$N_{normalized}$", va="center", rotation="vertical", fontsize=15)

        # Loop over the magnitude bins
        for rm in range(len(self.range_mag)):
            # No legend
            leg = 0

            if rm > 0:
                if rm == 2:
                    leg = 1  # legend in the first panel

                # Define the subplots and pass the panel
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # mag limits
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]

                # set the axis
                ax.axis([-2, 2.9, 0, 1.49])

                # new condition with the magnitude range for galaxies, chi gal< chi star
                conda = (
                    self.cond
                    & self.condgal
                    & (self.mag > magmin)
                    & (self.mag < magmax)
                    & (self.nbFilt > 2)
                    & (self.chi > 0)
                )
                if len(self.zp[conda]) > 1:
                    chireduit = self.chi[conda] / (self.nbFilt[conda] - 1)
                    # Histogram of log (reduced chi2 gal)
                    ax.hist(
                        np.log10(chireduit),
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="b",
                        label=r"$\chi^2_{gal} \; if \; \chi^2_{gal}<\chi^2_{stars}$",
                    )

                # new condition with the magnitude range for stars, chi star < chi gal
                conda = (
                    self.cond
                    & self.condstar
                    & (self.mag > magmin)
                    & (self.mag < magmax)
                    & (self.nbFilt > 2)
                    & (self.chi > 0)
                )
                if len(self.zp[conda]) > 1:
                    chireduit = self.chis[conda] / (self.nbFilt[conda] - 1)
                    # Histogram of log (reduced chi2 gal)
                    ax.hist(
                        np.log10(chireduit),
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="r",
                        alpha=0.5,
                        label=r"$\chi^2_{star} \; if \;  \chi^2_{gal}>\chi^2_{stars}$",
                    )

                    # print the legend
                    if leg == 1:
                        ax.legend(prop={"size": 8})

                    # labels
                    ax.annotate(
                        "$" + str(round(magmin, 2)) + " < mag < " + str(round(magmax, 2)) + "$",
                        xy=(-1.8, 1.2),
                        color="black",
                        fontsize=15,
                    )

        return

    # Number of filters used in the fit
    def dist_filt(self):
        # redshift step for the histogram
        nstep = 10

        ## Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowZ, self.nbColZ, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "Number of filters", ha="center", fontsize=15)
        f.text(0.04, 0.5, "N", va="center", rotation="vertical", fontsize=15)

        # Loop over the redshift bins
        for rm in range(len(self.range_z)):
            if rm > 0:
                # Define the subplots and pass the panel
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # mag limits
                zmin = self.range_z[rm - 1]
                zmax = self.range_z[rm]

                # set the axis
                ax.axis([0, 50, 0, 3])

                # new condition with the magnitude range
                conda = self.cond & (self.zp > zmin) & (self.zp < zmax)

                # Histogram with the number of filters
                ax.hist(self.nbFilt[conda], bins=nstep, histtype="stepfilled", density=1, color="r")

                # labels
                ax.annotate(
                    "$" + str(round(zmin, 2)) + " < z < " + str(round(zmax, 2)) + "$",
                    xy=(0.1, 1),
                    color="black",
                    fontsize=15,
                )

        return

    # Template distribution
    def dist_model(self):
        # redshift step for the histogram
        nstep = 20

        ## Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowZ, self.nbColZ, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "best fit template", ha="center", fontsize=15)
        f.text(0.04, 0.5, "N", va="center", rotation="vertical", fontsize=15)

        # Loop over the redshift bins
        for rm in range(len(self.range_z)):
            if rm > 0:
                # Define the subplots and pass the panel
                zmin = self.range_z[rm - 1]
                zmax = self.range_z[rm]
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # set the axis
                ax.axis([0, 70, 0, 0.3])

                # new condition with the magnitude range
                conda = self.cond & self.condgal & (self.zp > zmin) & (self.zp < zmax)

                # check if some objects exist
                if len(self.mod[conda]) > 1:
                    # Histogram of the models
                    ax.hist(self.mod[conda], bins=nstep, histtype="stepfilled", density=1, color="b")

                    # labels
                    ax.annotate(
                        "$" + str(round(zmin, 2)) + " < z < " + str(round(zmax, 2)) + "$",
                        xy=(1, 0.15),
                        color="black",
                        fontsize=15,
                    )

        return

    # EBV distribution
    def dist_ebv(self):
        # number of steps for the histogram
        nstep = 10

        ## Create a figure with an array with nbRowZ*nbColZ subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowZ, self.nbColZ, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "E(B-V)", ha="center", fontsize=15)
        f.text(0.04, 0.5, "N", va="center", rotation="vertical", fontsize=15)

        # Loop over the redshift bins
        for rm in range(len(self.range_z)):
            if rm > 0:
                # Define the subplots and pass the panel
                zmin = self.range_z[rm - 1]
                zmax = self.range_z[rm]
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # set the axis
                ax.axis([0, 0.6, 0, 20])

                # new condition with the magnitude range
                conda = self.cond & self.condgal & (self.zp > zmin) & (self.zp < zmax)

                # check if some objects exist
                if len(self.ebv[conda]) > 1:
                    # Histogram with the E(B-V)
                    ax.hist(self.ebv[conda], bins=nstep, histtype="stepfilled", density=1, color="b")

                    # labels
                    ax.annotate(
                        "$" + str(round(zmin, 2)) + " < z < " + str(round(zmax, 2)) + "$",
                        xy=(0.02, 15),
                        color="black",
                        fontsize=15,
                    )

        return

        return

    # Redshift distribution of the second peak photo-z, zml and zbest of these sources
    def secondpeak(self):
        # redshift step for the histogram
        nstep = 20

        ## Create a figure with an array with nbRowM*nbColM subpanels
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        # No space between the figures
        plt.subplots_adjust(hspace=0, wspace=0)

        # label of the figure
        f.text(0.5, 0.04, "$Redshift$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$N_{normalized}$", va="center", rotation="vertical", fontsize=15)

        # Loop over the magnitude bins
        for rm in range(len(self.range_mag)):
            # No legend
            leg = 0

            if rm > 0:
                if rm == 4:
                    leg = 1  # legend in the first panel

                # Define the subplots and pass the panel
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]

                # mag limits
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]

                # set the axis
                ax.axis([0, self.z_max, 0, 5])

                # new condition with the magnitude range
                conda = self.cond & self.condgal & (self.mag > magmin) & (self.mag < magmax) & (self.zp2 > 0)

                # If some sources
                if len(self.zp2[conda]) > 1:
                    # Histogram with the photometric redshifts median PDF
                    ax.hist(
                        self.zp2[conda],
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="b",
                        label=r"$z_{phot}\; second\; peak$",
                    )

                    # Histogram with the photometric redshifts minimum chi2
                    ax.hist(
                        self.zp[conda],
                        bins=nstep,
                        histtype="stepfilled",
                        density=1,
                        color="r",
                        alpha=0.5,
                        label=r"$z_{phot}\; minimum \; \chi^2$",
                    )

                    # labels
                    ax.annotate(
                        "$" + str(round(magmin, 2)) + " < mag < " + str(round(magmax, 2)) + "$",
                        xy=(0.2, 3.5),
                        color="black",
                        fontsize=13,
                    )

                # print the legend
                if leg == 1:
                    ax.legend(prop={"size": 8})

        return

    # BzK per redshift bin
    def bzk(self):
        plt.clf()
        f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$B-z$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$z-K$", va="center", rotation="vertical", fontsize=15)

        bzk_conf = [[0.001, 1.4, 0, 0], [1.4, 2.7, 0, 1], [2.7, 6, 1, 0], [-0.001, 6, 1, 1]]

        for rm in range(4):
            zmin = bzk_conf[rm][0]
            zmax = bzk_conf[rm][1]
            ax = axarr[bzk_conf[rm][2], bzk_conf[rm][3]]
            ax.axis([-0.5, 5.5, -1.0, 4.5])

            if rm < 3:
                conda = self.cond & self.condgal & (self.zp > zmin) & (self.zp < zmax)
            else:
                conda = self.cond & self.condstar

            if (self.bFilt >= 0) and (self.zFilt >= 0) and (self.kFilt >= 0):
                colx = self.magb - self.magz
                coly = self.magz - self.magk
                ax.scatter(colx[conda], coly[conda], s=1, color="b", alpha=0.2, marker="s")

            ax.plot([-1, 6], [-0.5, 6])

            if rm < 3:
                ax.annotate(f"${zmin:.2f} < z < {zmax:.2f}$", xy=(0.1, 3), color="black", fontsize=15)
            else:
                ax.annotate(r"$\chi_s < \chi_g$", xy=(0.1, 3), color="black", fontsize=15)

        return

    # Absolute magnitude versus redshift
    def absmag_z(self):
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$Redshift$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$M_B$", va="center", rotation="vertical", fontsize=15)

        for rm in range(len(self.range_mag)):
            if rm > 0:
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]

                ax.axis([self.z_min, self.z_max, -16, -25.9])

                conda = self.cond & self.condgal & (self.mag > magmin) & (self.mag < magmax)
                ax.scatter(self.zp[conda], self.mabsb[conda], s=1, color="b", alpha=0.2, marker="s")

                ax.annotate(f"${magmin:.2f} < mag < {magmax:.2f}$", xy=(0.1, -25), color="black", fontsize=15)

    # Rest-frame colors versus absolute magnitude
    def rf_color(self):
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$M_R$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$M_U-M_R$", va="center", rotation="vertical", fontsize=15)

        for rm in range(len(self.range_z)):
            if rm > 0:
                zmin = self.range_z[rm - 1]
                zmax = self.range_z[rm]
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]
                ax.axis([-24.9, -16, -0.7, 2.4])

                conda = self.cond & self.condgal & (self.zp > zmin) & (self.zp < zmax)

                if self.uFilt >= 0 and self.rFilt >= 0:
                    magx = self.mabsr
                    coly = self.mabsu - self.mabsr

                    ax.scatter(magx[conda], coly[conda], s=1, color="b", alpha=0.2, marker="s")

                    ax.annotate(f"${zmin:.2f} < z < {zmax:.2f}$", xy=(-20, 1), color="black", fontsize=15)

        return

    # Williamrest-frame color-color plot
    def william(self):
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$M_R-M_J$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$M_U-M_R$", va="center", rotation="vertical", fontsize=15)

        for rm in range(len(self.range_z)):
            if rm > 0:
                zmin = self.range_z[rm - 1]
                zmax = self.range_z[rm]
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]
                ax.axis([-2.1, 1.9, 0, 2.5])

                conda = self.cond & self.condgal & (self.zp > zmin) & (self.zp < zmax)

                if self.uFilt >= 0 and self.rFilt >= 0 and self.jFilt >= 0:
                    colx = self.mabsr - self.mabsj
                    coly = self.mabsu - self.mabsr

                    ax.scatter(colx[conda], coly[conda], s=1, color="b", alpha=0.2, marker="s")

                    ax.annotate(f"${zmin:.2f} < z < {zmax:.2f}$", xy=(0.1, 2), color="black", fontsize=15)

        return

    # Check that 68% of the spec-z fall in the 68% error
    def cumulative68(self):
        nstep = 40
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$abs((z_{spec}-z_{phot})/z_{phot} uncertainties)$", ha="center", fontsize=15)
        f.text(0.04, 0.5, "$cumulative(N)$", va="center", rotation="vertical", fontsize=15)

        for rm in range(len(self.range_mag)):
            if rm > 0:
                ax = axarr[int(ceil(rm / 2.0) - 1), int(ceil(rm % 2) - 1)]
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]
                ax.axis([0, 9.9, 0, 0.99])

                diffzml = abs(self.zml - self.zs) / np.maximum(
                    abs(self.zmlu68 - self.zml), abs(self.zmll68 - self.zml)
                )
                diffzp = abs(self.zp - self.zs) / np.maximum(
                    abs(self.zu68 - self.zp), abs(self.zl68 - self.zp)
                )

                condazml = (
                    self.cond
                    & self.condspec
                    & (self.mag > magmin)
                    & (self.mag < magmax)
                    & (diffzml < 10)
                    & (self.zml > 0)
                )
                condazp = (
                    self.cond
                    & self.condspec
                    & (self.mag > magmin)
                    & (self.mag < magmax)
                    & (diffzp < 10)
                    & (self.zp > 0)
                )

                if len(self.zml[condazml]) > 1 and len(self.zp[condazp]) > 1:
                    ax.hist(
                        diffzp[condazp],
                        bins=nstep,
                        histtype="step",
                        density=1,
                        color="b",
                        cumulative=True,
                        label=r"$z_{phot}\; minimum \; \chi^2$",
                    )
                    ax.hist(
                        diffzml[condazml],
                        bins=nstep,
                        histtype="step",
                        density=1,
                        color="r",
                        cumulative=True,
                        label=r"$z_{phot}\; median \; PDF$",
                    )

                    ax.axvline(x=1, color="r", linestyle="--")
                    ax.axhline(y=0.68, color="r", linestyle="--")

                    ax.legend(prop={"size": 8})

                ax.annotate(
                    f"${magmin:.2f} < mag < {magmax:.2f}$", xy=(0.1, 0.80), color="black", fontsize=15
                )

        return

    # Check that 68% of the spec-z fall in the 68% errorPlot error in photo-z as a function of magnitude
    def check_error(self):
        # Clear the figure
        plt.clf()
        # Set the axis
        plt.axis([15, 27, -1, 1])

        # 1 sigma error on the photo-z
        diffu = self.zu68 - self.zp
        diffl = self.zl68 - self.zp

        conda = self.cond & self.condgal

        # Error versus mag
        plt.scatter(self.mag[conda], diffu[conda], s=1, color="b", alpha=0.5, marker="s")
        plt.scatter(self.mag[conda], diffl[conda], s=1, color="r", alpha=0.5, marker="s")

        # 68%
        plt.plot([15, 27], [0, 0], color="r")  # Adapte les limites x à ton graphique

        plt.xlabel("magnitude", fontsize=18, labelpad=13)
        plt.ylabel(r"$68\% \; photo-z \; error$", fontsize=18, labelpad=13)

        return

    # Photo-z errors versus mag
    def errormag(self):
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$mag$", ha="center", fontsize=15)
        f.text(0.04, 0.5, r"$z_{phot} \; uncertainties$", va="center", rotation="vertical", fontsize=15)

        for rm in range(len(self.range_z)):
            if rm > 0:
                zmin = self.range_z[rm - 1]
                zmax = self.range_z[rm]
                ax = axarr[(rm // 2) % self.nbRowM, rm % 2]

                mstep = 0.2
                mstart = 16
                mnbstep = 100

                medpvec = []
                mednvec = []
                magvec = []

                ax.axis([15, 27, -0.99, 1])

                for im in range(mnbstep):
                    diffp = self.zmlu68 - self.zml
                    diffn = self.zmll68 - self.zml

                    mlimi = mstart + mstep * im
                    mlims = mstart + mstep * (im + 1)

                    conda = (
                        self.cond
                        & self.condgal
                        & (self.mag > mlimi)
                        & (self.mag <= mlims)
                        & (self.zml > zmin)
                        & (self.zml < zmax)
                    )

                    if len(self.zml[conda]) > 3:
                        medp = np.median(diffp[conda])
                        medn = np.median(diffn[conda])

                        medpvec.append(medp)
                        mednvec.append(medn)
                        magvec.append(mstart + mstep * (im + 0.5))

                ax.plot(magvec, medpvec, color="b", linestyle="-")
                ax.plot(magvec, mednvec, color="r", linestyle="-")

                ax.annotate(f"${zmin:.2f} < z < {zmax:.2f}$", xy=(mstart, 0.60), color="black", fontsize=15)

        return

    # Photo-z errors versus z
    def errorz(self):
        plt.clf()
        f, axarr = plt.subplots(self.nbRowM, self.nbColM, sharex=True, sharey=True, figsize=(12, 8))
        plt.subplots_adjust(hspace=0, wspace=0)

        f.text(0.5, 0.04, "$z_{phot}$", ha="center", fontsize=15)
        f.text(0.04, 0.5, r"$z_{phot} \; uncertainties$", va="center", rotation="vertical", fontsize=15)

        for rm in range(len(self.range_mag)):
            if rm > 0:
                ax = axarr[(rm // 2) % self.nbRowM, rm % 2]
                magmin = self.range_mag[rm - 1]
                magmax = self.range_mag[rm]

                zstep = 0.1
                zstart = 0
                znbstep = 60

                medpvec = []
                mednvec = []
                zvec = []

                ax.axis([0, 5.9, -0.99, 1])

                for iz in range(znbstep):
                    diffp = self.zmlu68 - self.zml
                    diffn = self.zmll68 - self.zml

                    zlimi = zstart + zstep * iz
                    zlims = zstart + zstep * (iz + 1)

                    conda = (
                        self.cond
                        & self.condgal
                        & (self.zml > zlimi)
                        & (self.zml <= zlims)
                        & (self.mag > magmin)
                        & (self.mag < magmax)
                    )

                    if len(self.zml[conda]) > 3:
                        medp = np.median(diffp[conda])
                        medn = np.median(diffn[conda])

                        medpvec.append(medp)
                        mednvec.append(medn)
                        zvec.append(zstart + zstep * (iz + 0.5))

                ax.plot(zvec, medpvec, color="b", linestyle="-")
                ax.plot(zvec, mednvec, color="r", linestyle="-")

                ax.annotate(
                    f"${magmin:.2f} < mag < {magmax:.2f}$",
                    xy=(zstart + 0.2, 0.60),
                    color="black",
                    fontsize=15,
                )

        return
