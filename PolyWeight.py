#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Created on Tue Jul  5 15:32:27 2022

@author: Atilio Minotto


"""

from tkinter.filedialog import asksaveasfilename
from itertools import islice
import time
from PIL import ImageTk, Image
from scipy.integrate import quad
import scipy.interpolate
from idlelib.tooltip import Hovertip
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
from tkinter.messagebox import askokcancel, showerror, showinfo, showwarning
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import tkinter as tk
import numpy as np
import os

script_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(script_dir)

from help_frame import *
from about_frame import *
from scipy.optimize import curve_fit

from scipy.integrate import quad
from scipy.special import gamma, gammaincc

import scipy.special as sc

from lmfit import  Model, Minimizer, Parameters, fit_report

from uncertainties import ufloat,correlated_values, unumpy
from uncertainties.umath import *

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                            DECLARATION OF GLOBAL VARIABLES: GAR
                            
################################################################################################"""

global delay_value
delay_value = 500

global data_flag, gpc_flag, mwd_flag
global data, tau, H, new_H, errorH, path_spec
global data_gpc, mass, distrib, path_gpc
global a, B, k, M0, G0, Me, opt
global Mn, Mw, Mz, MwMn
global M, mwd
global opt1
global clip_text


global un
un = " g/mol"

global optRouse, limRouse, slopeRouse, rouseMult

optRouse = 0
limRouse = 22
slopeRouse = 0
rouseMult = 0

global H2
global rouseAdj, popt, new_tau

plot_size_x = 6.1
plot_size_y = 4
plot_pos_x = 305
plot_pos_y = 485

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                            DECLARATION OF GLOBAL VARIABLES: GEX
                            
################################################################################################"""

global data_mod, omega, G1, G2, path_spec, mod_flag
global a_GEX, B_GEX, k_GEX, M0_GEX, G0_GEX, Me_GEX
global Mn_GEX, Mw_GEX, Mz_GEX, MwMn_GEX
global data_gpc_GEX, mass_GEX, distrib_GEX, path_gpc_GEX, gpc_flag_GEX, mwd_flag_GEX
global M_GEX, mwd_GEX
global opt2
global freqLim, weight, init_a, init_b, init_m0
global new_omega, new_G1, new_G2, comboY
global a_par, b_par, m0_par
global rTime
global G1_GEX_mod, G2_GEX_mod, omega_GEX_mod
global report
global fitting_report

freqLim = 100
weight = 1.0
init_a = 1.0
init_b = 1.0
init_m0 = 1.0

"""------------------------------------------------------------------------------------------------
###################################################################################################

                    CREATING AND CONFIGURING THE MAIN WINDOW: root
                    
################################################################################################"""

root = tk.Tk()
root.title('PolyWeight v0.5')

window_width = 940
window_height = 625

# get the screen dimension
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# find the center point
center_x = int(screen_width/2 - window_width / 2)
center_y = int(screen_height/2 - window_height / 2)

# set the position of the window to the center of the screen
root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')

#root.geometry("%dx%d" % (screen_width-2, screen_height-10))

root.resizable(False, False)

#root.attributes('-zoomed', True)

open_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-opened-folder-32.png"))
saveAs_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-save-as-32.png"))
#save_icon = ImageTk.PhotoImage(Image.open(script_dir+"/icons/icons8-save-32.png"))
material_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-molecule-32.png"))
spectrum_icon = ImageTk.PhotoImage(
    Image.open(script_dir+"/icons/icons8-sine-32.png"))
gpc_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-normal-distribution-histogram-32.png"))
help_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-help-32.png"))
info_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-info-32.png"))

run_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-circled-play-32.png"))
clear_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-broom-32.png"))
export_icon = ImageTk.PhotoImage(Image.open(
    script_dir+"/icons/icons8-export-32.png"))
console_icon = tk.PhotoImage(file=script_dir+"/icons/icons8-open-file-under-cursor-32.png")

main_icon = tk.PhotoImage(file=script_dir+"/icons/icons8-polymer-64.png")

root.iconphoto(False, main_icon)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

            CLEAR_ALL: CLEAR ALL PROGRAM VARIABLES AND GRAPHICS (ON GEX TAB)
            
################################################################################################"""


def clear_all_GEX():
    
    global data_mod, omega, G1, G2, path_spec, mod_flag
    global a_GEX, B_GEX, k_GEX, M0_GEX, G0_GEX, Me_GEX
    global Mn_GEX, Mw_GEX, Mz_GEX, MwMn_GEX
    global data_gpc_GEX, mass_GEX, distrib_GEX, path_gpc_GEX, gpc_flag_GEX, mwd_flag_GEX
    global M_GEX, mwd_GEX
    global opt2
    global new_omega, new_G1, new_G2, comboY
    global report, fitting_report
    global rTime
    global G1_GEX_mod, G2_GEX_mod, omega_GEX_mod
    
    data_mod, omega, G1, G2, path_spec, mod_flag = [], [], [], [], [], []
    a_GEX, B_GEX, k_GEX, M0_GEX, G0_GEX, Me_GEX = [], [], [], [], [], []
    Mn_GEX, Mw_GEX, Mz_GEX, MwMn_GEX = [], [], [], []
    data_gpc_GEX, mass_GEX, distrib_GEX, path_gpc_GEX, gpc_flag_GEX = [], [], [], [], []
    M_GEX, mwd_GEX = [], []
    opt2 = []
    new_omega, new_G1, new_G2, comboY = [], [], [], []
    report, fitting_report = [], []
    rTime = []
    G1_GEX_mod, G2_GEX_mod, omega_GEX_mod = [], [], []
    mwd_flag_GEX = 0
    
    selectMat2.set("Select a material")
    Hovertip(selectMat2, 'Select pre-saved material and reptation parameters',
             hover_delay=delay_value)

    rept_a_GEX.delete(0, 'end')
    rept_B_GEX.delete(0, 'end')
    mat_k_GEX.delete(0, 'end')
    mat_M0_GEX.delete(0, 'end')
    mat_G0_GEX.delete(0, 'end')
    mat_Me_GEX.delete(0, 'end')

    resMn_GEX.config(state='active')
    resMw_GEX.config(state='active')
    resMz_GEX.config(state='active')
    resMwMn_GEX.config(state='active')

    resMn_GEX.delete(0, 'end')
    resMw_GEX.delete(0, 'end')
    resMz_GEX.delete(0, 'end')
    resMwMn_GEX.delete(0, 'end')

    resMn_GEX.config(state='readonly')
    resMw_GEX.config(state='readonly')
    resMz_GEX.config(state='readonly')
    resMwMn_GEX.config(state='readonly')

    plot_clear2()

    menu4.entryconfig('Plot estimated distribution', state='disabled')
    menu4.entryconfig(
        'Plot MWD + estimated distribution', state='disabled')
    menu4.entryconfig('Plot dynamic moduli', state='disabled')
    menu4.entryconfig('Plot MWD data only', state='disabled')
    menu4.entryconfig('Plot dynamic moduli + fitted moduli', state='disabled')

    menu5.entryconfig('Copy to clipboard', state='disabled')
    
"""------------------------------------------------------------------------------------------------ 
###################################################################################################

            CLEAR_ALL: CLEAR ALL PROGRAM VARIABLES AND GRAPHICS (ON GAR TAB)
            
################################################################################################"""


def clear_all_GAR():
    global data_flag, gpc_flag, mwd_flag
    global data, tau, H, new_H, errorH, path_spec
    global data_gpc, mass, distrib, path_gpc
    global Mn, Mw, Mz, MwMn
    global a, B, k, M0, G0, Me, opt
    global M, mwd

    data_flag, gpc_flag = [], []
    a, B, k, M0, G0, Me, opt = [], [], [], [], [], [], []
    data, tau, H, new_H, errorH, path_spec = [], [], [], [], [], []
    data_gpc, mass, distrib, path_gpc = [], [], [], []
    Mn, Mw, Mz, MwMn = [], [], [], []
    M, mwd = [], []
    
    mwd_flag = 0
    
    selectMat1.set("Select a material")
    Hovertip(selectMat1, 'Select pre-saved material and reptation parameters',
             hover_delay=delay_value)

    rept_a_GAR.delete(0, 'end')
    rept_B_GAR.delete(0, 'end')
    mat_k_GAR.delete(0, 'end')
    mat_M0_GAR.delete(0, 'end')
    mat_G0_GAR.delete(0, 'end')
    mat_Me_GAR.delete(0, 'end')

    resMn.config(state='active')
    resMw.config(state='active')
    resMz.config(state='active')
    resMwMn.config(state='active')

    resMn.delete(0, 'end')
    resMw.delete(0, 'end')
    resMz.delete(0, 'end')
    resMwMn.delete(0, 'end')

    resMn.config(state='readonly')
    resMw.config(state='readonly')
    resMz.config(state='readonly')
    resMwMn.config(state='readonly')

    plot_clear1()

    menu1.entryconfig('Plot estimated distribution', state='disabled')
    menu1.entryconfig(
        'Plot MWD + estimated distribution', state='disabled')
    menu1.entryconfig('Plot relaxation time spectrum', state='disabled')
    menu1.entryconfig('Plot MWD data only', state='disabled')

    menu2.entryconfig('Copy to clipboard', state='disabled')


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

    CONFIRM_CLEAR: REQUESTS CONFIRMATION TO CLEAR ALL VARIABLES AND GRAPHICS
    
################################################################################################"""


def confirm_clear_GAR():
    answer = askokcancel(title='Clear all the data?',
                         message='All the data relative to the current analysis will be deleted.\
                         Do you wish to proceed?')
    if answer:
        clear_all_GAR()
        
def confirm_clear_GEX():
    answer = askokcancel(title='Clear all the data?',
                         message='All the data relative to the current analysis will be deleted.\
                         Do you wish to proceed?')
    if answer:
        clear_all_GEX()


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                    OPEN_SPECTRUM: OPENS RELAXATION SPECTRUM
                    
################################################################################################"""


def open_spectrum():

    global data, tau, H, errorH, path_spec, data_flag, mwd_flag

    path_spec = askopenfilename(
        filetypes=[("Text files", "*.txt *.dat *.sol "), ("All Files", "*.*")])

    if not path_spec:
        return
    with open(path_spec, mode="r", encoding="utf-8") as input_file:
        data, tau, H, errorH, path_spec = [], [], [], [], []
        data = np.loadtxt(input_file)

        tau = data[:, 0]
        H = data[:, 1]
        #errorH = data[:, 2]
        data_flag = 1
        menu1.entryconfig('Plot relaxation time spectrum', state='active')
        mwd_flag = 0
        
"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                            OPEN_MOD: OPENS DYNAMIC MODULI
                    
################################################################################################"""


def open_mod():

    global data_mod, omega, G1, G2, path_spec, mod_flag, mwd_flag_GEX

    path_spec = askopenfilename(
        filetypes=[("Text files", "*.txt *.dat *.sol "), ("All Files", "*.*")])

    if not path_spec:
        return
    with open(path_spec, mode="r", encoding="utf-8") as input_file:
        data_mod, omega, G1, G2, path_spec = [], [], [], [], []
        data_mod = np.loadtxt(input_file)

        omega = data_mod[:, 0]
        G1 = data_mod[:, 1]
        G2 = data_mod[:, 2]
        mod_flag = 1
        menu4.entryconfig('Plot dynamic moduli', state='active')
        mwd_flag_GEX = 0


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

            OPEN_GPC: OPENS MWD FILE OBTAINED WITH SEC/GPC (ON GAR TAB)
                        
################################################################################################"""


def open_gpc():

    global data_gpc, mass, distrib, path_gpc, gpc_flag

    path_gpc = askopenfilename(
        filetypes=[("Text files", "*.txt *.dat *.sol "), ("All Files", "*.*")])

    if not path_gpc:
        return
    with open(path_gpc, mode="r", encoding="utf-8") as input_file:
        data_gpc, mass, distrib, path_gpc = [], [], [], []
        data_gpc = np.loadtxt(input_file)

        mass = data_gpc[:, 0]
        distrib = data_gpc[:, 1]
        menu1.entryconfig('Plot MWD data only', state='active')
        gpc_flag = 1

        if mwd_flag == 1:
            menu1.entryconfig(
                'Plot MWD + estimated distribution', state='active')

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                OPEN_GPC: OPENS MWD FILE OBTAINED WITH SEC/GPC (ON GEX TAB)
                        
################################################################################################"""


def open_gpc_GEX():

    global data_gpc_GEX, mass_GEX, distrib_GEX, path_gpc_GEX, gpc_flag_GEX, mwd_flag_GEX

    path_gpc_GEX = askopenfilename(
        filetypes=[("Text files", "*.txt *.dat *.sol "), ("All Files", "*.*")])

    if not path_gpc_GEX:
        return
    with open(path_gpc_GEX, mode="r", encoding="utf-8") as input_file:
        data_gpc_GEX, mass_GEX, distrib_GEX, path_gpc_GEX = [], [], [], []
        data_gpc_GEX = np.loadtxt(input_file)

        mass_GEX = data_gpc_GEX[:, 0]
        distrib_GEX = data_gpc_GEX[:, 1]
        menu4.entryconfig('Plot MWD data only', state='active')
        gpc_flag_GEX = 1

        if mwd_flag_GEX == 1:
            menu4.entryconfig(
                'Plot MWD + estimated distribution', state='active')


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                PLOT_CLEAR: PLOTS BLANK FIGURE IN THE GRAPHICS AREA
            
################################################################################################"""


def plot_clear1():
    figure = Figure(figsize=(6.1, 4), dpi=90)

    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame1)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame1)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

def plot_clear2():
    figure = Figure(figsize=(6.1, 4), dpi=90)

    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame2)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame2)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                            PLOT_SPEC: PLOTS RELAXATION SPECTRUM
                    
################################################################################################"""


def plot_spec():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.loglog(tau, H)
    #axes.semilogx(tau, H)
    axes.set_xlabel(r'$\tau$ [s]')
    axes.set_ylabel(r'H($\tau$)')
    axes.set_title('Relaxation time spectrum')
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame1)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame1)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                        PLOT_MOD: PLOTS DYNAMIC MODULI
                    
################################################################################################"""


def plot_mod():
    figure = Figure(figsize=(plot_size_x, plot_size_y), dpi=90)
    axes = figure.add_subplot()
    axes.loglog(omega, G1,label='G\' - Storage modulus')
    axes.loglog(omega, G2,label='G\'\' - Loss modulus')
    axes.set_xlabel(r'$\omega$ [rad/s]')
    axes.set_ylabel(r'G*($\omega$)')
    axes.set_title('Dynamic moduli')
    axes.legend(loc='best')
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame2)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame2)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=plot_pos_x, y=plot_pos_y)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

            PLOT_MOD_GEX: PLOT FITTED DYNAMIC MODULI (GEX MODEL)
                    
################################################################################################"""


def plot_mod_GEX():
    figure = Figure(figsize=(plot_size_x, plot_size_y), dpi=90)
    axes = figure.add_subplot()
    axes.loglog(omega, G1,label='G\' - Storage modulus')
    axes.loglog(omega, G2,label='G\'\' - Loss modulus')
    axes.loglog(omega_GEX_mod, G1_GEX_mod,label='G\' model')
    axes.loglog(omega_GEX_mod, G2_GEX_mod,label='G\'\' model')
    axes.set_xlabel(r'$\omega$ [rad/s]')
    axes.set_ylabel(r'G*($\omega$)')
    axes.set_title('Dynamic moduli')
    axes.legend(loc='best')
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame2)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame2)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=plot_pos_x, y=plot_pos_y)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                    PLOT_GPC: PLOTS MWD OBTAINED WITH SEC/GPC (ON GAR TAB)
                        
################################################################################################"""


def plot_gpc():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.semilogx(mass, distrib)
    axes.set_xlabel('Molecular weight [g/mol]')
    axes.set_ylabel('w(M)')
    axes.set_title('MWD obtained by SEC/GPC')
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame1)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame1)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                PLOT_GPC: PLOTS MWD OBTAINED WITH SEC/GPC (ON GEX TAB)
                        
################################################################################################"""


def plot_gpc_GEX():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.semilogx(mass_GEX, distrib_GEX)
    axes.set_xlabel('Molecular weight [g/mol]')
    axes.set_ylabel('w(M)')
    axes.set_title('MWD obtained by SEC/GPC')
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame2)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame2)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                        PLOT_MWD: PLOTS ESTIMATED MWD (ON GAR TAB)
                
################################################################################################"""


def plot_mwd():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.semilogx(M, mwd)
    axes.set_xlabel('Molecular weight [g/mol]')
    axes.set_ylabel('w(M)')
    axes.set_title('Estimated MWD')
    axes.set_xlim(min(M), max(M))
    # axes.set_ylim(0,max(mwd)*1.1)
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame1)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame1)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                        PLOT_MWD: PLOTS ESTIMATED MWD (ON GEX TAB)
                
################################################################################################"""


def plot_mwd_GEX():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.semilogx(M_GEX, mwd_GEX)
    axes.set_xlabel('Molecular weight [g/mol]')
    axes.set_ylabel('w(M)')
    axes.set_title('Estimated MWD')
    axes.set_xlim(min(M_GEX), max(M_GEX))
    # axes.set_ylim(0,max(mwd)*1.1)
    # create FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(figure, plotsFrame2)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame2)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                    PLOT_BOTH: PLOTS GPC + ESTIMATED MWD (ON GAR TAB)
                
################################################################################################"""


def plot_both():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.semilogx(M, mwd, label='Estimated MWD')
    axes.semilogx(mass, distrib, label='GPC/SEC')
    axes.set_xlabel('Molecular weight [g/mol]')
    axes.set_ylabel('w(M)')
    axes.set_title('Estimated MWD + SEC/GPC data')
    axes.legend(loc='best')
    # create FigureCanvasTkAgg object

    canvas = FigureCanvasTkAgg(figure, plotsFrame1)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame1)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                PLOT_BOTH: PLOTS GPC + ESTIMATED MWD (ON GEX TAB)
                
################################################################################################"""


def plot_both_GEX():
    figure = Figure(figsize=(6.1, 4), dpi=90)
    axes = figure.add_subplot()
    axes.semilogx(M_GEX, mwd_GEX, label='Estimated MWD')
    axes.semilogx(mass_GEX, distrib_GEX, label='GPC/SEC')
    axes.set_xlabel('Molecular weight [g/mol]')
    axes.set_ylabel('w(M)')
    axes.set_title('Estimated MWD + SEC/GPC data')
    axes.legend(loc='best')
    # create FigureCanvasTkAgg object

    canvas = FigureCanvasTkAgg(figure, plotsFrame2)

    # create the toolbar
    tool = NavigationToolbar2Tk(canvas, frame2)
    canvas.get_tk_widget().place(x=30, y=15)
    tool.place(x=305, y=485)

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

        Extend scipy interp1d to interpolate/extrapolate per axis in log space              
        
################################################################################################"""


class interpolate1d(scipy.interpolate.interp1d):

    def __init__(self, x, y, *args, xspace='linear', yspace='linear', **kwargs):
        self.xspace = xspace
        self.yspace = yspace
        if self.xspace == 'log':
            x = np.log10(x)
        if self.yspace == 'log':
            y = np.log10(y)
        super().__init__(x, y, *args, **kwargs)

    def __call__(self, x, *args, **kwargs):
        if self.xspace == 'log':
            x = np.log10(x)
        if self.yspace == 'log':
            return 10**super().__call__(x, *args, **kwargs)
        else:
            return super().__call__(x, *args, **kwargs)


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

            GET_VALUES: READ AND VALIDATE VARIABLE VALUES (GAR MODEL)
            
################################################################################################"""


def get_values():
    global a, B, k, M0, G0, Me
    a, B, k, M0, G0, Me = [], [], [], [], [], []

    try:
        a = float(rept_a_GAR.get())
        B = float(rept_B_GAR.get())
        k = float(mat_k_GAR.get())
        M0 = float(mat_M0_GAR.get())
        G0 = float(mat_G0_GAR.get())
        Me = float(mat_Me_GAR.get())

        ok_flag = 1

        if not data_flag:
            showerror(title="Invalid RTS",
                      message="Please select a valid relaxation time spectrum file. " +
                      "For more info about the RTS file, check the \"Help\" option on \
                                the top menu.")
            ok_flag = 0

    except ValueError:
        showerror(title="Invalid values",
                  message="Please check if the values of the typed material and reptation model\
                        parameters are valid!")
        ok_flag = 0

    return ok_flag

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

            GET_VALUES: READ AND VALIDATE VARIABLE VALUES (GEX MODEL)
            
################################################################################################"""


def get_values_GEX():
    global a_GEX, B_GEX, k_GEX, M0_GEX, G0_GEX, Me_GEX
    a_GEX, B_GEX, k_GEX, M0_GEX, G0_GEX, Me_GEX = [], [], [], [], [], []

    try:
        a_GEX = float(rept_a_GEX.get())
        B_GEX = float(rept_B_GEX.get())
        k_GEX = float(mat_k_GEX.get())
        M0_GEX = float(mat_M0_GEX.get())
        G0_GEX = float(mat_G0_GEX.get())
        Me_GEX = float(mat_Me_GEX.get())

        ok_flag_GEX = 1

        if not mod_flag:
            showerror(title="Invalid dynamic modulii",
                      message="Please select a valid dynamic moduli file. " +
                      "For more info about the dynamic moduli file, check the \"Help\" option on \
                                the top menu.")
            ok_flag_GEX = 0

    except ValueError:
        showerror(title="Invalid values",
                  message="Please check if the values of the typed material and reptation model\
                        parameters are valid!")
        ok_flag_GEX = 0

    return ok_flag_GEX

def get_options():
    return os.listdir(script_dir+"/materials")


def opt_changed1(event):

    global opt1
    opt1 = []
    #showinfo(title='Result', message=selected_value.get())
    opt1 = selected_value1.get()

    mat_path = script_dir+"/materials/"+opt1

    with open(mat_path, mode="r", encoding="utf-8") as input_file:
        #temp = np.loadtxt(input_file)
        temp = input_file.readlines()

        rept_a_GAR.delete(0, 'end')
        rept_a_GAR.insert(0, float(temp[0]))

        rept_B_GAR.delete(0, 'end')
        rept_B_GAR.insert(0, float(temp[1]))

        mat_k_GAR.delete(0, 'end')
        mat_k_GAR.insert(0, float(temp[2]))

        mat_M0_GAR.delete(0, 'end')
        mat_M0_GAR.insert(0, float(temp[3]))

        mat_G0_GAR.delete(0, 'end')
        mat_G0_GAR.insert(0, float(temp[4]))

        mat_Me_GAR.delete(0, 'end')
        mat_Me_GAR.insert(0, float(temp[5]))

        note_text = temp[6:]

        Hovertip(selectMat1, note_text, hover_delay=delay_value)

        input_file.close()
        
def opt_changed2(event):

    global opt2
    opt2 = []
    #showinfo(title='Result', message=selected_value.get())
    opt2 = selected_value2.get()

    mat_path = script_dir+"/materials/"+opt2

    with open(mat_path, mode="r", encoding="utf-8") as input_file:
        #temp = np.loadtxt(input_file)
        temp = input_file.readlines()

        rept_a_GEX.delete(0, 'end')
        rept_a_GEX.insert(0, float(temp[0]))

        rept_B_GEX.delete(0, 'end')
        rept_B_GEX.insert(0, float(temp[1]))

        mat_k_GEX.delete(0, 'end')
        mat_k_GEX.insert(0, float(temp[2]))

        mat_M0_GEX.delete(0, 'end')
        mat_M0_GEX.insert(0, float(temp[3]))

        mat_G0_GEX.delete(0, 'end')
        mat_G0_GEX.insert(0, float(temp[4]))

        mat_Me_GEX.delete(0, 'end')
        mat_Me_GEX.insert(0, float(temp[5]))

        note_text = temp[6:]

        Hovertip(selectMat2, note_text, hover_delay=delay_value)

        input_file.close()
        
"""------------------------------------------------------------------------------------------------ 
###################################################################################################

        FUNÇÃO ROUSE_CORRECT: CALLS ROUSE SPECTRUM CORRECTION FUNCTION
    
################################################################################################"""

def rouse_correct():
    global rouseAdj, popt, new_tau, slopeRouse, rouseMult
    def func1(x, kRouse):
        return kRouse*x**(slopeRouse)
        #return kRouse*np.exp((slopeRouse)*np.log(x))
    
    def func2(x, kRouse, slope):
        return kRouse*x**(slope)
        #return kRouse*np.exp((slope)*np.log(x))
        
    lim = (k*(Me)**a)
    tau_aux = []
    H_aux = []
    idx = 0
    
    for i in tau:
        if i <= lim:
            tau_aux.append(i)
            H_aux.append(H[idx])
            idx = idx + 1
            
    new_tau = []
    lim2 = (k*(limRouse*Me)**a)
    idx = 0
    
    for i in tau:
        if i <= lim2:
            new_tau.append(i)
            idx = idx + 1
    
    popt, pcov = [], []
    rouseAdj = []
    
    if optRouse == 1 or optRouse == 3:
        popt, pcov = curve_fit(func1, tau_aux, H_aux, p0=1)
        rouseAdj = func1(np.array(new_tau),popt[0])
        rouseMult = popt[0]
    elif optRouse == 2:
        popt, pcov = curve_fit(func2, tau_aux, H_aux, p0=(1, -0.5))
        rouseAdj = func2(np.array(new_tau),popt[0],popt[1])
        rouseMult = popt[0]
        slopeRouse = popt[1]
    
    H_new = []
    H_new = H - np.append(rouseAdj,np.zeros(len(H)-len(rouseAdj)))
    
    for i in range(len(H_new)):
        if (H_new[i]) < 0:
            H_new[i] = 0

    for i in range(len(H_new)):
        if ((tau[i]/k)**(1/a) - Me) < Me:
            H_new[i] = 0
            
  
    return H_new

"""------------------------------------------------------------------------------------------------        
###################################################################################################

                            CREATING THE ROUSE PARAMETERS WINDOW
                    
################################################################################################"""


def configRouse():
    
    if checkRouseVal.get() == 0:
        menu3.entryconfig('Open Rouse configurations', state='disabled')
        return
    else:
        def close_window3():
            if checkRouseVal.get() == 1:
                root3.destroy()
            else:
                checkRouseVal.set(0)
                root3.destroy()

            
        def get_state():
            if selected.get() == '3':
                slope.config(state='active')
            else:
                slope.config(state='disabled')
                
        def get_Rouse_values():
            global optRouse, limRouse, slopeRouse
            
            optRouse = int(selected.get())
            
            try:
                limRouse = float(limEntry.get())
                if optRouse == 3:    
                    try: 
                        slopeRouse = float(slope.get())
                        root3.destroy()
                    
                    except ValueError:
                        showerror(title="Invalid values",
                                  message="Please check if the typed values are valid!")
                    return
                elif optRouse == 1:
                    slopeRouse = -0.5
                    root3.destroy()
                elif optRouse == 2:
                    root3.destroy()
                    
            except ValueError:
                showerror(title="Invalid values",
                          message="Please check if the typed values are valid!")
                
            
                   
        def create_controls_frame3(container):
            frame = ttk.Frame(container)
            
            OK_button = ttk.Button(frame, text='  OK  ', command=get_Rouse_values)
            OK_button.grid(column=0, row=0)

            cancel_button = ttk.Button(frame, text=' Cancel ', command=close_window3)
            cancel_button.grid(column=1, row=0)

            for widget in frame.winfo_children():
                widget.grid(padx=5, pady=10)

            return frame  
        
        selected = tk.StringVar()
                    
        root3 = tk.Toplevel(root)
        root3.wm_transient(root)
    
        #root2.attributes('-topmost', 'true')
    
        root3.title('Rouse model settings')
    
        window_width = 300
        window_height = 285
    
        # get the screen dimension
        screen_width = root3.winfo_screenwidth()
        screen_height = root3.winfo_screenheight()
    
        # find the center point
        center_x = int(screen_width/2 - window_width / 2)
        center_y = int(screen_height/2 - window_height / 2)
    
        # set the position of the window to the center of the screen
        root3.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
    
        #root2.geometry("%dx%d" % (screen_width-2, screen_height-10))
    
        root3.resizable(False, False)
        
        radioFrame = ttk.LabelFrame(root3, text=' Rouse spectrum slope ', height=150, width=260)
        radioFrame.place(x=20, y=10)
        
        r1 = ttk.Radiobutton(root3, text='Use standard -1/2 slope', value=1, variable=selected, command=get_state)
        r1.place(x=35, y=30)
        r2 = ttk.Radiobutton(root3, text='Use best numerical fitting', value=2, variable=selected, command=get_state)
        r2.place(x=35, y=60)
        r3 = ttk.Radiobutton(root3, text='Use given value:', value=3, variable=selected, command=get_state)            
        r3.place(x=35, y=90)
        slope = ttk.Entry(root3, width=32, state='disabled')
        slope.place(x=35, y=125)
        
        limLabel = ttk.Label(root3, text='Rouse spectrum subtraction limit:')
        limLabel.place(x=20, y=185)
        Hovertip(limLabel, 'Value in multiples of Me',hover_delay=delay_value)
        limEntry = ttk.Entry(root3, width=36,)
        limEntry.place(x=22, y=210)
        limEntry.insert(0,'22')
        Hovertip(limEntry, 'Value in multiples of Me',hover_delay=delay_value)
        
        controls_frame3 = create_controls_frame3(root3)
        controls_frame3.place(x=65,y=235)
        
        if not optRouse or optRouse == 1:
            r1.invoke()
        elif optRouse == 2:
            r2.invoke()
        elif optRouse == 3:
            r3.invoke()
            slope.insert(0,slopeRouse)
        
        limEntry.delete(0,'end')
        limEntry.insert(0,limRouse)
        
        menu3.entryconfig('Open Rouse configurations', state='active')
           
        root3.protocol("WM_DELETE_WINDOW", close_window3)
        root3.mainloop()

"""------------------------------------------------------------------------------------------------        
###################################################################################################

                            CREATION OF THE GEX CONFIGURATION WINDOW
                    
################################################################################################"""
def closing_configLim():
    
    if checkLim.get() == 1:
        configGEX()
    elif checkLim.get() == 0 and checkWeight.get() == 0 and checkGuess.get() == 0:
        menu6.entryconfig('Open configurations window', state='disabled')

def closing_configWeight():
    if checkWeight.get() == 1:
        configGEX()
    elif checkLim.get() == 0 and checkWeight.get() == 0 and checkGuess.get() == 0:
        menu6.entryconfig('Open configurations window', state='disabled')

def closing_configGuess():
    if checkGuess.get() == 1:
        configGEX()
    elif checkLim.get() == 0 and checkWeight.get() == 0 and checkGuess.get() == 0:
        menu6.entryconfig('Open configurations window', state='disabled')


def configGEX():
    
    def close_window4():
        root4.destroy()
    
    def create_controls_frame4(container):
        frame = ttk.Frame(container)
        
        OK_button = ttk.Button(frame, text='  OK  ', command=get_GEX_values)
        OK_button.grid(column=0, row=0)

        cancel_button = ttk.Button(frame, text=' Cancel ', command=close_window4)
        cancel_button.grid(column=1, row=0)

        for widget in frame.winfo_children():
            widget.grid(padx=5, pady=10)

        return frame
    
    def get_GEX_values():
        global freqLim, weight, init_a, init_b, init_m0
        
        if checkLim.get() == 1:
            try:
                freqLim = float(limEntry1.get())
                #print(freqLim)
            except ValueError:
                showerror(title="Invalid values",
                          message="Please check if the typed values are valid!")
                
        if checkWeight.get() == 1:
            try:
                weight = float(weightEntry1.get())

                #print(wG1)
                #print(wG2)
            except ValueError:
                showerror(title="Invalid values",
                          message="Please check if the typed values are valid!")
                
        if checkGuess.get() == 1:
            try:
                init_a = float(guessEntry1.get())
                init_b = float(guessEntry2.get())
                init_m0 = float(guessEntry3.get())
                #print(init_a)
                #print(init_b)
                #print(init_m0)
            except ValueError:
                showerror(title="Invalid values",
                          message="Please check if the typed values are valid!")
    
        root4.destroy() 
                
    root4 = tk.Toplevel(root)
    root4.wm_transient(root)


    root4.title('GEX settings')

    window_width = 300
    window_height = 295

    # get the screen dimension
    screen_width = root4.winfo_screenwidth()
    screen_height = root4.winfo_screenheight()

    # find the center point
    center_x = int(screen_width/2 - window_width / 2)
    center_y = int(screen_height/2 - window_height / 2)

    # set the position of the window to the center of the screen
    root4.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
    root4.resizable(False, False)
    
    limFrame = ttk.LabelFrame(root4, text=' Upper frequency window limit ', height=55, width=260)
    limFrame.place(x=20, y=10)
    Hovertip(limFrame, 'Value in rad/s',hover_delay=delay_value)
    limEntry1 = ttk.Entry(root4, width=32)
    limEntry1.place(x=30, y=35)
    if checkLim.get() == 0:
        limEntry1.config(state='readonly')
    else:
        limEntry1.delete(0,'end')
        limEntry1.insert(0,freqLim)
    
    weightFrame = ttk.LabelFrame(root4, text=' Relative weight ', height=55, width=260)
    weightFrame.place(x=20, y=75)
    Hovertip(weightFrame, 'Values between 0 and 1',hover_delay=delay_value)
    weightLabel1 = ttk.Label(root4, text=' w:')
    weightLabel1.place(x=30, y=100)
    weightEntry1 = ttk.Entry(root4, width=27)
    weightEntry1.place(x=64, y=100)

    if checkWeight.get() == 0:
        weightEntry1.config(state='readonly')
    else:
        weightEntry1.delete(0,'end')
        weightEntry1.insert(0,weight)

    
    guessFrame = ttk.LabelFrame(root4, text=' Initial guesses ', height=105, width=260)
    guessFrame.place(x=20, y=140)
    Hovertip(guessFrame, 'For more info about this section, check the "Help" option on the top menu',
             hover_delay=delay_value)
    guessLabel1 = ttk.Label(root4, text=' a:')
    guessLabel1.place(x=30, y=165)
    guessEntry1 = ttk.Entry(root4, width=27)
    guessEntry1.place(x=64, y=165)
    guessLabel2 = ttk.Label(root4, text=' b:')
    guessLabel2.place(x=30, y=190)
    guessEntry2 = ttk.Entry(root4, width=27)
    guessEntry2.place(x=64, y=190)
    guessLabel3 = ttk.Label(root4, text='m0:')
    guessLabel3.place(x=30, y=215)
    Hovertip(guessLabel3, 'Value multiplied by 1e5',
             hover_delay=delay_value)
    guessEntry3 = ttk.Entry(root4, width=27)
    guessEntry3.place(x=64, y=215)
    if checkGuess.get() == 0:
        guessEntry1.config(state='readonly')
        guessEntry2.config(state='readonly')
        guessEntry3.config(state='readonly')
    else:
        guessEntry1.delete(0,'end')
        guessEntry1.insert(0,init_a)
        guessEntry2.delete(0,'end')
        guessEntry2.insert(0,init_b)
        guessEntry3.delete(0,'end')
        guessEntry3.insert(0,init_m0)
    
    controls_frame4 = create_controls_frame4(root4)
    controls_frame4.place(x=65,y=250)
    
    if (checkLim.get() == 1 or checkWeight.get() == 1 or checkGuess.get() == 1):
        menu6.entryconfig('Open settings window', state='active')
       
    root4.protocol("WM_DELETE_WINDOW", close_window4)
    root4.mainloop()


def consoleWindow():
                
    root5 = tk.Toplevel(root)
    root5.wm_transient(root)


    root5.title('GEX console')

    window_width = 400
    window_height = 325

    # get the screen dimension
    screen_width = root5.winfo_screenwidth()
    screen_height = root5.winfo_screenheight()

    # find the center point
    center_x = int(screen_width/2 - window_width / 2)
    center_y = int(screen_height/2 - window_height / 2)

    # set the position of the window to the center of the screen
    root5.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
    root5.resizable(False, False)
    
    consoleBox = tk.Text(root5)
    consoleBox.pack(pady=5, padx=5, fill='both')
    consoleBox.delete(1.0, 'end')
    consoleBox.insert(1.0, report)
    consoleBox.config(state='disabled')  
    
    #root4.protocol("WM_DELETE_WINDOW", close_window5)
    root5.mainloop()
    
"""------------------------------------------------------------------------------------------------ 
###################################################################################################

    RUN_LIKE_HELL: CALLS VARIABLE READING FUNCTION AND CALCULATES MWD WITH GAR MODEL
    
################################################################################################"""

def run_like_hell():
    st = time.time()
    global mwd_flag, new_H
    global Mn, Mw, Mz, MwMn
    global M, mwd
    global H2

    is_ok = get_values()

    if is_ok == 1:
        
        H2 = np.zeros(len(H))
        
        if checkRouseVal.get() == 1:
            H2 = rouse_correct()
        else:
            H2 = H.copy()
            
        aux_H = np.zeros(len(H2))

        if checkMe1.get() == 1:
            Me_aux1 = Me
        else:
            Me_aux1 = 0

        for i in range(len(H2)):
            if (((tau[i]/k)**(1/a) - Me_aux1) > Me):
                aux_H[i] = H2[i]
                
        new_H = []
        M_aux = []
        M = []
        nz = []

        nz = np.nonzero(aux_H)[0]

        new_H = aux_H[min(nz):max(nz)]
        
        M_aux = ((tau[min(nz):max(nz)]/(k))**(1/a))

        #if checkMe2.get() == 0:
        #    M_aux = M_aux - Me
        
        m = M_aux/M0

        M = M_aux[0:len(M_aux)-1]
        
        interp_func = interpolate1d(m, new_H/m, xspace='log', fill_value='extrapolate')
        
        
        H_integ = []

        for i in range(len(m)-1):
            H_integ.append(quad(interp_func, m[i], max(m), limit=250)[0])
        
        mwd_temp = []
        mwd_temp = (1/B)*((a/G0)**(1/B))*(np.array(H_integ)
                                     ** (1/B-1))*new_H[0:len(new_H)-1]
        
        mwd_interp = interpolate1d(M, mwd_temp/M, xspace='log', fill_value='extrapolate')
        norm_value_mwd = quad(mwd_interp,min(M),max(M),limit=250)[0]
        
        mwd = mwd_temp/norm_value_mwd

        plot_mwd()

        menu1.entryconfig('Plot estimated distribution', state='active')

        mwd_flag = 1

        if gpc_flag == 1:
            menu1.entryconfig(
                'Plot MWD + estimated distribution', state='active')

        Mn, Mw, Mz, MwMn = [], [], [], []
        mwd_vr1, mwd_vr2, mwd_vr3, mwd_vr4 = [], [], [], []
        int_mwd_vr1, int_mwd_vr2, int_mwd_vr3, int_mwd_vr4 = [], [], [], []

        mwd_vr1 = interpolate1d(M, M*mwd, xspace='log',
                                fill_value='extrapolate')
        mwd_vr2 = interpolate1d(M, mwd, xspace='log', fill_value='extrapolate')
        mwd_vr3 = interpolate1d(M, mwd/M, xspace='log',
                                fill_value='extrapolate')
        mwd_vr4 = interpolate1d(
            M, (M**2)*mwd, xspace='log', fill_value='extrapolate')

        int_mwd_vr1 = quad(mwd_vr1, min(M), max(M), limit=250)[0]
        int_mwd_vr2 = quad(mwd_vr2, min(M), max(M), limit=250)[0]
        int_mwd_vr3 = quad(mwd_vr3, min(M), max(M), limit=250)[0]
        int_mwd_vr4 = quad(mwd_vr4, min(M), max(M), limit=250)[0]

        Mn = int_mwd_vr2/int_mwd_vr3
        Mw = int_mwd_vr1/int_mwd_vr2
        Mz = int_mwd_vr4/int_mwd_vr1
        MwMn = Mw/Mn

        resMn.config(state='active')
        resMw.config(state='active')
        resMz.config(state='active')
        resMwMn.config(state='active')

        resMn.delete(0, 'end')
        resMw.delete(0, 'end')
        resMz.delete(0, 'end')
        resMwMn.delete(0, 'end')

        resMn.insert(0, str("{:.0f}".format(Mn)))
        resMw.insert(0, str("{:.0f}".format(Mw)))
        resMz.insert(0, str("{:.0f}".format(Mz)))
        resMwMn.insert(0, str("{:.2f}".format(MwMn)))

        resMn.config(state='readonly')
        resMw.config(state='readonly')
        resMz.config(state='readonly')
        resMwMn.config(state='readonly')

        menu2.entryconfig('Copy to clipboard', state='active')

    else:
        pass

    et = time.time()
    elapsed_time = et - st
    #print('Execution time:', elapsed_time, 'seconds')
    #print('Normalization value:',norm_value_mwd)
    

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

    RUN_FOR_YOUR_LIFE: CALLS VARIABLE READING FUNCTION AND CALCULATES MWD WITH GEX MODEL
    
################################################################################################"""

def run_for_your_life():
    st = time.time()
    global mwd_flag_GEX, gpc_flag_GEX
    global Mn_GEX, Mw_GEX, Mz_GEX, MwMn_GEX
    global M_GEX, mwd_GEX
    global new_omega, new_G1, new_G2, comboY
    global a_par, b_par, m0_par
    global rTime
    global G1_GEX_mod, G2_GEX_mod, omega_GEX_mod
    
    global report, fitting_report

    is_ok_GEX = get_values_GEX()

    if is_ok_GEX == 1:
        
        Mn_GEX = []
        Mw_GEX = []
        Mz_GEX = []
        MwMn_GEX = []
        
        new_omega = []
        new_G1 = []
        new_G2 = []
        
        comboY = []
        report, fitting_report = [], []
        
        if checkLim.get() == 1:
            for i in range(len(omega)-1):
                if omega[i] <= freqLim:
                    new_omega = np.append(new_omega, omega[i])
                    new_G1 = np.append(new_G1, G1[i])
                    new_G2 = np.append(new_G2, G2[i])
            comboY = np.append(new_G1,new_G2)
            #comboY = np.add(new_G1,new_G2)
        else:
            comboY = np.append(G1,G2)
            #comboY = np.add(G1,G2)
            new_omega = omega.copy()

        
        M_GEX, rTime = [], []
        M_GEX = np.logspace(4,7,100)
        rTime = k_GEX*M_GEX**a_GEX
        #tau_GEX = k_GEX*M_GEX**a_GEX
        
        def Hgex(rTime,params):
            a = params['a'].value
            b = params['b'].value
            m0 = 1e5*params['m0'].value 
            
            mass_temp = np.array((rTime/k_GEX)**(1/a_GEX))
            
            return G0_GEX*(B_GEX/a_GEX)*(b/((sc.gamma((a+1)/b))**B_GEX))*((mass_temp/m0)**(a+1))*np.exp(-(mass_temp/m0)**b)*(sc.gamma((a+1)/b)*(sc.gammaincc((a+1)/b,(mass_temp/m0)**b))**(B_GEX-1))
        
#-------------------------------------------------------------------------------------------------------------------------------

        def storage_mod(params, omega): 
            a = params['a'].value
            b = params['b'].value
            m0 = params['m0'].value 
            
            def integrand(rTime):
                return (Hgex(rTime,params)/rTime)*(((w*rTime)**2)/(1+((w*rTime)**2)))
                
            model = []

            for w in omega:
                model.append(quad(integrand,0,np.inf,limit=250)[0])

            return model

#-------------------------------------------------------------------------------------------------------------------------------

        def loss_mod(params, omega): 
            a = params['a'].value
            b = params['b'].value
            m0 = params['m0'].value 
            
            def integrand(rTime):
                return (Hgex(rTime,params)/rTime)*(((w*rTime))/(1+((w*rTime)**2)))
    
            model = []
 
            for w in omega:
                model.append(quad(integrand,0,np.inf,limit=250)[0])
            
            return model
        
#-------------------------------------------------------------------------------------------------------------------------------

        def func2min1(params, omega, dataY): 
            a = params['a'].value
            b = params['b'].value
            m0 = params['m0'].value 
            
            def integrand(rTime):
                return (Hgex(rTime,params)/rTime)*(((w*rTime)**2)/(1+((w*rTime)**2)))
                
            
            model = []

            for w in omega:
                model.append(quad(integrand,0,np.inf,limit=250)[0])
             
            
            resids = np.array(model) - dataY
            weighted = np.sqrt(resids**2 / dataY**2)
            #weighted = resids**2 / dataY
            
            return weighted

#-------------------------------------------------------------------------------------------------------------------------------

        def func2min2(params, omega, dataY): 
            a = params['a'].value
            b = params['b'].value
            m0 = params['m0'].value 
            
            def integrand(rTime):
                return (Hgex(rTime,params)/rTime)*(((w*rTime))/(1+((w*rTime)**2)))
    
            model = []
 
            for w in omega:
                model.append(quad(integrand,0,np.inf,limit=250)[0])
               
            resids = np.array(model) - dataY
            weighted = np.sqrt(resids**2 / dataY**2)
            #weighted = resids**2 / dataY
            
            return weighted
        
#-------------------------------------------------------------------------------------------------------------------------------
         
        def func2min3(params, omega, dataY): 
            a = params['a'].value
            b = params['b'].value
            m0 = params['m0'].value 
            
            def integrand(rTime):
                return (Hgex(rTime,params)/rTime)*((((w*rTime)**2) + (w*rTime))/(1+((w*rTime)**2)))
            
            model = []

            for w in omega:
                model.append(quad(integrand,0,np.inf,limit=250)[0])
               
            resids = np.array(model) - dataY
            weighted = np.sqrt(resids**2 / dataY**2)
            
            return weighted

#-------------------------------------------------------------------------------------------------------------------------------
        
        
        def function(params, dataX, comboY, factor):
            result1 = np.array(func2min1(params, dataX, comboY[:len(dataX)]))
            result2 = np.array(func2min2(params, dataX, comboY[len(dataX):]))
            #return np.add(factor1*result1, factor2*result2)
            return np.append(result1, factor*result2)
        
        """        
        def function(params,dataX,comboY):
            result = func2min3(params, dataX, comboY)
            return result
        """  
        
        def per_iteration(pars, iteration, resid, *args, **kws):
            print(" ITER ", iteration, " RESIDUAL ", np.sum(resid), [f"{p.name} = {p.value:.5f}" for p in pars.values()])
        

        params, fitter, result_ampgo = [], [], []
        
        params = Parameters()
        params.add('a', value=init_a, min=0.1, max=init_a+2, vary=True)
        params.add('b', value=init_b, min=0.1, max=init_b+2, vary=True)
        params.add('m0', value=init_m0, min=0.01, max=init_m0+2, vary=True)
        
        kws={'disp':False, 'totaliter':10}
        
        fitter = Minimizer(function, params, fcn_args=(new_omega, comboY, weight), iter_cb=per_iteration)
        
        result_ampgo = fitter.minimize() #method='ampgo',**kws)
        print(fit_report(result_ampgo,min_correl=0.9))
        
        a_nom = result_ampgo.params['a'].value
        b_nom = result_ampgo.params['b'].value
        m0_nom = result_ampgo.params['m0'].value
        
        a_dev = result_ampgo.params['a'].stderr
        b_dev = result_ampgo.params['b'].stderr
        m0_dev = result_ampgo.params['m0'].stderr
        
        a_par = ufloat(a_nom, a_dev)
        b_par = ufloat(b_nom, b_dev)
        m0_par = ufloat(m0_nom, m0_dev)
        
        Mn_GEX = 1e5*m0_par*(gamma((a_par+1)/b_par)/gamma(a_par/b_par))
        Mw_GEX = 1e5*m0_par*(gamma((a_par+2)/b_par)/gamma((a_par+1)/b_par))
        Mz_GEX = 1e5*m0_par*(gamma((a_par+3)/b_par)/gamma((a_par+2)/b_par))
        MwMn_GEX = Mw_GEX/Mn_GEX
        
        resMn_GEX.config(state='active')
        resMw_GEX.config(state='active')
        resMz_GEX.config(state='active')
        resMwMn_GEX.config(state='active')

        resMn_GEX.delete(0, 'end')
        resMw_GEX.delete(0, 'end')
        resMz_GEX.delete(0, 'end')
        resMwMn_GEX.delete(0, 'end')

        resMn_GEX.insert(0, str("{:.0f}".format(Mn_GEX)))
        resMw_GEX.insert(0, str("{:.0f}".format(Mw_GEX)))
        resMz_GEX.insert(0, str("{:.0f}".format(Mz_GEX)))
        resMwMn_GEX.insert(0, str("{:.2f}".format(MwMn_GEX)))

        resMn_GEX.config(state='readonly')
        resMw_GEX.config(state='readonly')
        resMz_GEX.config(state='readonly')
        resMwMn_GEX.config(state='readonly')
        
        mwd_GEX, mwd_GEX_temp, mwd_GEX_interp, norm_value_mwd_GEX  = [], [], [], []
        
        mwd_GEX_temp = (b_nom/gamma((a_nom+1)/b_nom))*((M_GEX/(1e5*m0_nom))**(a_nom+1))*np.exp(-(M_GEX/(1e5*m0_nom))**b_nom)
        
        mwd_GEX_interp = interpolate1d(M_GEX, mwd_GEX_temp/M_GEX, xspace='log', fill_value='extrapolate')
        norm_value_mwd_GEX = quad(mwd_GEX_interp,min(M_GEX),max(M_GEX),limit=250)[0]
        
        mwd_GEX = mwd_GEX_temp/norm_value_mwd_GEX
        
        plot_mwd_GEX()

        menu4.entryconfig('Plot estimated distribution', state='active')
        
        omega_GEX_mod = []
        G1_GEX_mod, G2_GEX_mod = [], []
        
        omega_GEX_mod = new_omega.copy()
        G1_GEX_mod = storage_mod(result_ampgo.params, new_omega)
        G2_GEX_mod = loss_mod(result_ampgo.params, new_omega)
        
        
        menu4.entryconfig('Plot dynamic moduli + fitted moduli', state='active')
        
        mwd_flag_GEX = 1

        if gpc_flag_GEX == 1:
            menu4.entryconfig(
                'Plot MWD + estimated distribution', state='active')

        menu5.entryconfig('Copy to clipboard', state='active')

    else:
        pass

    et = time.time()
    elapsed_time = et - st
    #print('Execution time:', elapsed_time, 'seconds')
    #print('Normalization value:',norm_value_mwd_GEX)

    #report = ("### DONE! ###" + '\n' + str(fit_report(result_ampgo)) + '\n' + time_message)
    
    if mwd_flag_GEX == 1:
        fitting_report = fit_report(result_ampgo)
        report = '### DONE! ###'+'\n'+'\n'+fitting_report+'\n'+'\n'+'Execution time: '+str(elapsed_time)+' seconds'
        consoleWindow()


"""------------------------------------------------------------------------------------------------        
###################################################################################################

                COPY_VALUES: COPY NUMERIC RESULTS TO THE CLIPBOARD
                    
################################################################################################"""


def copy_values():
    clip_Mn, clip_Mw, clip_Mz, clip_MwMn = [], [], [], []

    root.clipboard_clear()
    clip_Mn = resMn.get()
    clip_Mw = resMw.get()
    clip_Mz = resMz.get()
    clip_MwMn = resMwMn.get()
    root.clipboard_append(clip_Mn+'\n'+clip_Mw+'\n'+clip_Mz+'\n'+clip_MwMn)
    
def copy_values_GEX():
    clip_Mn_GEX, clip_Mw_GEX, clip_Mz_GEX, clip_MwMn_GEX = [], [], [], []

    root.clipboard_clear()
    clip_Mn_GEX = resMn_GEX.get()
    clip_Mw_GEX = resMw_GEX.get()
    clip_Mz_GEX = resMz_GEX.get()
    clip_MwMn_GEX = resMwMn_GEX.get()
    root.clipboard_append(clip_Mn_GEX+'\n'+clip_Mw_GEX+'\n'+clip_Mz_GEX+'\n'+clip_MwMn_GEX)


"""------------------------------------------------------------------------------------------------        
###################################################################################################

                        EXPORT_MWD: EXPORT MWD TO TEXT FILE
                    
################################################################################################"""


def export_mwd():
    if mwd_flag == 1:
        path_to_pref = asksaveasfilename(defaultextension='.txt', filetypes=[
                                         ("txt files", '*.txt')], title="Choose filename")
        if not path_to_pref:
            return
        with open(path_to_pref, 'w') as f:
            for i in range(len(mwd)-1):
                d2w = []
                d2w = str(M[i])+' '+str(mwd[i])+'\n'
                f.write(str(d2w))
        f.close
    else:
        showwarning(title="Nothing to export!",
                    message="There is no calculated MWD to export.")
        
def export_mwd_GEX():
    if mwd_flag_GEX == 1:
        path_to_pref = asksaveasfilename(defaultextension='.txt', filetypes=[
                                         ("txt files", '*.txt')], title="Choose filename")
        if not path_to_pref:
            return
        with open(path_to_pref, 'w') as f:
            for i in range(len(mwd_GEX)-1):
                d2w = []
                d2w = str(M_GEX[i])+' '+str(mwd_GEX[i])+'\n'
                f.write(str(d2w))
        f.close
    else:
        showwarning(title="Nothing to export!",
                    message="There is no calculated MWD to export.")
        
"""------------------------------------------------------------------------------------------------        
###################################################################################################

    hidden function >>> EXPORT_MODULI: EXPORT FITTED G*(w) TO TEXT FILE
                    
################################################################################################"""       
        
def export_moduli():
    if mwd_flag_GEX == 1:
        path_to_pref = asksaveasfilename(defaultextension='.txt', filetypes=[
                                         ("txt files", '*.txt')], title="Choose filename")
        if not path_to_pref:
            return
        with open(path_to_pref, 'w') as f:
            for i in range(len(omega_GEX_mod)-1):
                d2w = []
                d2w = str(omega_GEX_mod[i])+' '+str(G1_GEX_mod[i])+' '+str(G2_GEX_mod[i])+'\n'
                f.write(str(d2w))
        f.close
    else:
        showwarning(title="Nothing to export!",
                    message="There is no fitted moduli to export.")


"""------------------------------------------------------------------------------------------------        
###################################################################################################

                    SAVE_LOG: SAVE CURRENT ANALYSIS LOG FILE
                    
################################################################################################"""


def save_log():
    if mwd_flag == 1:
        path_to_pref = asksaveasfilename(defaultextension='.txt', filetypes=[
                                         ("txt files", '*.txt')], title="Choose filename")
        if not path_to_pref:
            return
        with open(path_to_pref, mode="w", encoding="utf-8") as f:
            np.savetxt(f, [a, B, k, M0, G0, Me])
            np.savetxt(f, [checkMe1.get(), checkRouseVal.get(), checkMe2.get()])
            np.savetxt(f, [optRouse, limRouse, slopeRouse, rouseMult])
            np.savetxt(f, [Mn, Mw, Mz, MwMn])
            np.savetxt(f, [len(data)])
            np.savetxt(f, data)
            np.savetxt(f, np.c_[M, mwd])
        f.close

    else:
        showwarning(title="No current analysis!",
                    message="There is no current analysis to save or the loaded relaxation time spectrum does not match to the current analysis")

def save_log_GEX():
    if mwd_flag_GEX == 1:
        path_to_pref = asksaveasfilename(defaultextension='.txt', filetypes=[
                                         ("txt files", '*.txt')], title="Choose filename")
        if not path_to_pref:
            return
        with open(path_to_pref, mode="w", encoding="utf-8") as f:
            np.savetxt(f, [a_GEX, B_GEX, k_GEX, M0_GEX, G0_GEX, Me_GEX])
            np.savetxt(f, [checkLim.get(), checkWeight.get(), checkGuess.get()])
            np.savetxt(f, [freqLim, weight, init_a, init_b, init_m0])
            np.savetxt(f, [Mn_GEX.nominal_value, Mn_GEX.std_dev, Mw_GEX.nominal_value, Mw_GEX.std_dev, Mz_GEX.nominal_value, Mz_GEX.std_dev, MwMn_GEX.nominal_value, MwMn_GEX.std_dev])
            np.savetxt(f, [a_par.nominal_value, a_par.std_dev, b_par.nominal_value, b_par.std_dev, m0_par.nominal_value, m0_par.std_dev])
            np.savetxt(f, [len(data_mod)])
            np.savetxt(f, data_mod)
            np.savetxt(f, np.c_[M_GEX, mwd_GEX])
            #np.savetxt(f, str(fitting_report))
        f.close

    else:
        showwarning(title="No current analysis!",
                    message="There is no current analysis to save or the loaded relaxation time spectrum does not match to the current analysis")


"""------------------------------------------------------------------------------------------------        
###################################################################################################

                            OPEN_LOG: OPEN ANALYSIS LOG FILE
                    
################################################################################################"""


def open_log():

    global tau, H, errorH, data_flag, data
    global M, mwd, mwd_flag
    global optRouse, limRouse, slopeRouse, rouseMult

    answer = askokcancel(title='Open log file',
                         message='If you open a log file, all your current data will be lost.\
                         Do you wish to proceed?')

    if answer:

        path_to_pref = askopenfilename(
            filetypes=[("Text files", "*.txt *.dat *.sol "), ("All Files", "*.*")])
        if not path_to_pref:
            return
        with open(path_to_pref, mode="r", encoding="utf-8") as f:
            clear_all_GAR()
            data_flag = 1

            temp, dist = [], []
            
            temp = np.loadtxt(islice(f, 18))
            data = np.loadtxt(islice(f, int(temp[17])))
            dist = np.loadtxt(f)

            rept_a_GAR.insert(0,temp[0])
            rept_B_GAR.insert(0,temp[1])
            mat_k_GAR.insert(0,temp[2])
            mat_M0_GAR.insert(0,temp[3])
            mat_G0_GAR.insert(0,temp[4])
            mat_Me_GAR.insert(0,temp[5])
            
            checkMe1.set(int(temp[6]))
            checkRouseVal.set(int(temp[7]))
            checkMe2.set(int(temp[8]))
            
            optRouse = int(temp[9])
            limRouse = float(temp[10])
            slopeRouse = float(temp[11])
            rouseMult = float(temp[12])
            
            resMn.config(state='active')
            resMw.config(state='active')
            resMz.config(state='active')
            resMwMn.config(state='active')
            
            resMn.insert(0,str("{:.0f}".format(float(temp[13])))+un)
            resMw.insert(0,str("{:.0f}".format(float(temp[14])))+un)
            resMz.insert(0,str("{:.0f}".format(float(temp[15])))+un)
            resMwMn.insert(0,format(float(temp[16]), '.2f'))
            
            resMn.config(state='readonly')
            resMw.config(state='readonly')
            resMz.config(state='readonly')
            resMwMn.config(state='readonly')
            
            menu2.entryconfig('Copy to clipboard', state='active')
            
            tau = data[:,0]
            H = data[:,1]
            errorH = data[:,2]
            
            menu1.entryconfig('Plot relaxation time spectrum', state='active')
            
            M = dist[:,0]
            mwd = dist[:,1]
            
            plot_mwd()

            menu1.entryconfig('Plot estimated distribution', state='active')

            mwd_flag = 1

            if gpc_flag == 1:
                menu1.entryconfig(
                    'Plot MWD + estimated distribution', state='active')

        f.close

def open_log_GEX():

    global omega, G1, G2, mod_flag, data_mod
    global M_GEX, mwd_GEX, mwd_flag_GEX
    global Mn_GEX, Mw_GEX, Mz_GEX, MwMn_GEX
    global freqLim, weight, init_a, init_b, init_m0
    global a_par, b_par, m0_par
    global report

    answer = askokcancel(title='Open log file',
                         message='If you open a log file, all your current data will be lost.\
                         Do you wish to proceed?')

    if answer:

        path_to_pref = askopenfilename(
            filetypes=[("Text files", "*.txt *.dat *.sol "), ("All Files", "*.*")])
        if not path_to_pref:
            return
        with open(path_to_pref, mode="r", encoding="utf-8") as f:
            clear_all_GEX()
            mod_flag = 1

            temp, dist, data_mod = [], [], []
            
            temp = np.loadtxt(islice(f, 29))
            data_mod = np.loadtxt(islice(f, int(temp[28])))
            dist = np.loadtxt(f)

            rept_a_GEX.insert(0,temp[0])
            rept_B_GEX.insert(0,temp[1])
            mat_k_GEX.insert(0,temp[2])
            mat_M0_GEX.insert(0,temp[3])
            mat_G0_GEX.insert(0,temp[4])
            mat_Me_GEX.insert(0,temp[5])
            
            checkLim.set(int(temp[6]))
            checkWeight.set(int(temp[7]))
            checkGuess.set(int(temp[8]))
            
            if (checkLim.get() == 1 or checkWeight.get() == 1 or checkGuess.get() == 1):
                menu6.entryconfig('Open settings window', state='active')
            
            freqLim = float(temp[9])
            weight = float(temp[10])
            init_a = float(temp[11])
            init_b = float(temp[12])
            init_m0 = float(temp[13])
            
            resMn_GEX.config(state='active')
            resMw_GEX.config(state='active')
            resMz_GEX.config(state='active')
            resMwMn_GEX.config(state='active')
            
            resMn_GEX.insert(0,str("{:.0f}".format(ufloat(temp[14],temp[15]))))
            resMw_GEX.insert(0,str("{:.0f}".format(ufloat(temp[16],temp[17]))))
            resMz_GEX.insert(0,str("{:.0f}".format(ufloat(temp[18],temp[19]))))
            resMwMn_GEX.insert(0,format(ufloat(temp[20],temp[21]), '.2f'))
            
            resMn_GEX.config(state='readonly')
            resMw_GEX.config(state='readonly')
            resMz_GEX.config(state='readonly')
            resMwMn_GEX.config(state='readonly')
            
            menu5.entryconfig('Copy to clipboard', state='active')
            
            a_par = ufloat(temp[22],temp[23])
            b_par = ufloat(temp[24],temp[25])
            m0_par = ufloat(temp[26],temp[27])
            
            omega = data_mod[:,0]
            G1 = data_mod[:,1]
            G2 = data_mod[:,2]
            
            menu4.entryconfig('Plot dynamic moduli', state='active')
            
            M_GEX = dist[:,0]
            mwd_GEX = dist[:,1]
            
            plot_mwd_GEX()

            menu4.entryconfig('Plot estimated distribution', state='active')

            mwd_flag_GEX = 1

            if gpc_flag_GEX == 1:
                menu4.entryconfig(
                    'Plot MWD + estimated distribution', state='active')
                
            report = 'a = '+str(a_par)+'\n'+'b = '+str(b_par)+'\n'+'m0 = '+str(m0_par)
            consoleWindow()
            
        f.close

"""------------------------------------------------------------------------------------------------        
###################################################################################################

                        CREATING A NEW MATERIAL ADD WINDOW
                    
################################################################################################"""


def newMaterial():

    def close_window2():
        root2.destroy()

    def get_values_newMat():

        try:
            str(mat_title.get())
            float(rept_a2.get())
            float(rept_B2.get())
            float(mat_k2.get())
            float(mat_M02.get())
            float(mat_G02.get())
            float(mat_Me2.get())

            mat_path = script_dir+"/materials/"

            data2w = rept_a2.get()+'\n'+rept_B2.get()+'\n'+mat_k2.get()+'\n' + \
                mat_M02.get()+'\n'+mat_G02.get()+'\n'+mat_Me2.get()

            if len(textBox.get(1.0, 'end')) > 1:
                data2w += '\n'+' '+textBox.get(1.0, 'end')+' '

            file = open(os.path.join(mat_path, str(mat_title.get())), "w+")
            file.write(data2w)

            file.close()

            options1 = get_options()
            options2 = get_options()
            
            selectMat1['values'] = [options1[m] for m in range(len(options1))]
            selectMat2['values'] = [options2[m] for m in range(len(options2))]

            showinfo(title="Done!",
                     message="A new material was added to your options.")

            mat_title.delete(0, 'end')
            rept_a2.delete(0, 'end')
            rept_a2.delete(0, 'end')
            rept_B2.delete(0, 'end')
            mat_k2.delete(0, 'end')
            mat_M02.delete(0, 'end')
            mat_G02.delete(0, 'end')
            mat_Me2.delete(0, 'end')
            textBox.delete(1.0, 'end')

        except ValueError:
            showerror(title="Invalid values",
                      message="Please check if the typed values are valid!")
            # close_window2()

    def create_controls_frame2(container):
        frame = ttk.Frame(container)

        save_button = ttk.Button(
            frame, text='  Save  ', command=get_values_newMat)
        save_button.grid(column=0, row=0)

        cancel_button = ttk.Button(
            frame, text=' Cancel ', command=close_window2)
        cancel_button.grid(column=1, row=0)

        for widget in frame.winfo_children():
            widget.grid(padx=5, pady=10)

        return frame

    root2 = tk.Toplevel(root)
    root2.wm_transient(root)

    #root2.attributes('-topmost', 'true')

    root2.title('Create new material parameters')

    window_width = 470
    window_height = 440

    # get the screen dimension
    screen_width = root2.winfo_screenwidth()
    screen_height = root2.winfo_screenheight()

    # find the center point
    center_x = int(screen_width/2 - window_width / 2)
    center_y = int(screen_height/2 - window_height / 2)

    # set the position of the window to the center of the screen
    root2.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')

    #root2.geometry("%dx%d" % (screen_width-2, screen_height-10))

    root2.resizable(False, False)

    label_title = ttk.Label(root2, text=" Material name ")
    label_title.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    mat_title = ttk.Entry(root2, width=45)
    mat_title.pack(pady=5, padx=20, anchor='w', fill='x')

    # ------------------------------------------------------------------------------------

    label_a2 = ttk.Label(root2, text=" a: Power law index ")
    label_a2.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    rept_a2 = ttk.Entry(root2, width=45)
    rept_a2.pack(pady=5, padx=20, anchor='w', fill='x')

    # ------------------------------------------------------------------------------------

    label_B2 = ttk.Label(root2, text=u" \u03B2: Reptation parameter ")
    label_B2.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    rept_B2 = ttk.Entry(root2)
    rept_B2.pack(pady=5, padx=20, anchor='w', fill='x')

    # ------------------------------------------------------------------------------------

    label_k2 = ttk.Label(
        root2, text=" k: Molecular weight dependent time parameter [s.(mol/g)^a] ")
    label_k2.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    mat_k2 = ttk.Entry(root2)
    mat_k2.pack(pady=5, padx=20, anchor='w', fill='x')

    # ------------------------------------------------------------------------------------

    label_M02 = ttk.Label(root2, text=" M0: Monomer molecular weight [g/mol] ")
    label_M02.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    mat_M02 = ttk.Entry(root2, width=45)
    mat_M02.pack(pady=5, padx=20, anchor='w', fill='x')

    # ------------------------------------------------------------------------------------

    label_G02 = ttk.Label(root2, text=" G0: Plateau modulus [Pa] ")
    label_G02.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    mat_G02 = ttk.Entry(root2, width=45)
    mat_G02.pack(pady=5, padx=20, anchor='w', fill='x')

    # ------------------------------------------------------------------------------------

    label_Me2 = ttk.Label(
        root2, text=" Me: Molecular weight between entanglements [g/mol] ")
    label_Me2.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    mat_Me2 = ttk.Entry(root2, width=45)
    mat_Me2.pack(pady=5, padx=20, anchor='w', fill='x')

    label_text = ttk.Label(root2, text=" Notes: ")
    label_text.pack(pady=(5, 0), padx=20, fill='x', anchor='w')

    textBox = tk.Text(root2, height=5)
    textBox.pack(pady=5, padx=20, fill='x')

    controls_frame2 = create_controls_frame2(root2)
    controls_frame2.pack(padx=10)

    root2.mainloop()
    
"""------------------------------------------------------------------------------------------------        
###################################################################################################

                    CREATION OF THE FRAME OF THE UPPER BUTTONS MENU: GAR
                    
################################################################################################"""


def create_button_frame1(container, img1, img2, img3, img4, img5, img6, img7):

    frame = ttk.Frame(container)

    open_button = ttk.Button(frame, text=' Open file ',
                             command=open_log, compound=tk.LEFT, image=img1)
    open_button.grid(column=0, row=0)
    Hovertip(open_button, 'Open molecular weigth distribution log file',
             hover_delay=delay_value)

    saveAs_button = ttk.Button(
        frame, text=' Save ', command=save_log, compound=tk.LEFT, image=img2)
    saveAs_button.grid(column=1, row=0)
    Hovertip(saveAs_button, 'Save molecular weigth distribution log file',
             hover_delay=delay_value)

    material_button = ttk.Button(frame, text=' New material ', command=newMaterial,
                                 compound=tk.LEFT, image=img3)
    material_button.grid(column=3, row=0)
    Hovertip(material_button, ' Add new material parameters option',
             hover_delay=delay_value)

    spectrum_button = ttk.Button(frame, text=' Open RTS ', command=open_spectrum,
                                 compound=tk.LEFT, image=img4)
    spectrum_button.grid(column=4, row=0)
    Hovertip(spectrum_button, 'Open relaxation time spectrum file',
             hover_delay=delay_value)

    gpc_button = ttk.Button(frame, text=' Open MWD ', command=open_gpc,
                            compound=tk.LEFT, image=img5)
    gpc_button.grid(column=5, row=0)
    Hovertip(gpc_button, 'Open molecular weight distribution file', hover_delay=delay_value)

    help_button = ttk.Button(frame, text=' Help ', command=create_help, compound=tk.LEFT, image=img6)
    help_button.grid(column=6, row=0)
    Hovertip(help_button, 'Help', hover_delay=delay_value)

    info_button = ttk.Button(frame, text='About', command=create_about, compound=tk.LEFT, image=img7)
    info_button.grid(column=7, row=0)
    Hovertip(info_button, 'About', hover_delay=delay_value)

    for widget in frame.winfo_children():
        widget.grid(padx=5, pady=10)

    return frame

"""------------------------------------------------------------------------------------------------        
###################################################################################################

                    CREATION OF THE FRAME OF THE UPPER BUTTONS MENU: GEX
                    
################################################################################################"""


def create_button_frame2(container, img1, img2, img3, img4, img5, img6, img7):

    frame = ttk.Frame(container)

    open_button = ttk.Button(frame, text=' Open file ',
                             command=open_log_GEX, compound=tk.LEFT, image=img1)
    open_button.grid(column=0, row=0)
    Hovertip(open_button, 'Open molecular weigth distribution log file',
             hover_delay=delay_value)

    saveAs_button = ttk.Button(
        frame, text=' Save ', command=save_log_GEX, compound=tk.LEFT, image=img2)
    saveAs_button.grid(column=1, row=0)
    Hovertip(saveAs_button, 'Save molecular weigth distribution log file',
             hover_delay=delay_value)

    material_button = ttk.Button(frame, text=' New material ', command=newMaterial,
                                 compound=tk.LEFT, image=img3)
    material_button.grid(column=3, row=0)
    Hovertip(material_button, ' Add new material parameters option',
             hover_delay=delay_value)

    mod_button = ttk.Button(frame, text=' Open G*(\u03C9) ', command=open_mod,
                                 compound=tk.LEFT, image=img4)
    mod_button.grid(column=4, row=0)
    Hovertip(mod_button, 'Open dynamic moduli (G` and G``) file',
             hover_delay=delay_value)

    gpc_button = ttk.Button(frame, text=' Open MWD ', command=open_gpc_GEX,
                            compound=tk.LEFT, image=img5)
    gpc_button.grid(column=5, row=0)
    Hovertip(gpc_button, 'Open molecular weight distribution file', hover_delay=delay_value)

    help_button = ttk.Button(frame, text=' Help ', command=create_help, compound=tk.LEFT, image=img6)
    help_button.grid(column=6, row=0)
    Hovertip(help_button, 'Help', hover_delay=delay_value)

    info_button = ttk.Button(frame, text='About', command=create_about, compound=tk.LEFT, image=img7)
    info_button.grid(column=7, row=0)
    Hovertip(info_button, 'About', hover_delay=delay_value)

    for widget in frame.winfo_children():
        widget.grid(padx=4, pady=10)

    return frame



"""------------------------------------------------------------------------------------------------
###################################################################################################

                    CREATION OF THE PROGRAM COMMAND BUTTONS FRAME: GAR
                    
################################################################################################"""

def create_controls_frame_GAR(container, img8, img9, img10):
    frame = ttk.Frame(container)

    run_button = ttk.Button(frame, text='  RUN   ', image=img8, command=run_like_hell,
                            compound=tk.LEFT)
    run_button.grid(column=0, row=0)
    Hovertip(run_button, ' Run the program with the current parameters',
             hover_delay=delay_value)

    clear_button = ttk.Button(frame, text=' CLEAR  ', image=img9, command=confirm_clear_GAR,
                              compound=tk.LEFT)
    clear_button.grid(column=1, row=0)
    Hovertip(clear_button, ' Clear all the current parameters and plots',
             hover_delay=delay_value)

    export_button = ttk.Button(frame, text=' EXPORT ', image=img10, command=export_mwd,
                               compound=tk.LEFT)
    export_button.grid(column=2, row=0)
    Hovertip(export_button, ' Export the MWD data to a text file',
             hover_delay=delay_value)

    for widget in frame.winfo_children():
        widget.grid(padx=5, pady=10)

    return frame

"""------------------------------------------------------------------------------------------------
###################################################################################################

                    CREATION OF THE PROGRAM COMMAND BUTTONS FRAME: GEX
                    
################################################################################################"""

def create_controls_frame_GEX(container, img8, img9, img10, img11):
    frame = ttk.Frame(container)

    run_button = ttk.Button(frame, text='  RUN   ', image=img8, command=run_for_your_life,
                            compound=tk.LEFT)
    run_button.grid(column=0, row=0)
    Hovertip(run_button, ' Run the program with the current parameters',
             hover_delay=delay_value)

    clear_button = ttk.Button(frame, text=' CLEAR  ', image=img9, command=confirm_clear_GEX,
                              compound=tk.LEFT)
    clear_button.grid(column=1, row=0)
    Hovertip(clear_button, ' Clear all the current parameters and plots',
             hover_delay=delay_value)

    export_button = ttk.Button(frame, text=' EXPORT ', image=img10, command=export_mwd_GEX,
                               compound=tk.LEFT)
    export_button.grid(column=2, row=0)
    Hovertip(export_button, ' Export the MWD data to a text file',
             hover_delay=delay_value)
    
    console_button = ttk.Button(frame, text=' CONSOLE ', image=img11, command=consoleWindow,
                               compound=tk.LEFT)
    console_button.grid(column=3, row=0)
    Hovertip(export_button, 'Open console',
             hover_delay=delay_value)

    for widget in frame.winfo_children():
        widget.grid(padx=5, pady=10)

    return frame


"""################################################################################################

                            CREATION OF TABS FOR DIFFERENT MODELS
                            
################################################################################################"""

mwd_notebook = ttk.Notebook(root)
mwd_notebook.place(x=0, y=10)

frame1 = ttk.Frame(mwd_notebook, width=screen_width, height=screen_height)
frame1.place()
frame2 = ttk.Frame(mwd_notebook, width=screen_width, height=screen_height)
frame2.place()
#frame3 = ttk.Frame(mwd_notebook, width=screen_width, height=screen_height)
# frame3.place()
#frame4 = ttk.Frame(mwd_notebook, width=screen_width, height=screen_height)
#frame4.place()

mwd_notebook.add(frame1, text=" Generalized Analytical Relation ")
mwd_notebook.add(frame2, text=" Generalized Exponential Distribution ")
#mwd_notebook.add(frame3,text=" Another Model ")
#mwd_notebook.add(frame4, text=" Models Comparison ")

#mwd_notebook.tab(2, state="disabled")

button_frame1 = create_button_frame1(frame1, open_icon, saveAs_icon, material_icon,
                                   spectrum_icon, gpc_icon, help_icon, info_icon)

button_frame1.place(x=15, y=10)

button_frame2 = create_button_frame2(frame2, open_icon, saveAs_icon, material_icon,
                                   spectrum_icon, gpc_icon, help_icon, info_icon)

button_frame2.place(x=15, y=10)

controls_frame_GAR = create_controls_frame_GAR(
    frame1, run_icon, clear_icon, export_icon)
controls_frame_GAR.place(x=560, y=525)

controls_frame_GEX = create_controls_frame_GEX(
    frame2, run_icon, clear_icon, export_icon, console_icon)
controls_frame_GEX.place(x=440, y=525)

"""------------------------------------------------------------------------------------------------
###################################################################################################

CREATION OF THE FRAME OF THE PARAMETERS OF THE MATERIAL AND THE REPTATION MODEL: GAR
        
################################################################################################"""


reptFrame1 = ttk.LabelFrame(
    frame1, text=' Reptation model parameters ', height=210, width=250)
reptFrame1.place(x=20, y=85)
Hovertip(reptFrame1, 'For more info about this section, check the "Help" option on the top menu',
         hover_delay=delay_value)

selected_value1 = tk.StringVar()
selectMat1 = ttk.Combobox(frame1, text='Select material',
                         width=29, textvariable=selected_value1)
selectMat1.set("Select a material")
options1 = get_options()
selectMat1['values'] = [options1[m] for m in range(len(options1))]
#selectMat['state'] = 'readonly'

selectMat1.place(x=35, y=110)
Hovertip(selectMat1, 'Select pre-saved material and reptation parameters',
         hover_delay=delay_value)

#opt = selected_value.get()
selectMat1.bind('<<ComboboxSelected>>', opt_changed1)

rept_a_GAR = ttk.Entry(frame1, width=15)
rept_a_GAR.place(x=60, y=135)
label_a_GAR = ttk.Label(frame1, text='a:')
Hovertip(label_a_GAR, 'Power law index', hover_delay=delay_value)
label_a_GAR.place(x=35, y=135)
label_un0_GAR = ttk.Label(frame1, text='----')
label_un0_GAR.place(x=175, y=135)

rept_B_GAR = ttk.Entry(frame1, width=15)
rept_B_GAR.place(x=60, y=160)
label_B_GAR = ttk.Label(frame1, text=u"\u03B2:")
Hovertip(label_B_GAR, 'Reptation parameter', hover_delay=delay_value)
label_B_GAR.place(x=30, y=160)
label_un00_GAR = ttk.Label(frame1, text='----')
label_un00_GAR.place(x=175, y=160)

mat_k_GAR = ttk.Entry(frame1, width=15)
mat_k_GAR.place(x=60, y=185)
label_k_GAR = ttk.Label(frame1, text='k:')
Hovertip(label_k_GAR, 'Molecular weight dependent time parameter',
         hover_delay=delay_value)
label_k_GAR.place(x=35, y=185)
label_un1_GAR = ttk.Label(frame1, text='s.(mol/g)^a')
label_un1_GAR.place(x=175, y=185)

mat_M0_GAR = ttk.Entry(frame1, width=15)
mat_M0_GAR.place(x=60, y=210)
label_M0_GAR = ttk.Label(frame1, text='M0:')
Hovertip(label_M0_GAR, 'Monomer molecular weight', hover_delay=delay_value)
label_M0_GAR.place(x=35, y=210)
label_un2_GAR = ttk.Label(frame1, text='g/mol')
label_un2_GAR.place(x=175, y=210)

mat_G0_GAR = ttk.Entry(frame1, width=15)
mat_G0_GAR.place(x=60, y=235)
label_G0_GAR = ttk.Label(frame1, text='G0:')
Hovertip(label_G0_GAR, 'Plateau modulus', hover_delay=delay_value)
label_G0_GAR.place(x=35, y=235)
label_un3_GAR = ttk.Label(frame1, text='Pa')
label_un3_GAR.place(x=175, y=235)

mat_Me_GAR = ttk.Entry(frame1, width=15)
mat_Me_GAR.place(x=60, y=260)
label_Me_GAR = ttk.Label(frame1, text='Me:')
Hovertip(label_Me_GAR, 'Molecular weight between entanglements',
         hover_delay=delay_value)
label_Me_GAR.place(x=35, y=260)
label_un4_GAR = ttk.Label(frame1, text='g/mol')
label_un4_GAR.place(x=175, y=260)

"""------------------------------------------------------------------------------------------------
###################################################################################################

CREATION OF THE FRAME OF THE PARAMETERS OF THE MATERIAL AND THE REPTATION MODEL: GEX
        
################################################################################################"""

reptFrame2 = ttk.LabelFrame(
    frame2, text=' Reptation model parameters ', height=210, width=250)
reptFrame2.place(x=20, y=85)
Hovertip(reptFrame2, 'For more info about this section, check the "Help" option on the top menu',
         hover_delay=delay_value)

selected_value2 = tk.StringVar()
selectMat2 = ttk.Combobox(frame2, text='Select material',
                         width=29, textvariable=selected_value2)
selectMat2.set("Select a material")
options2 = get_options()
selectMat2['values'] = [options2[m] for m in range(len(options1))]
#selectMat['state'] = 'readonly'

selectMat2.place(x=35, y=110)
Hovertip(selectMat2, 'Select pre-saved material and reptation parameters',
         hover_delay=delay_value)

#opt = selected_value.get()
selectMat2.bind('<<ComboboxSelected>>', opt_changed2)

rept_a_GEX = ttk.Entry(frame2, width=15)
rept_a_GEX.place(x=60, y=135)
label_a_GEX = ttk.Label(frame2, text='a:')
Hovertip(label_a_GEX, 'Power law index', hover_delay=delay_value)
label_a_GEX.place(x=35, y=135)
label_un0_GEX = ttk.Label(frame2, text='----')
label_un0_GEX.place(x=175, y=135)

rept_B_GEX = ttk.Entry(frame2, width=15)
rept_B_GEX.place(x=60, y=160)
label_B_GEX = ttk.Label(frame2, text=u"\u03B2:")
Hovertip(label_B_GEX, 'Reptation parameter', hover_delay=delay_value)
label_B_GEX.place(x=30, y=160)
label_un00_GEX = ttk.Label(frame2, text='----')
label_un00_GEX.place(x=175, y=160)

mat_k_GEX = ttk.Entry(frame2, width=15)
mat_k_GEX.place(x=60, y=185)
label_k_GEX = ttk.Label(frame2, text='k:')
Hovertip(label_k_GEX, 'Molecular weight dependent time parameter',
         hover_delay=delay_value)
label_k_GEX.place(x=35, y=185)
label_un1_GEX = ttk.Label(frame2, text='s.(mol/g)^a')
label_un1_GEX.place(x=175, y=185)

mat_M0_GEX = ttk.Entry(frame2, width=15)
mat_M0_GEX.place(x=60, y=210)
label_M0_GEX = ttk.Label(frame2, text='M0:')
Hovertip(label_M0_GEX, 'Monomer molecular weight', hover_delay=delay_value)
label_M0_GEX.place(x=35, y=210)
label_un2_GEX = ttk.Label(frame2, text='g/mol')
label_un2_GEX.place(x=175, y=210)

mat_G0_GEX = ttk.Entry(frame2, width=15)
mat_G0_GEX.place(x=60, y=235)
label_G0_GEX = ttk.Label(frame2, text='G0:')
Hovertip(label_G0_GEX, 'Plateau modulus', hover_delay=delay_value)
label_G0_GEX.place(x=35, y=235)
label_un3_GEX = ttk.Label(frame2, text='Pa')
label_un3_GEX.place(x=175, y=235)

mat_Me_GEX = ttk.Entry(frame2, width=15)
mat_Me_GEX.place(x=60, y=260)
label_Me_GEX = ttk.Label(frame2, text='Me:')
Hovertip(label_Me_GEX, 'Molecular weight between entanglements',
         hover_delay=delay_value)
label_Me_GEX.place(x=35, y=260)
label_un4_GEX = ttk.Label(frame2, text='g/mol')
label_un4_GEX.place(x=175, y=260)

"""------------------------------------------------------------------------------------------------
###################################################################################################

            CREATING THE ADJUSTMENT FRAME OF ROUSE MODES AND OTHER SETTINGS
            
################################################################################################"""

adjustFrame = ttk.LabelFrame(
    frame1, text=' Settings ', height=110, width=250)
adjustFrame.place(x=20, y=315)
Hovertip(adjustFrame, 'For more info about this section, check the "Help" option on the top menu',
         hover_delay=delay_value)


checkMe1 = tk.IntVar()
checkRouseVal = tk.IntVar()
checkMe2 = tk.IntVar()

checkMe_spec = ttk.Checkbutton(frame1, text=' Subtract Me ', variable=checkMe1,
                               onvalue=1, offvalue=0)
checkMe_spec.place(x=35, y=345)
checkMe1.set(1)

checkRouse = ttk.Checkbutton(frame1, text=' Subtract Rouse modes', variable=checkRouseVal,
                             onvalue=1, offvalue=0, command=configRouse)
checkRouse.place(x=35, y=385)
checkRouseVal.set(0)

#checkMe_dist = ttk.Checkbutton(frame1, text=' Add Me on the estimated MWD', variable=checkMe2,
#                               onvalue=1, offvalue=0)
#checkMe_dist.place(x=35, y=390)
#checkMe2.set(1)

"""------------------------------------------------------------------------------------------------
###################################################################################################

                CREATING THE ADJUSTMENT FRAME AND OTHER SETTINGS: GEX
            
################################################################################################"""

configFrame = ttk.LabelFrame(
    frame2, text=' Settings ', height=110, width=250)
configFrame.place(x=20, y=315)
Hovertip(configFrame, 'For more info about this section, check the "Help" option on the top menu',
         hover_delay=delay_value)


checkLim = tk.IntVar()
checkWeight = tk.IntVar()
checkGuess = tk.IntVar()

checkLim_mod = ttk.Checkbutton(frame2, text=' Change frequency window', variable=checkLim,
                               onvalue=1, offvalue=0, command=closing_configLim)
checkLim_mod.place(x=35, y=340)
checkLim.set(0)


checkWeight_mod = ttk.Checkbutton(frame2, text=' Change relative weight', variable=checkWeight,
                             onvalue=1, offvalue=0, command=closing_configWeight)
checkWeight_mod.place(x=35, y=365)
checkWeight.set(0)


checkGuess_adj = ttk.Checkbutton(frame2, text=' Change initial guesses', variable=checkGuess,
                               onvalue=1, offvalue=0, command=closing_configGuess)
checkGuess_adj.place(x=35, y=390)
checkGuess.set(0)


"""------------------------------------------------------------------------------------------------
###################################################################################################

CREATION OF THE FRAME OF NUMERICAL RESULTS (AVERAGES AND POLYDISPERSION): GAR
            
################################################################################################"""

resultsFrame = ttk.LabelFrame(frame1, text=' Average MW\'s [g/mol] and ratio  ', height=135, width=250)
resultsFrame.place(x=20, y=445)
Hovertip(resultsFrame, 'For more info about this section, check the "Help" option on the top menu',
         hover_delay=delay_value)


resMn = ttk.Entry(frame1, width=22, state='readonly')
#resMn = ttk.Label(frame1, text='                ')
resMn.place(x=85, y=470)
label_Mn = ttk.Label(frame1, text=' Mn:')
label_Mn.place(x=35, y=470)
Hovertip(label_Mn, 'Number average molecular weight', hover_delay=delay_value)
#label_un5 = ttk.Label(frame1, text='')
#abel_un5.place(x=210, y=470)


resMw = ttk.Entry(frame1, width=22, state='readonly')
#resMw = ttk.Label(frame1, text='                ')
resMw.place(x=85, y=495)
label_Mw = ttk.Label(frame1, text=' Mw:')
label_Mw.place(x=35, y=495)
Hovertip(label_Mw, 'Weight average molecular weight', hover_delay=delay_value)
#label_un6 = ttk.Label(frame1, text='')
#label_un6.place(x=210, y=495)


resMz = ttk.Entry(frame1, width=22, state='readonly')
#resMz = ttk.Label(frame1, text='                ')
resMz.place(x=85, y=520)
label_Mz = ttk.Label(frame1, text=' Mz:')
label_Mz.place(x=35, y=520)
Hovertip(label_Mz, 'z-Average molecular weight', hover_delay=delay_value)
#label_un7 = ttk.Label(frame1, text='')
#label_un7.place(x=210, y=520)


resMwMn = ttk.Entry(frame1, width=22, state='readonly')
#resMwMn = ttk.Label(frame1, text='                ')
resMwMn.place(x=85, y=545)
label_MwMn = ttk.Label(frame1, text='Mw/Mn:')
label_MwMn.place(x=35, y=545)
Hovertip(label_MwMn, 'Polydispersity index', hover_delay=delay_value)

"""------------------------------------------------------------------------------------------------
###################################################################################################

CREATION OF THE FRAME OF NUMERICAL RESULTS (AVERAGES AND POLYDISPERSION): GEX
            
################################################################################################"""

resultsFrame_GEX = ttk.LabelFrame(frame2, text=' Average MW\'s [g/mol] and ratio ', height=135, width=250)
resultsFrame_GEX.place(x=20, y=445)
Hovertip(resultsFrame_GEX, 'For more info about this section, check the "Help" option on the top menu',
         hover_delay=delay_value)

resMn_GEX = ttk.Entry(frame2, width=22, state='readonly')
resMn_GEX.place(x=85, y=470)
label_Mn_GEX = ttk.Label(frame2, text=' Mn:')
label_Mn_GEX.place(x=35, y=470)
Hovertip(label_Mn_GEX, 'Number average molecular weight', hover_delay=delay_value)

resMw_GEX = ttk.Entry(frame2, width=22, state='readonly')
resMw_GEX.place(x=85, y=495)
label_Mw_GEX = ttk.Label(frame2, text=' Mw:')
label_Mw_GEX.place(x=35, y=495)
Hovertip(label_Mw_GEX, 'Weight average molecular weight', hover_delay=delay_value)

resMz_GEX = ttk.Entry(frame2, width=22, state='readonly')
resMz_GEX.place(x=85, y=520)
label_Mz_GEX = ttk.Label(frame2, text=' Mz:')
label_Mz_GEX.place(x=35, y=520)
Hovertip(label_Mz_GEX, 'z-Average molecular weight', hover_delay=delay_value)

resMwMn_GEX = ttk.Entry(frame2, width=22, state='readonly')
resMwMn_GEX.place(x=85, y=545)
label_MwMn_GEX = ttk.Label(frame2, text='Mw/Mn:')
label_MwMn_GEX.place(x=35, y=545)
Hovertip(label_MwMn_GEX, 'Polydispersity index', hover_delay=delay_value)

"""------------------------------------------------------------------------------------------------
###################################################################################################

                                CREATION OF THE GRAPHICS FRAME
                                
################################################################################################"""

plotsFrame1 = ttk.LabelFrame(frame1, text='Graph', width=615, height=440)
plotsFrame1.place(x=300, y=85)
Hovertip(plotsFrame1, 'For more info about this section, check the "Help" option on the top menu')

plotsFrame2 = ttk.LabelFrame(frame2, text='Graph', width=615, height=440)
plotsFrame2.place(x=300, y=85)
Hovertip(plotsFrame2, 'For more info about this section, check the "Help" option on the top menu')

"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                    ADDITIONAL FUNCTIONS: KEYBOARD COMMANDS AND POPUP MENU
                        
################################################################################################"""

def run_key(e=None):
    run_like_hell()

def clear_key_GEX(e=None):
    confirm_clear_GEX()
    
def clear_key_GAR(e=None):
    confirm_clear_GAR()

#def clear_all_key(e=None):
#    clear_all_GAR()

def close_window(e=None):
    root.destroy()

def export_key(e=None):
    export_mwd()
    
def export_mod(e=None):
    export_moduli()

def my_popup(e):
    menu1.tk_popup(e.x_root, e.y_root)


menu1 = tk.Menu(root, tearoff=False)

menu1.add_command(label='Plot estimated distribution',
                  command=plot_mwd, state='disabled')
menu1.add_command(label='Plot MWD data only',
                  command=plot_gpc, state='disabled')
menu1.add_command(label='Plot MWD + estimated distribution',
                  command=plot_both, state='disabled')
menu1.add_command(label='Plot relaxation time spectrum',
                  command=plot_spec, state='disabled')



def my_popup2(e):
    menu2.tk_popup(e.x_root, e.y_root)

menu2 = tk.Menu(root, tearoff=False)
menu2.add_command(label='Copy to clipboard',command=copy_values, state='disabled')

def my_popup3(e):
    menu3.tk_popup(e.x_root, e.y_root)

menu3 = tk.Menu(root, tearoff=False)
menu3.add_command(label='Open Rouse settings',command=configRouse, state='disabled')


def my_popup4(e):
    menu4.tk_popup(e.x_root, e.y_root)
    
menu4 = tk.Menu(root, tearoff=False)

menu4.add_command(label='Plot estimated distribution',
                  command=plot_mwd_GEX, state='disabled')
menu4.add_command(label='Plot MWD data only',
                  command=plot_gpc_GEX, state='disabled')
menu4.add_command(label='Plot MWD + estimated distribution',
                  command=plot_both_GEX, state='disabled')
menu4.add_command(label='Plot dynamic moduli',
                  command=plot_mod, state='disabled')
menu4.add_command(label='Plot dynamic moduli + fitted moduli',
                  command=plot_mod_GEX, state='disabled')

def my_popup5(e):
    menu5.tk_popup(e.x_root, e.y_root)

menu5 = tk.Menu(root, tearoff=False)
menu5.add_command(label='Copy to clipboard',command=copy_values_GEX, state='disabled')

def my_popup6(e):
    menu6.tk_popup(e.x_root, e.y_root)

menu6 = tk.Menu(root, tearoff=False)
menu6.add_command(label='Open settings window',command=configGEX, state='disabled')

def on_closing():
    if mwd_flag == 1 or mwd_flag_GEX == 1: 
        if askokcancel("Quit", "Do you really want to quit? All your current analysis will be lost."):
            root.destroy()
    else:
        root.destroy()
            
# Chamada da função "clear_all()" na inicialização do programa
clear_all_GAR()
clear_all_GEX()

plotsFrame1.bind("<Button-3>", my_popup)
resultsFrame.bind("<Button-3>", my_popup2)
adjustFrame.bind("<Button-3>", my_popup3)
plotsFrame2.bind("<Button-3>", my_popup4)
resultsFrame_GEX.bind("<Button-3>", my_popup5)
configFrame.bind("<Button-3>", my_popup6)

#root.bind('<Control-n>', newMat_key)

#root.bind('<Control-r>', run_key)
#root.bind('<Control-l>', clear_key)
#root.bind('<Control-e>', export_key)

#root.bind('<Control-L>', clear_all_key)

root.bind('<Control-G>', export_mod)

#root.bind('<Escape>', close_window)


"""------------------------------------------------------------------------------------------------ 
###################################################################################################

                                 THIS IS THE END OF THE LINE

################################################################################################"""

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()
