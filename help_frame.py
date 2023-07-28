#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:22:49 2023

@author: atilio
"""

from tkinter import ttk
import tkinter as tk

def create_help():
    
    root = tk.Tk()
    root.title('Help!')
    
    window_width = 750
    window_height = 500
    
    # get the screen dimension
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    
    # find the center point
    center_x = int(screen_width/2 - window_width / 2)
    center_y = int(screen_height/2 - window_height / 2)
    
    # set the position of the window to the center of the screen
    root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
    
    #root.geometry("%dx%d" % (screen_width-2, screen_height-10))
    
    #root.resizable(False, False)
    
    main_frame = ttk.Frame(root)
    main_frame.pack(fill=tk.BOTH, expand=1)
    
    my_canvas = tk.Canvas(main_frame)
    my_canvas.pack(side='left',fill=tk.BOTH, expand=1)
    
    my_scrollbar = ttk.Scrollbar(main_frame, orient='vertical', command=my_canvas.yview)
    my_scrollbar.pack(side='right',fill=tk.Y,padx=10,pady=10)
    
    my_canvas.configure(yscrollcommand=my_scrollbar.set)
    my_canvas.bind('<Configure>', lambda e: my_canvas.configure(scrollregion=my_canvas.bbox("all")))
    
    second_frame = ttk.Frame(my_canvas)
    
    my_canvas.create_window((0,0),window=second_frame,anchor=tk.NW)
    
    intro = ttk.Label(second_frame, text='Hello there! \nBelow you will find a set of general instructions and information concerned to the software functioning. For more specific details, please check the recommended papers or feel free to contact us. \nCheers!',\
                           font=("Arial", 11), width=150,wraplength=650)
    intro.pack(padx=20, pady=10, anchor=tk.NW)
    
    about_title1 = ttk.Label(second_frame, 
                           text='>>> About the software',font=("Arial", 11, "bold"))
    about_title1.pack(padx=20, pady=(5,0),anchor=tk.NW)
    
    about_text1 = ttk.Label(second_frame, 
                           text='PolyWeight is a software dedicated to molecular weight determination of linear polymers. By utilizing dynamic moduli, users can calculate MWD as well as molecular weight averages such as Mn, Mw, and Mz with two solving methods: one method based on an analytical relation and one way based on a parametric model.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text1.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_title2 = ttk.Label(second_frame, 
                           text='>>> About the solution methods',font=("Arial", 11, "bold"))
    about_title2.pack(padx=20, pady=(5,0),anchor=tk.NW)
    
    about_solution1 = ttk.Label(second_frame, 
                           text='\u2022 Generalized Analytical Relation: Obtaining the MWD through analytical relationships is based on closed mathematical forms, with which it is possible to calculate the distribution without using recurrent methods or optimization processes. The generalized analytical relation incorporates a mathematical formulation that considers the relaxation spectrum, which is obtained from the dynamic moduli using the NLREG software. The numerical integration of the relaxation spectrum is performed utilizing the integrate.quad function from the SciPy library. ',
                               font=("Arial", 11), width=150,wraplength=650)
    about_solution1.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_solution2 = ttk.Label(second_frame, 
                           text='\u2022 Generalized Exponential Distribution: An optimization process is employed to perform a multiobjective fit of the viscoelastic models to the dynamic moduli, where this models are expressed as functions of the generalized exponential (GEX) distribution parameters The integrals in this case are also calculated using the integrate.quad package from the SciPy library. The optimization procedure utilizes the lmfit library, where a Minimizer object is created based on the optimization function and subsequently utilized with the minimize function, employing the default Levenberg-Marquardt algorithm for the fit.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_solution2.pack(padx=20, pady=5,anchor=tk.NW)
    
    
    about_title3 = ttk.Label(second_frame, 
                           text='>>> About the interface',font=("Arial", 11, "bold"))
    about_title3.pack(padx=20, pady=(20,0),anchor=tk.NW)
    
    
    about_text2 = ttk.Label(second_frame, 
                           text='\u2022 Folders: The software is divided into two tabs that can be seen at the top of the interface named "Generalized Analytical Relation" and “Generalized Exponential Distribution". In each of these tabs, sets of configurations and fields referring to the different implemented methods can be found.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text2.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text3 = ttk.Label(second_frame, 
                           text='\u2022 Upper button menu: The button menu located at the top of the interface, just below the tabs, has 7 buttons: ',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text3.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text4 = ttk.Label(second_frame, 
                           text='\u2713 Open file: Allows the user to open a .txt file created with the “Save” function (see next item) and containing the data of a previous analysis made with PolyWeight. When clicking on this button, a message is displayed that informs the user that, if he opens a log file, all the information of the current analysis will be overwritten, so that the file will only be opened upon acceptance;',
                               font=("Arial", 11), width=150,wraplength=625)
    about_text4.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text5 = ttk.Label(second_frame, 
                           text='\u2713 Save: Creates a file with all of the current analysis data (the resulting molecular weight averages and ratio, estimated MWD, settings, etc.) which can be later opened with the “Open file” function. If there is no current analysis, the software will inform the user that there is no data available to be saved;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text5.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text6 = ttk.Label(second_frame, 
                           text='\u2713 New material: Opens a new window for the user to create a new set of parameters for a given polymer. This set can be saved and accessed with a dropdown menu located in the “Reptation model parameters” field, which is detailed below. All values must be numerical, with the exception of the name of the new material created and the notes at the end. Ex.: 1, 3.14, 13e-5, etc;',
                               font=("Arial", 11), width=150,wraplength=625)
    about_text6.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text7 = ttk.Label(second_frame, 
                           text='\u2713 Open RTS/Open G*(\u03C9): In the “Generalized Analytical Relation” tab, the user will find an “Open RTS” button, whichallows opening the relaxation time spectrum of a given polymer, while in the “Generalized Exponential Distribution” tab the user will find an “Open G*(\u03C9)”, which allows opening a file containing the dynamic moduli G*(\u03C9), i.e. the storage modulus, G′(\u03C9), and the loss modulus, G′′(\u03C9);',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text7.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text8 = ttk.Label(second_frame, 
                           text='\u2713 Open MWD: With this function, the user can open a file containing a MWD obtained with any method (e.g., GPC) to compare graphically with the current analysis;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text8.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text9 = ttk.Label(second_frame, 
                           text='\u2713 Help: Opens a window with general information about the software’s controls, settings, and instructions of use, along with some references;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text9.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text10 = ttk.Label(second_frame, 
                           text='\u2713 About: General and contact information about the software’s developers.',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text10.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text11 = ttk.Label(second_frame, 
                           text='\u2022 Reptation model parameters: The reptation model is the liaison between the rheology of a linear amorphous polymer and its molecular weight distribution. In this field, the user should define the values of the parameters used in the calculations. For more details, consult the sections "Using PolyWeight", the articles indicated in the section "Selected Bibliography" or the USER MANUAL.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text11.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text111 = ttk.Label(second_frame, 
                           text='\u2713 Dropdown menu: The user can choose a set of parameters previously created with the “New material” option;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text111.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text12 = ttk.Label(second_frame, 
                           text='\u2713 Power law index (a): Experimental parameter in the scaling law that determines the relaxation time of a polymer sample of average molecular weight M. Dimensionless constant;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text12.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text13 = ttk.Label(second_frame, 
                           text='\u2713 Reptation parameter (\u03B2): Parameter defined explicitly in the reptation model. Dimensionless constant;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text13.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text14 = ttk.Label(second_frame, 
                           text='\u2713 Molecular weight dependent time parameter (k): Experimental parameter in the scaling law that determines the relaxation time of a polymer sample of average molecular weight M. Constant in s.(mol/g)^a;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text14.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text15 = ttk.Label(second_frame, 
                           text='\u2713 Monomer molecular weight (M0): Molar weight of the polymer monomer used in the current analysis, Value in g/mol;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text15.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text16 = ttk.Label(second_frame, 
                           text='\u2713 Plateau modulus (G0): Value obtained in the region of the viscoelastic plateau displayed in the storage modulus, G′ (\u03C9), for which there is a local minimum in the loss factor, tan \u03B4. Value in Pa;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text16.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text17 = ttk.Label(second_frame, 
                           text='\u2713 Molecular weight between entanglements (Me): Average molecular weight between topological constraints, as defined in the tube model. This value depends on the polymer used in the analysis in question. Value in g/mol.',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text17.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text181 = ttk.Label(second_frame, 
                           text='\u2022 GAR settings: ',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text181.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text1811 = ttk.Label(second_frame, 
                           text='\u2713 Subtract Me: This option allows the user to shift the entire relaxation spectrum to the left by the factor \u03C4(Me), consequently, the entire resulting MWD will be shifted to the left by the factor Me;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text1811.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text1812 = ttk.Label(second_frame, 
                           text='\u2713 Subtract Rouse modes: Even if the time range below \u03C4(Me) is disregarded the Rouse modes still influence on the relaxation processes in the adjacent time range. These modes can be subtracted with this option which, when checked, opens a parameter definition window where the user can define the slope of the curve referring to the Rouse spectrum: the original theory assumes −0.5, but the user can input any value or let the software define a value according to the curve of the relaxation spectrum. The user can also define in what time range of the spectrum the subtraction will occur. The Rouse spectrum subtraction limit is given in multiples of Me;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text1812.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text182 = ttk.Label(second_frame, 
                           text='\u2022 GEX settings: When selecting one of the options, the settings window is opened, and the selected option is enabled.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text182.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text1821 = ttk.Label(second_frame, 
                           text='\u2713 Change frequency window: This option allows one to select only a part of the experimental dynamic moduli within a frequency window by defining a maximum limit frequency (i.e., the upper limit of the window) for which the GEX model will apply. If this option is unchecked, the complete frequency range of the experimental moduli is used;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text1821.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text1822 = ttk.Label(second_frame, 
                           text='\u2713 Change relative weight: The relative weight weighs the influence of the loss modulus fitting in relation to the storage modulus fitting in the objective function. The user can determine this value by analyzing the contribution of each portion of the optimization function. If this option is not checked, the software assumes that both fittings have the same weight;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text1822.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text1823 = ttk.Label(second_frame, 
                           text='\u2713 Change initial guesses: The GEX distribution parameters are initialized with default values (a = b = m0 = 1). This option enables the fields to change the initial guesses of the parameters and, consequently, the search in the solution space. In the case of the parameter m0, the initial default value and the value entered in its respective field in the settings window are multiplied by a factor 10^5;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text1823.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text19 = ttk.Label(second_frame, 
                           text='\u2022 Average MW’s [g/mol] and ratio: In this field are presented the numerical values of the number average molecular weight (Mn), the weight average molecular weight (Mw), the z-average molecular weight (Mz), and the polydispersity index, represented by the ratio Mw/Mn. When the calculated values are presented, the user can select them with the “Copy to clipboard” function, which appears by right-clicking the mouse within the average MW’s area;',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text19.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text20 = ttk.Label(second_frame, 
                           text='\u2022 Graph: A few graphs can be plotted, and the options are selectable via a pop-up menu that is opened by right-clicking the mouse within the graphics area. The options are:',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text20.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text201 = ttk.Label(second_frame, 
                           text='\u2713 Plot estimated distribution: Shows the MWD calculated with the respective method. When the calculation is finished, the distribution is plotted automatically, and if the visualization is switched to another dataset, it can be accessed again using this option;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text201.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text202 = ttk.Label(second_frame, 
                           text='\u2713 Plot MWD data only: Plots the molecular weight distribution graph of the data opened with the “Open MWD” option;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text202.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text203 = ttk.Label(second_frame, 
                           text='\u2713 Plot MWD + estimated distribution: Plots the distribution opened with the “Open MWD” option along with the calculated distribution on the same plane;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text203.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text204 = ttk.Label(second_frame, 
                           text='\u2713 Plot relaxation time spectrum/dynamic moduli: Plots the input datasets for the generalized analytical relation and the generalized exponential distribution, respectively;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text204.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text205 = ttk.Label(second_frame, 
                           text='\u2713 Plot dynamic moduli + fitted moduli: This option is available only for the generalized exponential distribution method, and it shows the dynamic moduli resulting from the fitting procedure made to the input dataset;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text205.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text21 = ttk.Label(second_frame, 
                           text='\u2022 Lower button menu: It has three common buttons with the functions of running the program (“RUN”), clearing all the variables and graphs (“CLEAR”), and exporting a file with the calculated MWD (“EXPORT”). The “Generalized exponential distribution” tab has an additional button, “CONSOLE”, which opens a new window with information and parameters of the fitting procedure made in this method.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text21.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_title4 = ttk.Label(second_frame, 
                           text='>>> Using PolyWeight',font=("Arial", 11, "bold"))
    about_title4.pack(padx=20, pady=(5,0),anchor=tk.NW)
    
    about_text31 = ttk.Label(second_frame, 
                           text='\u2022 Loading input datasets: The user must load a file containing the relaxation time spectrum in case of using the Generalized Analytical Relation method, or a file containing the dynamic mduli in case of using the Generalized Exponential Distribution method. In these files, the data must be arranged in 3 columns: relaxation time (in seconds), relaxation spectrum and error; frequency (in rad/s), G\' (in Pascal) and G\" (in Pascal). The files must be in one of the allowed formats (.txt, .dat or .sol) and the data must be organized by increasing frequency (from lowest to highest). If the user tries to run the software without a properly loaded file, the program will issue a warning;',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text31.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text32 = ttk.Label(second_frame, 
                           text='\u2022 Choosing set of reptation parameters: Values can be typed in specific fields, observing the units of each parameter. These choices must be made based on the user\'s previous knowledge based on the rheology of the polymer in question and on the specific literature. It is also possible to create a set of parameters with the “New material” option, which after being saved, becomes available for selection in the “Dropdown menu”. If there is a problem with the value of any parameter (for example, value not typed, wrong value format, etc) the program will issue a warning and will only be executed after correction;',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text32.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text33 = ttk.Label(second_frame, 
                           text='\u2022 Running the software: With the input file (relaxation spectrum or dynamic modules) properly loaded and with the values entered correctly, it is possible to run the software by clicking on the “RUN” button. It is important to point out that, when running the program (in either of the two methods), all interface controls are blocked and are released only at the end of the execution, when the molecular weight distribution graph is displayed. The execution time varies according to the user\'s hardware and the length of the input files;',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text33.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text34 = ttk.Label(second_frame, 
                           text='\u2022 Loading MWD file: The MWD file must be in one of the permitted formats (.txt, .dat or .sol), and the data must be arranged in two columns: molar mass (in g/mol) and, distribution. The data must be organized by increasing molar mass (from lowest to highest);',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text34.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text35 = ttk.Label(second_frame, 
                           text='\u2022 Clearing and exporting results: The user can clear all current analysis data, including input files, results, settings and graphs, by clicking on the “CLEAR” button. The user will be informed that, if he accepts to carry out the cleaning, all data from the current analysis will be lost. If no analysis has been performed, the program will show a message stating that there is no data to clean. The user can also export the MWD resulting from an analysis to a .txt file by clicking on the “EXPORT” button. Likewise, if no analysis has been carried out, the program will show a message stating that there is no data to be exported;',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text35.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text36 = ttk.Label(second_frame, 
                           text='\u2022 Saving and opening log files: The “Save” function creates a log file in .txt format with the values of the reptation parameters, the selected settings and their respective values, the average molecular weights obtained in the current analysis, the MWD resulting from the analysis and the input data set (relaxation spectrum or dynamic moduli). If no analysis has been carried out, the program will show a message stating that there is no data to be saved. These log files can later be opened with the “Open file” function, so that new analyses can be carried out on the data loaded in this way. When clicking on this button, the program informs that loading a log file will overwrite all the current analysis data, so this will only be carried out upon acceptance.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text36.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text37 = ttk.Label(second_frame, 
                           text='\u2022 Additional controls and options',
                               font=("Arial", 11), width=150,wraplength=650)
    about_text37.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_text371 = ttk.Label(second_frame, 
                           text='\u2713 Copy to clipboard: Average molecular weights can be copied to the clipboard once calculated. As the values are not selectable, this can be done using the “Copy to clipboard” option, which appears when right-clicking inside the numerical values area. This option is only activated when there are values shown and they can be inserted elsewhere with the shortcut Ctrl+V;',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text371.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_text372 = ttk.Label(second_frame, 
                           text='\u2713 Plots toolbar: An image toolbar appears below the graphics area. It can be used to interact (zoom in and out) with the graphics and, if the user wishes, to save them in .png format.',\
                               font=("Arial", 11), width=150,wraplength=625)
    about_text372.pack(padx=50, pady=5,anchor=tk.NW)
    
    about_title5 = ttk.Label(second_frame, 
                           text='>>> Selected Bibliography',font=("Arial", 11, "bold"))
    about_title5.pack(padx=20, pady=(5,0),anchor=tk.NW)
    
    about_ref1 = ttk.Label(second_frame, 
                           text='ANDERSSEN, R. S.; MEAD, D. W.; DRISCOLL IV, J. J. On the recovery of molecular weight functionals from the double reptation model. Journal of Non-Newtonian Fluid Mechanics, v. 68, n. 2-3, p. 291-301, 1997.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref1.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref2 = ttk.Label(second_frame, 
                           text='CARROT, Christian; GUILLET, Jacques. From dynamic moduli to molecular weight distribution: A study of various polydisperse linear polymers. Journal of Rheology, v. 41, n. 5, p. 1203-1220, 1997.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref2.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref3 = ttk.Label(second_frame, 
                           text='DEALY, John M.; READ, Daniel J.; LARSON, Ronald G. Structure and rheology of molten polymers: from structure to flow behavior and back again. Carl Hanser Verlag GmbH Co KG, 2018.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref3.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref4 = ttk.Label(second_frame, 
                           text='ELSTER, C.; HONERKAMP, J.; WEESE, J. Using regularization methods for the determination of relaxation and retardation spectra of polymeric liquids. Rheologica acta, v. 31, p. 161-174, 1992.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref4.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref5 = ttk.Label(second_frame, 
                           text='HONERKAMP, J.; WEESE, Jürgen. Determination of the relaxation spectrum by a regularization method. Macromolecules, v. 22, n. 11, p. 4372-4377, 1989.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref5.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref6 = ttk.Label(second_frame, 
                           text='HONERKAMP, Josef; WEESE, Jürgen. Tikhonovs regularization method for ill-posed problems: A comparison of different methods for the determination of the regularization parameter. Continuum Mechanics and Thermodynamics, v. 2, p. 17-30, 1990.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref6.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref7 = ttk.Label(second_frame, 
                           text='HONERKAMP, Josef; WEESE, Jurgen. A nonlinear regularization method for the calculation of relaxation spectra. Rheologica acta, v. 32, p. 65-73, 1993.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref7.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref8 = ttk.Label(second_frame, 
                           text='LÉONARDI, Frédéric; ALLAL, Ahmed; MARIN, Gérard. Molecular weight distribution from viscoelastic data: The importance of tube renewal and Rouse modes. Journal of Rheology, v. 46, n. 1, p. 209-224, 2002.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref8.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref9 = ttk.Label(second_frame, 
                           text='LIU, Chenyang et al. Evaluation of different methods for the determination of the plateau modulus and the entanglement molecular weight. Polymer, v. 47, n. 13, p. 4461-4479, 2006.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref9.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref10 = ttk.Label(second_frame, 
                           text='MAIER, D. et al. Evaluation of models combining rheological data with the molecular weight distribution. Journal of Rheology, v. 42, n. 5, p. 1153-1173, 1998.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref10.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref11 = ttk.Label(second_frame, 
                           text='MEAD, D. W. Determination of molecular weight distributions of linear flexible polymers from linear viscoelastic material functions. Journal of Rheology, v. 38, n. 6, p. 1797-1827, 1994.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref11.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref12 = ttk.Label(second_frame, 
                           text='NOBILE, Maria Rossella; COCCHINI, F.; LAWLER, John V. On the stability of molecular weight distributions as computed from the flow curves of polymer melts. Journal of Rheology, v. 40, n. 3, p. 363-382, 1996.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref12.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref13 = ttk.Label(second_frame, 
                           text='NOBILE, Maria Rosella; COCCHINI, Franco. Predictions of linear viscoelastic properties for polydisperse entangled polymers. Rheologica acta, v. 39, p. 152-162, 2000.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref13.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref14 = ttk.Label(second_frame, 
                           text='NOBILE, Maria Rossella; COCCHINI, Franco. Evaluation of molecular weight distribution from dynamic moduli. Rheologica acta, v. 40, p. 111-119, 2001.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref14.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref15 = ttk.Label(second_frame, 
                           text='NOBILE, Maria Rossella; COCCHINI, Franco. A generalized relation between MWD and relaxation time spectrum. Rheologica acta, v. 47, p. 509-519, 2008.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref15.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref16 = ttk.Label(second_frame, 
                           text='SHANBHAG, Sachin. Analytical rheology of polymer melts: State of the art. International Scholarly Research Notices, v. 2012, 2012.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref16.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref17 = ttk.Label(second_frame, 
                           text='THIMM, Wolfgang et al. An analytical relation between relaxation time spectrum and molecular weight distribution. Journal of Rheology, v. 43, n. 6, p. 1663-1672, 1999.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref17.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref18 = ttk.Label(second_frame, 
                           text='THIMM, Wolfgang et al. Determination of molecular weight distributions from Rheological data: an application to polystyrene, polymethylmethacrylate and isotactic polypropylene. Applied Rheology, v. 9, n. 4, p. 150-157, 1999.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref18.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref19 = ttk.Label(second_frame, 
                           text='THIMM, Wolfgang et al. On the Rouse spectrum and the determination of the molecular weight distribution from rheological data. Journal of Rheology, v. 44, n. 2, p. 429-438, 2000.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref19.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref20 = ttk.Label(second_frame, 
                           text='TUMINELLO, William H. Molecular weight and molecular weight distribution from dynamic measurements of polymer melts. Polymer Engineering & Science, v. 26, n. 19, p. 1339-1347, 1986.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref20.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref21 = ttk.Label(second_frame, 
                           text='VAN RUYMBEKE, Evelyne; KEUNINGS, Roland; BAILLY, Christian. Determination of the molecular weight distribution of entangled linear polymers from linear viscoelasticity data. Journal of non-newtonian fluid mechanics, v. 105, n. 2-3, p. 153-175, 2002.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref21.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref22 = ttk.Label(second_frame, 
                           text='WASSERMAN, S. H.; GRAESSLEY, W. W. Effects of polydispersity on linear viscoelasticity in entangled polymer melts. Journal of Rheology, v. 36, n. 4, p. 543-572, 1992.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref22.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref23 = ttk.Label(second_frame, 
                           text='WASSERMAN, Scott H. Calculating the molecular weight distribution from linear viscoelastic response of polymer melts. Journal of Rheology, v. 39, n. 3, p. 601-625, 1995.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref23.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref24 = ttk.Label(second_frame, 
                           text='WEESE, Jürgen. A regularization method for nonlinear ill-posed problems. Computer Physics Communications, v. 77, n. 3, p. 429-440, 1993.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref24.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref25 = ttk.Label(second_frame, 
                           text='WU, Souheng. Polymer molecular‐weight distribution from dynamic melt viscoelasticity. Polymer Engineering & Science, v. 25, n. 2, p. 122-128, 1985.',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref25.pack(padx=20, pady=5,anchor=tk.NW)
    
    about_ref26 = ttk.Label(second_frame, 
                           text='WU, Souheng. Characterization of polymer molecular weight distribution by transient viscoelasticity: polytetrafluoroethylenes. Polymer Engineering & Science, v. 28, n. 8, p. 538-543, 1988.\n',
                               font=("Arial", 11), width=150,wraplength=650)
    about_ref26.pack(padx=20, pady=5,anchor=tk.NW)
    
    def OnMouseWheel(self,event):
        self.scrollbar1.yview("scroll",event.delta,"units")
        return "break" 
    
    
    root.bind("<MouseWheel>", OnMouseWheel)
    
    
    
    
    root.mainloop()
