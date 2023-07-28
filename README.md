# PolyWeight 0.5.0
Open source program in Python for determination of MWD of polymers based on analytical rheology. 
It features a user-friendly graphical user interface (GUI),
which offers two distinct approaches for MWD determination: an analytical relation-based method and a parametric
model-based method. By utilizing dynamic moduli, users can calculate MWD as well as molecular weight averages
such as M<sub>n</sub>, M<sub>w</sub>, and M<sub>z</sub>.


## Setting things up

It is recommended to use Anaconda (or Miniconda) to install the required packages and dependencies. Installation instructions for Miniconda (recommended) can be found [here](https://docs.conda.io/en/latest/miniconda.html). It is advisable creating and using a [`conda` environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

* Get the code and all the necessary folders and files:
```
git clone https://github.com/a-minotto/PolyWeight.git
```

* Open Anaconda Prompt (Windows) or a terminal window (Linux or MacOS) and navigate to the PolyWeight directory:
```
cd /path/to/PolyWeight
```

* Create and activate the conda environment:
```
conda create -n PolyWeight-env python=3.10
```

* Activate the environment:
```
conda activate PolyWeight-env
```

* Install the required packages:
```
pip install -r requirements.txt
```
## Directory structure

* **data**
  - **experimental:** Experimentally obtained data from PS extracted from the literature. Dynamic moduli, MWD and relaxation spectrum obtained with NLREG; 
  - **synthetic**
    - **dataset1:** MWD and dynamic moduli from generalized exponential (GEX) distribution. Relaxation spectrum obtained with NLREG;
    - **dataset2:** MWD and dynamic moduli from dataset1 with gaussian noise added. Relaxation spectrum obtained with NLREG;
    - **dataset3:** MWD and dynamic moduli from log-normal (LN) distribution. Relaxation spectrum obtained with NLREG;
    - **Synthetic MWD generation.ipynb:** Jupyter Notebook for synthetic datasets generation;
* **icons:** Windows 10 style .png icons used in the software. Free images available at https://icons8.com/ were used;  
* **logs:** Example log files;
* **materials**
  - **Experimental:** Set of reptation parameter values ​​for the experimental data; 
  - **Synthetic:** Set of reptation parameter values ​​for the synthetic data;
* **about_frame.py:** Code for generating the `About` window/function;
* **help_frame.py:** Code for generating the `Help` window/function;
* **PolyWeight.py:** Main code. Execute this file to run the program and open the graphical user interface;
* **requirements.txt:** Text file that lists the required dependencies for the project;
* **User Manual.pdf:** Set of general instructions and information concerned to the software functioning.
  
## Running a few examples

Run the program:
```
python3 PolyWeight.py
```

For more details about the software, its functionalities, commands, files, etc., open the _"Help!"_ window of the program by clicking on the `Help` button on the top button menu or check the User Manual.

Three synthetic datasets are available: 

* **dataset1:** MWD and dynamic moduli generated from generalized exponential (GEX) distribution; 
* **dataset2:** Dynamic moduli from dataset1 with gaussian noise added; 
* **dataset3:** MWD and dynamic moduli generated from log-normal (LN) distribution.

These datasets were generated in a [Jupyter Notebook](https://jupyter.org/) (available in the `data ⇨ synthetic` folder, within the PolyWeight directory) based on monodispersed linear polymers, closely resembling the characteristics of experimental data for polystyrene. The relaxation spectra of these datasets were obtained using the NLREG program and are also available.

The "Generalized Analytical Relation" method is here referred to as GAR model, while the "Generalized Exponential Distribution"
method is referred as the GEX model.

### Example #1: Synthetic dataset
#### Synthetic data with GAR model

1. In the GAR model tab, the relaxation spectrum must be used as an input file. Click on the `Open RTS` button in the top button menu and navigate to the appropriate directory to select the file:

    `data ⇨ synthetic ⇨ dataset1 ⇨ H_dataset1.dat`

2. If you want to visualize the relaxation spectrum, right-click inside the plot area (around the white area) and select the `Plot relaxation time spectrum` option from the popup menu;

3. Select the set of values titled _"Synthetic"_ from the _"Reptation model parameters"_ field using the dropdown menu and then click the `RUN` button. All program commands will be blocked until the end of execution, when the resulting MWD will be automatically displayed on the screen;

4. To graphically compare the resulting model distribution with the initially generated MWD, click on the `Open MWD` button and navigate to the appropriate directory to select the file:

    `data ⇨ synthetic ⇨ dataset1 ⇨ MWD_dataset1.txt`

    You can choose to view just this distribution (`Plot MWD data only`) or both MWD's (`Plot MWD + estimated distribution`) via the popup menu in the graph area;

5. To copy the resulting values from the _"Average MW's [g/mol] and ratio"_ field, right-click inside this area (close to the edges) and select the `Copy to clipboard` function in the popup menu that will open.

#### Synthetic data with GEX model

1. In the GEX model tab, the dynamic moduli - the storage modulus, G'(ω), and the loss modulus, G''(ω) - must be provided in an input file. Click on the `Open G*(ω)` button in the top button menu and navigate to the appropriate directory to select the file:

    `data ⇨ synthetic ⇨ dataset1 ⇨ G12_dataset1.dat`

2. If you want to visualize the dynamic moduli, right-click inside the plot area (around the white area) and select the `Plot dynamic moduli` option from the popup menu;

3. Select the set of values titled _"Synthetic"_ from the _"Reptation model parameters"_ field using the dropdown menu and then click the `RUN` button. During execution, the parameter values ​​resulting from the optimization process at each step are printed on the terminal. All program commands will be blocked until the end of execution, when the resulting MWD will be automatically displayed on the screen. The console window will also open on the screen, containing information regarding the adjustment of the moduli. This window can be closed and accessed later using the `CONSOLE` button;

   - **NOTE:** This step may take longer compared to the GAR model because it is an optimization process, in which a multiobjective fit is performed;
   
5. To graphically compare the resulting model distribution with the initially generated MWD, click on the `Open MWD` button and navigate to the appropriate directory to select the file:

    `data ⇨ synthetic ⇨ dataset1 ⇨ MWD_dataset1.txt`

    You can choose to view just this distribution (`Plot MWD data only`) or both MWD's (`Plot MWD + estimated distribution`) via the popup menu in the graph area;

6. It is also possible to visualize graphically the adjustment of the dynamic moduli. To do so, right-click in the graphics area and select the option `Plot dynamic moduli + fitted moduli`;

7. To copy the resulting values from the _"Average MW's [g/mol] and ratio"_ field, right-click inside this area (close to the edges) and select the `Copy to clipboard` function in the popup menu that will open.

### Example #2: Experimental dataset

The experimentally obtained data were extracted from the literature and correspond to a monodisperse PS sample. Both the dynamic moduli and the MWD determined via GPC are available. The relaxation spectrum of this dataset was obtained with the NLREG software.

#### Experimental data with GAR model

Follow the same steps described above in the section [Synthetic data with GAR model](https://github.com/a-minotto/PolyWeight_test/tree/main#synthetic-data-with-gar-model). By clicking on the `Open RTS` button, open the appropriate file:

`data ⇨ experimental ⇨ H_experimental.dat`

Both in this case and in the tab referring to the GEX model, if you want to load the MWD determined via GPC for visualization, open the following file with the `Open MWD` function:

`data ⇨ experimental ⇨ MWD_experimental.txt`

#### Experimental data with GEX model

Follow the same steps described above in the section [Synthetic data with GEX model](https://github.com/a-minotto/PolyWeight_test/tree/main#synthetic-data-with-gex-model). By clicking on the `Open G*(ω)` button, open the appropriate file:

`data ⇨ experimental ⇨ G12_experimental.dat`

In this case, it is recommended to use a smaller frequency window. To do so, select the `Change frequency window` option in the _"Settings"_ area. In the window that will open, enter a value for the upper limit frequency (e.g. 5 rad/s - enter numerical value only!) and click on `OK`. Run the program by clicking on the `RUN` button.   

### Other example datasets

The files referring to dataset2 and dataset3 can be found in the `data ⇨ synthetic` folder and can be used following the same steps described in the examples above.

### Log files

Within the PolyWeight directory, some log files can be found in the `logs` folder. These files contain previous analyzes carried out with the data available here and can be opened with the `Open file` function.

The log file must be opened in the tab corresponding to the method used for the analysis. The file name specifies the dataset and resolution method used. For example:

* **dataset3_GAR.txt:** dataset3 data analysis with the GAR model;
* **experimental_GEX.txt:** Analysis of experimental data with the GEX model.
