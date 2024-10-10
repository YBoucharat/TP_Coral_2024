# TP coral reef dynamics

## Description 
Here is given a coral reef construction model (new version of REEF; Husson et al., 2018; Pastier et al., 2019), considering past eustatic sea-level, reef growth, marine erosion and clastic sedimentation of the eroded materials, modelling on an initial linear slope.

Husson, L., Pastier, A. M., Pedoja, K., Elliot, M., Paillard, D., Authemayou, C., ... & Cahyarini, S. Y. (2018). Reef carbonate productivity during quaternary sea level oscillations. Geochemistry, Geophysics, Geosystems, 19(4), 1148-1164.

Pastier, A. M., Husson, L., Pedoja, K., Bézos, A., Authemayou, C., Arias‐Ruiz, C., & Cahyarini, S. Y. (2019). Genesis and architecture of sequences of Quaternary coral reef terraces: Insights from numerical models. Geochemistry, Geophysics, Geosystems, 20(8), 4248-4272.

## Installation

### Plug and play solution with binder :
You can test the model easily on binder using the following link : 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/YBoucharat/TP_Coral_2024.git/HEAD) <br>
This is the easiest way as all the environment is already charged on binder, but if you quit the page, nothing will be saved. To save, you have to download the file, and upload it again. <br>

### For UGA students, you can install the model on gricad UGA :
With this method, the installation is a bit longer (but really easy, just copy/paste the commands below), but every changes will be saved on your account.<br>
Follow this link and log in with agalan : https://gricad-jupyter.univ-grenoble-alpes.fr/hub/login <br>
Open a terminal window and go to notebooks folder : `cd notebooks` <br> 
Git clone this repository : `git clone https://github.com/YBoucharat/TP_Coral_2024.git` <br>
Go into the repository : `cd TP_Coral_2024`<br>
Create a python virtual environment : `python -m venv TP_venv` (TP_venv will be the name of the virtual environment). <br>
Activate it : `source TP_venv/bin/activate` <br>
Install the libraries : `pip install -r requirements.txt` <br>
Install ipykernel : `pip install ipykernel` <br>
Add the environment to jupyter as a kernel : `python -m ipykernel install --user --name=TP_venv --display-name "Python (TP_venv)"`

### Git clone on your computer
Git clone this repository to have it on your computer : `git clone https://github.com/YBoucharat/TP_Coral_2024.git` <br>
If you choose to clone the repository and work on local, you will have to run `pip install -r requirements.txt` in your python environment. 


## Utilisation
REEF model is available in the REEF folder, and given with a step-by-step jupyter notebook, 'Run_REEF.ipynb', with some exercices. REEF source code is in the Library folder. <br>
If using binder, just follow the instructions in the notebooks.<br>
If using gricad, open 'Run_REEF.ipynb', click on 'Python 3 (ipykernel)' in the top right corner of the window and select 'Python (TP_venv)'. You should be able to run the first cell without error messages, and then continue.
