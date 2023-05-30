# PyPhyPred
A python toolbox for predicting phenotypes using features from brain imaging. 

# Installation
on a computer that has python3 and conda set up (cluster or laptop), clone the git directory:

$ git clone https://github.com/rezvanfarahi/PyPhyPred.git

- cd to the PyPhyPred directory on your system.
- the environment.yml file includes dependencies, run the following command to make a conda environment
called pyphypred with all dependencies installed:

$ conda env create -f environment.yml

$ conda activate pyphypred (on some systems (e.g. BMRC clusters) you'd need to do source activate instead of conda activate)

- now you're all set up, start python/ipython and carry on!

# Data
See data folder in this Git repository. There is a file called notes.txt that explains what each file contains.


# Interacting with Git and sharing your updated code with everyone:
Important note: please submit any changes to the branch called 'develop'

1- First way: 
After cloning the repository (see Installation above), you can contribute to it:
- cd to the PyPhyPred directory on your system. and do the following steps
- $ git checkout develop
- make any changes you want to the the repository files, locally on your computer
- $ git add .
- $ git commit -m "some description of changes that you made"
- $ git push origin develop

2- Second way: 
You can use the Web IDE in GitLab to edit code, and save changes within the develop branch. 

3- Third way:
For those who do not have Git set up (or don't use it often), you can send me your codes by email/google doc, or path to where you saved it on BMRC/jalapeno clusters and I'll add it to the repository. 

# Receiving latest changes made by other groups:
The following command will allow you to pull the latests updates:
- cd to the PyPhyPred directory on your system.
- $ git pull origin develop  





