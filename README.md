### Data Analysis Experience
# Gravitational Waves

## Introduction
This set of Python Notebooks introduces you the material
used for the experience on data analysis of gravitational waves.

## Structure of the package
The package contains the following directories

* *data*: Directory to store the data used for this experience;
* *tutorials*: Contains Python Juyter notebooks used as a tutorial for the analysis;
* *code*: Contains a sample Python Juyter notebook to be used for theanalysis;
* *report*: Directory to put the report, such as the Latex code and compiled report in PDF format;

## Getting the package
In order to get the packege with the tutorial and the exercise, you should run a git clone.
For instance 
```
git clone https://github.com/mmphyslab-pi/exp-gw.git
```

## Working with the package
Of course, since this is a git repo, you can do what do you want (commit, push, create branches, etc).
However, if you want to modify and play with the tutorial code and not change the original, you can create a branch:
```
git branch gw-mybranch
git checkout gw-mybranch
git branch 
```
**This package contains the individual work related to the experience, so if you want you can work directly on the master branch.*


## At the end?
When you have completed your taks, please edit the evaluation.md file filling with your name, surname and date.
And of course, please remember to commit and push adding a meaningful comment.
```
git commit -a -m "A comment" 
git push origin master gw-mybranch
```