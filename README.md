# PSA-ReD: An R script for generating PSA-ReD plots. 
# Built using R 3.5.2
# Described in the publication titled "Increasing the information provided by probabilistic sensitivity analysis: The relative density plot"


Required software: 
- R
- RStudio

Tutorial: 
- Download 'PSA ReD script - v1.0.2.R' and 'PSA ReD project - v1.0.2.Rproj'
- Make sure that they are in the same folder, preferably a dedicated folder (ie, not the generic 'download' folder)
- Open 'PSA-ReD project - v1.0.2.Rproj'
- Open 'PSA-ReD script - v1.0.2.R'

Now, the script is opened in RStudio and ready to use. 
- Information on using the script is provided in the following places:
    - In the peer-reviewed publication, available at: <link to paper>.
    - In the script, we provide step-by step guidance on how to use it. 
- Lines 1 to 406 should be ran without user input. From line 406, we provide guidance on what user input is required. 
- Parts in between "Do not change the part in between / below" and  "Do not change the part in between / above" should not be altered by the user, these should simply be ran. 
- Datafiles should be stored in the same folder as 'PSA-ReD project - v1.0.2.Rproj' and 'PSA-ReD script - v1.0.2.R'.
- The PSA-ReD images will automatically be saved in this same folder. 

We also provide an example datafile with a premade script that is ready to use with this particular datafile. This script is setup so that it generates Figure 2b from the publication. 

Tutorial for the premade example: 
- Download 'PSA-ReD example data script.R' and 'eHTA case study data, used for figure 2.csv'
- Make sure they are in the same folder. 
- Open 'PSA ReD project - v1.0.2.Rproj'
- Open 'PSA-ReD example data script.R'
- Run lines 1 to 550, runtime on our typical consumer grade desktop was approximately 33 seconds.
- In the same folder as where you put the previous files, you will now find "ExampleFigure.PNG"
