# MediAnex
The Mediator Analysis Nexus (MediAnex) application was built for the lipid mediator unit of Queen Mary University with the aim of putting multiple lipid mediator analysis tools into one user friendly application.
***
## **Table of contents**
1. General Information
2. Installed Libraries
3. Installation Process
4. Application use
5. Limitations and future development
6. Project Contribution
7. Acknowledgements
***
## General Information
This Application was built for the lipid mediator Unit of Queen Mary university of London, under the supervision of Professor Jesmond Dali and Dr Esteban Alberto Gomez Cifuentes. The lipid mediator unit had developed tools to assist in the analysis of lipid mediator data. For multivariate statistics the tools developed invovled PCA and PLS-DA analysis and Differential analysis among other. Furthermore, the lipid mediator unit has also utilised a variety of libraries to build machine learning (ML) models which could be trained on Specialised Pro-Resolving Mediator (SPMs) concentrations to differentiate between lipid profiles. SPMs are lipid mediators derived from fatty acids and the lipid mediator is part of a community which has found SPMs important in inflammation resolution. The fatty acids studied by the lipid mediator unit are the arachidonic acid (AA), eicosaptaenoic acid (EPA), docosahexaenoic acid (DHA) and docosapentaenoic acid (n-3 DPA). While the lipid mediator unit has created these tools, they exist as separate base R scripts. Furthermore, there was no dedicated tool for pathway analysis and to generate networks of fatty acids and the enzyme pathways leading to the constituent SPMs. Instead the lipid mediator unit has to rely on other platforms such as Cytoscape to build these networks. As such, this project was created with the aim of automating the network generation process allowing user to upload files and instantly show changes in regulation of SPMs. All the tools of the lipid mediator unit were combined in a single user-friendly application to complete lipid mediator data analysis.
***
## Installed Libraries
Below is a table of all the libraried used in this project. The included libraries are required to build the shiny application, with the relevant interfaces, as well as to complete analysis, build ML models an create plots and tables of results.

if (!require('shiny')) install.packages('shiny'); library('shiny')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('plotly')) install.packages('plotly'); library('plotly')
if (!require('corrplot')) install.packages('plotly'); library('corrplot')
if (!require("BiocManager")) install.packages("BiocManager"); library('BiocManager')
if (!require('mixOmics')) BiocManager::install('mixOmics'); library('mixOmics')
if (!require('ropls')) BiocManager::install('ropls'); library('ropls')
if (!require('DT')) install.packages('DT'); library('DT')
if (!require('MVN')) install.packages('MVN'); library('MVN')
if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('snowfall')) install.packages('snowfall')
if (!require('neldermead')) install.packages('neldermead')
if (!require('optimbase')) install.packages('optimbase')
if (!require('classyfire')) install.packages("classyfire_0.1-2.tar.gz", repos = NULL, type = "source"); library('classyfire')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('arm')) install.packages('arm'); library('arm')
if (!require('xgboost')) install.packages('xgboost'); library('xgboost')
if (!require('DiagrammeR')) install.packages('DiagrammeR'); library('DiagrammeR')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('Ckmeans.1d.dp')) install.packages('Ckmeans.1d.dp'); library('Ckmeans.1d.dp')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('tibble')) install.packages('tibble'); library('tibble')
if (!require('pROC')) install.packages('pROC'); library('pROC')
if (!require('shinythemes')) install.packages('shinythemes');library('shinythemes')
if (!require('bslib')) install.packages('bslib'); library('bslib')
if (!require('shinyWidgets')) install.packages('shinyWidgets'); library('shinyWidgets')
if (!require('shinydashboard')) install.packages('shinydashboard'); library('shinydashboard')
if (!require('colourpicker')) install.packages('colourpicker'); library('colourpicker')
if (!require('shinyalert')) install.packages('shinyalert'); library('shinyalert')
if (!require('reporter')) install.packages('reporter'); library('reporter')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('tidygraph')) install.packages('tidygraph'); library('tidygraph')
if (!require('visNetwork')) install.packages('visNetwork'); library('visNetwork')
if (!require('rintrojs')) install.packages('rintrojs'); library('rintrojs')
if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('pryr')) install.packages('pryr'); library('pryr')
