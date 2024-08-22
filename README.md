# MediAnex
The Mediator Analysis Nexus (MediAnex) application was built for the lipid mediator unit of Queen Mary University with the aim of putting multiple lipid mediator analysis tools into one user friendly application.
***
## **Table of contents**
1. General Information
2. Installed Libraries & Installation Process
3. Application use
4. Limitations and future development
5. Project Contribution
***
## General Information
This Application was built for the lipid mediator Unit of Queen Mary university of London, under the supervision of Professor Jesmond Dali and Dr Esteban Alberto Gomez Cifuentes. The lipid mediator unit had developed tools to assist in the analysis of lipid mediator data. For multivariate statistics the tools developed invovled PCA and PLS-DA analysis and Differential analysis among other. Furthermore, the lipid mediator unit has also utilised a variety of libraries to build machine learning (ML) models which could be trained on Specialised Pro-Resolving Mediator (SPMs) concentrations to differentiate between lipid profiles. SPMs are lipid mediators derived from fatty acids and the lipid mediator is part of a community which has found SPMs important in inflammation resolution. The fatty acids studied by the lipid mediator unit are the arachidonic acid (AA), eicosaptaenoic acid (EPA), docosahexaenoic acid (DHA) and docosapentaenoic acid (n-3 DPA). While the lipid mediator unit has created these tools, they exist as separate base R scripts. Furthermore, there was no dedicated tool for pathway analysis and to generate networks of fatty acids and the enzyme pathways leading to the constituent SPMs. Instead the lipid mediator unit has to rely on other platforms such as Cytoscape to build these networks. As such, this project was created with the aim of automating the network generation process allowing user to upload files and instantly show changes in regulation of SPMs. All the tools of the lipid mediator unit were combined in a single user-friendly application to complete lipid mediator data analysis.
***
## Installed Libraries & Installation process:
The installed libraried can be found in the MediAnex_final.R file but primarily we focus on 3 specific packages which are 'visNetwork', 'tidyverse' and 'tidygraph'. The packages were used to create the dataframes for the networks built in this application and also to create visulisations of them. These 3 packages work in convergence first by loading the data from csv files and performing formatting, filtering and creating dataframes using the tidyverse package. These dataframes contain information on the nodes (fatty acids, enzymes and SPMs) and edges (connections between all the different nodes) of each network and can be turned into network objects using the tidygraph package and later visualised using visNetwork. 

Below are images of networks prior and after the addition of lipid mediator data files which made automatic changes to the templates to show changes in regulation. The grey nodes show the template which has not been changed.
![n-3 prior](https://github.com/user-attachments/assets/554b3ffe-a08f-4839-b134-53cbca8c0a8f)

![n-3 DPA](https://github.com/user-attachments/assets/5199424d-f721-4749-b344-e2b9093aa110)

To install these packages and utilise them through an R script use the below code:

  if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
  
  if (!require('tidygraph')) install.packages('tidygraph'); library('tidygraph')
  
  if (!require('visNetwork')) install.packages('visNetwork'); library('visNetwork')

To simply install them through the console run:

  install.packages('tidyverse')
***
## Application use:
To use the application simply download the repository and keep all the files in the same folders. Then by opening R studio, open the script of the application and press the run application button. If the packages are not installed, they will begin installing if the user permits and this can take a while. Once all the packages have succesfully installed and the application is initiated, then the landing page is the PCA/PLS-DA analysis page. If users are unfamiliar with how to perform lipid mediator data analysis, the user can go to the TestFinal folder of this repository and there are csv files available to open, study their format and what data should look like. The application contains turorial buttons and extra information buttons to help make it easier to understand how it works, however some small knowledge of these analyses is required. The application has 7 sections in total which can be accessed by the navigation var next to the title and these are the PCA/PLS-DA, the Differential analysis sectios, the Network analysis and the Machine Learning. The machine learning has 4 sections, 3 of which are used to train ML models and download their files and the last one is to test the ML model build in actual datasets.

***
## Limitations and Future Development:

1. Error Catching:
Even though the application works as intended, there are some things to consider. At present only a small number of errors have been tested and only a small number of errors are informative to new users. The rest of the application utilises the trycatch function which is built into the shiny framework. This function works by preventing the aplpication from crushing and it outputs the error that would appear on the console onto the main panel. While this is still helpful, if a user who does not have bioinformatics experience tries to use this, then the error does not make much sense.

2. Networks displayed in Network analysis:
While the main aim of this project has been met and the networks are being displayed with all the changes being done automatically, the application still utilises csv files and dataframes to build templates of these networks. Some components cannot be changed permanently through the application. For example if a new pathway is discovered or to make new connection between node, changes will need to be made to the csv files. For major changes like those then the files are important, however the boxes which form the SPM families such as Protectins and Maresins in the above image cannot be made smaller or larger in any way. Instead the margins have to be changed through the csv files which is not user-friendly. At the moment there is no update from the visNetwork package (that I have found) which can assist in this issue, however it may be possible to do it programmatically by implementing input fields which connect to the dataframes and make changes in real time.

3. Machine Learning:
While all the ML methods work as intended, some of them take longer than other to complete analysis. The time taken for ML models to run can take between a few seconds to almost a couple of hours. Of course this will also depend on system resources available as well. The only model which takes a very long time to complete is the Extreme Gradient Boosting model and this can take up 1 hour and 40 minutes. It would be an interesting point for future development to change or update the XGBoost package used to run this ML method.

***
## Project contribution:
This project has contributed to making analysis of lipid mediator data an easier and faster process. Previously the tools existed as R scripts and as such, analysis would be conducted primarily by individuals with background knowledge in bioinformatics. Now that these scripts have been integrated into an application it is possible for more people to access and perform this analysis. Furthermore, by utilising this application which is free and built from the ground up using R, it will be easier to upgrade it, update the packages involved and add more to it to fit the needs of the lipid mediator unit. There is now no longer a need to rely on third party platforms or software which can sometimes cost a lot of money and also can be difficult to manage due to different requirements such as file formats. Finally, the network analysis tool built specifically for this project makes the use of software such as cytoscape redundnant for the purpose of analysis lipid mediator data networks. In cytoscape the user would have to change, move and rename each individual node while in the tool built these changes occur automatically and the user only has to select and slightly adjust the placement of nodes for aesthetic purposes. As such, the process of building these networks has turned form a very time consuming process to a much faster automated one. 
