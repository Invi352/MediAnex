# MediAnex
The Mediator Analysis Nexus (MediAnex) application was built for the lipid mediator unit of Queen Mary University with the aim of putting multiple lipid mediator analysis tools into one user friendly application.
***
## **Table of contents**
1. General Information
2. Installed Libraries & Installation Process
3. Application use
4. Limitations and future development
5. Project Contribution
6. Acknowledgements
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

  
