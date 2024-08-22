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

Attempt | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Seconds | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

