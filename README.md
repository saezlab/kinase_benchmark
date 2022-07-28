# kinase_review

In this project we will compare methods that estimate kinase activity from 
phosphoproteomic data.

These methods include:
- RoKAI
- KEA3
- INKA
- IKAP
- KSEA App
- PTM-SEA
- (KARP)
- (KinasePA)

KARP has no code available and KinasePA cannot be compared to the other
methods as it is intended for multiple condition analysis. 

We will prepare the input (CPTAC, perturbation data) for each methods, run them
and compare the results by correlating the activities and calculating a jaccard
index for the top differentially active kinases.
