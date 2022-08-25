# A Novel Bayesian Spatial-temporal Approach to Quantify SARS-CoV-2 Testing Disparities for Small Area Estimation
This GitHub repositipry provides data and codes for this paper.       

## Abstract
Objectives: We proposed a novel Bayesian spatial-temporal approach to identify and quantify SARS-CoV-2 testing disparities in the context of small area estimation.         
Methods: In step 1, we used a Bayesian inseparable space-time model framework to estimate the testing positivity rate (TPR) at geographical granular areas of the census block groups (CBGs). In Step 2, we adopted a rank-based approach to compare the estimated TPR and the testing rate, to identify areas with testing deficiency and quantify the number of needed tests.       
Results: We demonstrated the utility of the proposed approach using weekly SARS-CoV-2 infection and testing surveillance data from Cameron County, Texas. We identified the CBGs that had experienced substantial testing deficiency, quantify the number of tests that should have been conducted in these areas, and evaluate the short- and long-term testing disparities.        
Conclusions: Our proposed analytical framework offers policy makers and practitioners a tool for understanding SARS-CoV-2 testing disparities in geographically small communities. The analytical approach could also aid local public health departments in COVID-19 response planning, and inform the intervention programs to improve goal setting and strategy implementation related to SARS-CoV-2 testing uptake.       

## Data
- Simulated data on weekly tests and infection counts and population size for 222 GEOIDs: [COVID data](Data/weekly_data.xlsx)         
- 

## Codes
- Process and calculate testing gaps: [R script](01_TestingGap.R)
- Spatial inseparable modeling: [R script](02_SpatialSepModel.R)
- Prediton Positivity Rate: [R script](03_Predition.R)
