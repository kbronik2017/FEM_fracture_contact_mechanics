# Performance Analysis Tool

This software  calculates all the following important performance metrics 
(can be used to evaluate algorithmic performance against a manual segmentation):
- Dice Similarity Index (SI) 
- Sensitivity  
- Precision  
- Specificity 
- Jaccardindex  
- Conformity 
- HausdorffDistance
- Overlapcoefficient 
- Tanimotocoefficient
- Matthewscorrelationcoefficient  
- Accuracy 
- Voxel-level False Discovery Rate (FDR) 
- Voxel-level False Negative Ratio (FNR)
- Cluster-level FDR 
- Cluster-level FNR 
- Detection Error Rate (DER) 
- Outline Error Rate (OER)
- Mean Total Area (MTA)
- Volume of segmentation
- Volume of manual mask

# Running the Program!

First user will need to download and extract all the files at the root of folder where all testing folders previously placed, then user can use one of the following methods to run the program:


(1) Using shell script file "RunMe.sh" which can be run by the Unix shell

```sh
  - chmod 755 RunMe.sh  
  - ./RunMe.sh <lesionmask> <threshold> <manualmask> <foldershortname> <modelname> <distancefunction> <factor/opt> <symmetry> <weight> <p_norm> <saveoutput>
  or
  - bash  RunMe.sh <lesionmask> <threshold> <manualmask> <foldershortname> <modelname> <distancefunction> <factor/opt> <symmetry> <weight> <p_norm> <saveoutput>
```  
where
```sh
 - <lesionmask>: Segmented lesion mask using lesion prediction algorithm
 - <threshold>: Probability threshold that will be applied to  output before calculation of the overlap measures
 - <manualmask>: The manual segmentation 
 - <foldershortname>: The testing folders short name
 - <modelname>: The Model Name that has been used to  train the data 
 - <distancefunction>: Euclidean, minkowski, cityblock, cosine, hamming, seuclidean, sqeuclidean, correlation, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, wminkowski 
 - <factor/opt>: Mean, direct 
 - <symmetry>: Symmetric, nonsymmetric 
 - <weight>: Needed to compute the weighted Minkowski distance
 - <p_norm>: The order of norm
 - <saveoutput>: 0, 1
```
(2) Using shell script file "RunMeWithGui.sh" which can be also run by the Unix shell

```sh
  - chmod 755 RunMeWithGui.sh  
  - ./RunMeWithGui.sh
  or
  - bash  RunMeWithGui.sh
``` 
  Here a graphical user interface (GUI) will allow user easily to run the program
 
 
 ![Screenshot](GUI.jpg)
 
 
 and after Help Button Clicked, user will see the following help window which can help user easily fill in the required fields
 
 
 ![Screenshot](GUI-with_help.jpg)
 
 
 and after user will fill the form with necessary information and click OK button, user will see the following message
 
 <br>
 <img height="300" src="AfterClickOK.jpg" />
 </br>

finally after All processes done!, user can see the result in a file "final_result.txt",  saved at the root of folder where user previously has extracted all the program files.
 
