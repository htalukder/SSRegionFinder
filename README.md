SplinesMethylation
==================

Methylation Splines


The code performs Smoothing Splines ANOVA on methylation data. It uses 1000 permutations to produce null distrbutions of specified areas and then compares it with the calculated areas. 


Function inputs: Data matrix, column for response, column for class, column for position, column for individual ID,
                and number of permutations required.
                
                
Function Output: Matrix with Regions (start, end), calculated area under curve and p-value. 
