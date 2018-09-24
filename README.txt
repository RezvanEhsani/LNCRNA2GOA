LNCRNA2GOA

Description:

Long noncoding RNAs (lncRNAs) are non-protein-coding transcripts that play important roles in a wide range
of biological process and diseases. But in most cases their function is poorly understood. LNCRNA2GOA is 
a software written in R that predicts GO annotations for a given lncRNA, named g, using a hyper-geometric 
test on data from a list of the most related protein-coding genes with respect to g. To identify the related genes the
LNCRNA2GOA can use either one of two well-known statistical metrics (Pearson and Spearman), or two geometrical 
metrics (Soboloev and Fisher), or a combination method including all four metrics. 


Prerequisites:

R version >= 3.3.3
Plyr >= 1.8.4
GOSim >= 1.12.0


Running the software:

Please first download "ExpressionData" and "EnsemblID2GOID".
Use main function "Prediction_Function" where inputs and outputs are as following:

* Input:  
   1.GeneID; Ensemble ID of gene for annotation prediction, for example: ENSG00000228630.
   2.Onto; Ontology type, with two possible options: MF(Molecular Function) or BP(Biological Process).
   3.Method; with five possible options: Pearson, Spearman, Sobolev, Fisher, or combine.

* Output:
   1.GOID; Gene Ontology ID.
   2.Ontology; Ontology type (MF or BP).
   3.FDR; False Positive Rate.
   4.Term; description of GOID.

For example, to predict annotation for the lncRNA “ENGS00000228630” (HOTAIR) based on the combine method,
you should use the following command:

>Prediction_Function(GOID=“ENSG00000228630”, Onto=”BP”, Method=”combine”) 




