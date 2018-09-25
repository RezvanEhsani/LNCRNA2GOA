
######################################################################################################################################################################
#
#
#   - This code predicts functions of lncRNA (Protein-Coding) genes based on Gene Ontology dataset.
#   - Please first download "ExpressionData" and "EnsemblID2GOID" form supplementary and change blew addresses according to your system location of these files.
#   - Main function is Prediction Function;
#
#                       Input: 1- GeneID, Ensembl ID of gene that you want to predict of its function for example: ENSG00000228630
#                              2- Onto, Ontology type that have two options; MF(Molecular Function) or BP(Biological Process)
#                              3- Method, which has five options; Pearson or Spearman or Fisher or Sobolev or combine
#
#                       Output: the output is a list of ontology terms that ordered with respect to FDR values
#                              1- GOID, Gene Ontology ID
#                              2- Ontology, Ontology type (MF or BP)
#                              3- FDR, False Positive Rate
#                              4- Term, description of GOID
#    - There is an example for HOTAIR lncRNA in last line.
#
#
########################################################################################################################################################################

library(plyr)
library(GOSim)

options(stringsAsFactors = FALSE)
ExpressionData = read.table('C:/Program Files/R/lncRNA2function_data.txt', sep='\t', head=TRUE)
EnsemblID2GOID = read.table('C:/Program Files/R/EG2GO.txt', sep='\t',head=TRUE)

ExpressionData_PC <- ExpressionData[(ExpressionData$GeneType == "protein_coding"), ]

EnsemblID_PC <- ExpressionData_PC$GeneID
CutOff <- 250
l <- ncol(ExpressionData)

SobolevMetric <- function(x, y) {
            x1 <- x**2 / (sum(x**2))
            y1 <- y**2 / (sum(y**2))
	    z1 <- x1 - y1
	    FT <- fft(z1)
	    w <- 2*pi*(1:length(FT))/(length(FT))
	    s <- sum((1+w)*abs(FT)**2)**(1/2)
	    return(s)
	    }

FisherMetric <- function(x, y) {
	    x1 <- x**2 / (sum(x**2))
	    y1 <- y**2 / (sum(y**2))
	    t <- x1 * y1
	    return(acos(sum(sqrt(t))))
            }

Enrichment_func <- function(DF_, onto) {
    
	List_Top_Genes <- DF_[c(1:CutOff), 1]
	if(onto == "MF") EG_ <- EnsemblID2GOID[(EnsemblID2GOID[ ,3] == "molecular_function"), ]
	if(onto == "BP") EG_ <- EnsemblID2GOID[(EnsemblID2GOID[ ,3] == "biological_process"), ]

	EG_ <- EG_[(EG_$ensembl_gene_id %in% EnsemblID_PC), ]
	ListOfGos <- EG_[(EG_$ensembl_gene_id %in% List_Top_Genes),2]	
	ListOfGos <- unique(ListOfGos)
	ListOfGos <- ListOfGos[which(!is.na(ListOfGos))]

	TermID2ExtID <- ddply(EG_, .(go_id), function(x) Freq_Anno=nrow(x))
        colnames(TermID2ExtID)[2] <- "Freq_go"
	qTermID2ExtID <- TermID2ExtID[(TermID2ExtID$go_id %in% ListOfGos), ]
        qExtID2TermID <- ddply(EG_, .(go_id), function(x) Freq_go=nrow(x[(x$ensembl_gene_id %in% List_Top_Genes),]))
        colnames(qExtID2TermID)[2] <- "Freq_go"

	qExtID2TermID <- qExtID2TermID[(qExtID2TermID$go_id %in% ListOfGos),2]

	n1 = qExtID2TermID
	n2 = qTermID2ExtID$Freq_go-qExtID2TermID
	n3 = length(EnsemblID_PC) - CutOff - n2
	n4 = rep(CutOff, nrow(qTermID2ExtID))
	qTermID2ExtID <- cbind(qTermID2ExtID, n1, n2, n3, n4)

	qTermID2ExtID <- qTermID2ExtID[(qTermID2ExtID$Freq_go>=5),]

	args.df<-qTermID2ExtID[,c("n1", "n2", "n3", "n4")]
	pvalues <- apply(args.df, 1, function(n)
		     min(phyper(0:n[1]-1,n[2], n[3], n[4], lower.tail=FALSE)))

	GOID <- qTermID2ExtID$go_id
	Ontology <- rep(onto, nrow(args.df))
	Pvalue <- format(pvalues, scientific=TRUE, digits = 3)
	fdr  <- p.adjust(pvalues, method = "fdr", n = length(pvalues))
        TERM <- Term(qTermID2ExtID[ ,1])
        
	D_EN <- data.frame(GOID=GOID, Ontology=Ontology, Pvalue=Pvalue, FDR=format(fdr, scientific=TRUE, digits = 3),Term=TERM)
	D_EN <- na.omit(D_EN)
	D_EN <- D_EN[(order(as.numeric(D_EN$FDR))), ]
	D_EN <- D_EN[(as.numeric(D_EN$FDR)<.05),]
	return(D_EN)
        }

Prediction_Function<-function( GeneID, Onto, Method ){
	Target_EX <- ExpressionData[( ExpressionData$GeneID == GeneID ),]
        if( nrow(Target_EX) == 0) return("The GeneID is not in our database")
	Tareget_EX <- as.numeric( Target_EX[1,c(4:l)] )

	if( Method == "Pearson" )       SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX)))
        if( Method == "Spearman" )      SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX, method = "spearman")))
	if( Method == "Fisher" )  SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) SobolevMetric(as.numeric(x), Tareget_EX))
	if( Method == "Sobolev" ) SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) FisherMetric(as.numeric(x), Tareget_EX))
	if( Method == "combine" ) {
				  SCORE_Pearson  <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX)))
				  SCORE_Spearman <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX, method = "spearman")))
                                  SCORE_Fisher   <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) FisherMetric(as.numeric(x), Tareget_EX))
				  SCORE_Sobolev  <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) SobolevMetric(as.numeric(x), Tareget_EX))
				}

	if( Method == "combine" ) DFScore <- data.frame( EnsemblID_PC, SCORE_Pearson, SCORE_Spearman, SCORE_Sobolev, SCORE_Fisher )
	else DFScore <- data.frame( EnsemblID_PC, SCORE )

	DFScore <- na.omit(DFScore)

	DFScore <- DFScore[(DFScore$EnsemblID_PC!=GeneID),]
	if( Method == "Pearson" | Method == "Spearman")  EnrichResult<-Enrichment_func( DFScore[ rev(order(DFScore[,2])), ], Onto)
	if( Method == "Sobolev" | Method == "Fisher") EnrichResult<-Enrichment_func( DFScore[order(DFScore[,2]), ], Onto)
	if( Method == "combine" ) {
				EnrichPearson  <- Enrichment_func( DFScore[ rev(order(DFScore$SCORE_Pearson)), ], Onto)
				EnrichSpearman <- Enrichment_func( DFScore[ rev(order(DFScore$SCORE_Spearman)), ], Onto)
				EnrichSobolev  <- Enrichment_func( DFScore[     order(DFScore$SCORE_Sobolev) , ], Onto)
				EnrichFisher   <- Enrichment_func( DFScore[     order(DFScore$SCORE_Fisher) , ], Onto)
				EnrichCombine  <- rbind( EnrichPearson, EnrichSpearman, EnrichSobolev,  EnrichFisher )
				EnrichCombine  <- ddply(EnrichCombine, .(GOID), function(x) x[which.min(x$FDR),])
				EnrichResult   <- EnrichCombine[ order(as.numeric(EnrichCombine$FDR)), ]
				}
	if(nrow(EnrichResult)>0){
			rownames(EnrichResult) <- c(1:nrow(EnrichResult))
			return(EnrichResult)
			}
	else print("could not find anything!")

}						
# Example
# Prediction_Function( GeneID="ENSG00000228630", Onto="BP", Method="combine" )
