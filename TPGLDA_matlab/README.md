# TPGLDA


Novel prediction of associations between lncRNAs and diseases via lncRNA-disease-gene tripartite graph

======================= Instructions to TPGLDA software (version 1.0.0)

Developer: Liang,Ding(Ding520@mail.ustc.edu.cn) from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China

## **Requirement**

4GB memory
MATLAB R2015a or later

## **Related data information need to first load in tripartite graph** 

- /data/known_lncRNA_disease_interaction.txt
- /data/known_gene_disease_interaction.txt

The first text file _known_lncRNA_disease_interaction.txt_ is a table of konwn binary associations between diseases and lncRNAs. 
The second text file _known_gene_disease_interaction.txt_ is a table of konwn binary associations between genes and disease.

## **Run TPGLDA to infer potential associations between lncRNAs and diseases**

We provide relevant experimental validation data of the lncRNA_disease associations and  gene-disease associations. To analyze these data on TPGLDA to further infer potential associations between lncRNAs and diseases, you should input the appropriate code in the matlab Command Window:
	
    A=textread('known_lncRNA_disease_interaction.txt');
	B=textread('known_gene_disease_interaction.txt');
	TPGLDA(A,B,gama)

Then, the predicted results will be automatically saved in the excel table _./final prediction candidate pairs.xls/_.

## Configurations of TPGLDA
### Related configuration files
     ================================================================================================
    | FILE NAME            | DESCRIPTION                                                            |
    =================================================================================================
    |known_lncRNA_disease_ |The known lncRNA-disease associations is represented by adjacency matrix|
	|interaction.txt       | A in our method,which shows binary associations between diseases and   |
	|                      |lncRNAs.1 represents disease j is associated with lncRNA i,otherwise 0. |
    -------------------------------------------------------------------------------------------------
    |known_gene_disease_   |The gene-disease associations are represented by adjacency matrix B in  |
    |interaction.txt       |in our method,which shows binary associations between diseases and genes|
    |                      |.And 1 represents disease j is associated with gene i,otherwise 0.      |
    -------------------------------------------------------------------------------------------------
    |lncRNA_115.txt        |The corresponding names of 115 lncRNAs in adjacency matrix A.           |
    -------------------------------------------------------------------------------------------------
    |disease_178.txt       |The corresponding names of 178 diseases in adjacency matrix A.          |
    -------------------------------------------------------------------------------------------------
    |lncRNAsimilarity      |The lncRNA expression similarity as the Spearman correlation coefficient|
    |                      |between the expression profiles of each lncRNA pair                     |
    -------------------------------------------------------------------------------------------------
	|diseasesimilarity     | The semantic similarities among all diseases.                          |
    -------------------------------------------------------------------------------------------------

The configurations of TPGLDA can be changed in script file _./TPGLDA.m_, and the descriptions of these parameters are provided below:

    =================================================================================================
    | PARAMETER NAME       | DESCRIPTION                                                            |
    =================================================================================================
    | gama_γ               |parameter γ∈[0, 1] is tunable and used to balance the contribution     |
	|                      |between lncRNAs and genes). The default number is 0.6.                  |
    -------------------------------------------------------------------------------------------------
## **Mainly output variables of TPGLDA**

The descriptions of output variables of TPGLDA are provided below:

    ==================================================================================================
    | VARIABLE NAME        | DESCRIPTION                                                             |
    ==================================================================================================
    | predicted_results    |Predicted_results table shows the predicted results of potential lncRNA- |
    |                      |disease associations in descending order.In order to show the predictions|
    |                      |more clearly, we write the top 10000 potential candidate pairs in "final |
    |                      |prediction candidate pairs.xls".                                            |
    --------------------------------------------------------------------------------------------------
    |final_Rscore          |The final_Rscore vectors is the relevance score of disease-related lncRNA|
    |                      |computed by TPGLDA, and a higher value of a certain final_Rscore presents|
    |                      | a larger potential of the lncRNA to be the disease candidate.           |
    -------------------------------------------------------------------------------------------------
Other notes and output results can be available in the _./TPGLDA.m_.   

 ## TPGLDA for users without MATLAB licenses
For users without MATLAB licenses, we also offer R codes of TPGLDA.The detailed method can be seen in **R** package.

## **Contact**

Please feel free to contact us if you need any help: Ding520@mail.ustc.edu.cn

