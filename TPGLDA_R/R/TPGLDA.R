TPGLDA <- function(known_lncRNA_disease_interaction,known_gene_disease_interaction,gama=0.6){
#TPGLDA constructs lncRNA_disease_gene tripartite graph to infer potential lncRNA-disease associations
#The best predictive results achieved by TPGLDA when gama=0.6
#Adjacency matrix A: Konwn binary associations between 178 diseases and 115 lncRNAs
  A=known_lncRNA_disease_interaction
#known_lncRNA_disease_interaction: adjacency matrix for the lncRNA_disease associations
#known_lncRNA_disease_interaction(i,j) means lncRNA i associated with disease j
  A <- as.matrix(A)
# nl:the number of lncRNAs in A
# nd:the number of diseases in A
  nl = dim(A)[1]
  nd = dim(A)[2]
#Initialize matrix Rscore1_lnc_dis and the corresponding weight matrix W
  W = matrix(0,nrow = nl,ncol = nl)
  Rsocre1_lnc_dis = matrix(0,nrow = nl,ncol = nd)

#calculate the corresponding weight matrix W of lncRNAs
  for(j in 1:nl){
    q1=repmat(A[j,],nl,1)*A/repmat(colSums(A),nl,1)
    W[j, ]=1/rowSums(A)*rowSums(q1)
  }
#calculate the level of consistency between the contribution of resource moved in both directions
  Wji=t(W)
  Wj_i=Wji/repmat(colSums(Wji),nl,1)
#obtain the first level of resource score about lncRNA_disease_associations
  Rsocre1_lnc_dis = (W+Wj_i)%*%A
  #return (Rsocre1_lnc_dis)

#Adjacency matrix B:Konwn binary associations between 178 diseases and 1415 genes
  B=known_gene_disease_interaction
  B<-as.matrix(B)
#known_gene_disease_interaction(i,j) means gene i associated with disease j
  # A=known_lncRNA_disease_interaction
  # A <- as.matrix(A)
  # nl = dim(A)[1]
  ng = dim(B)[1]
#Initialize matrix Rscore2_lnc_dis and the corresponding weight matrix W_d
  W_d = matrix(0,nrow = nl,ncol = ng)
  Rsocre2_lnc_dis = matrix(0,nrow = nl,ncol = ng)

#calculate the corresponding weight matrix W of diseases
  for(j in 1:nl){
    q2=repmat(A[j,],ng,1)*B/repmat(colSums(B),ng,1)
    W_d[j,]=1/rowSums(B)*rowSums(q2)
  }
#obtain the first level of resource score about lncRNA_disease_associations
  Rsocre2_lnc_dis = W_d%*%B
  
#calculate the final predictive revelance socre (final_Rscore) to infer potenial lncRNA-disease associations  
  final_Rscore=gama*Rsocre1_lnc_dis+(1-gama)*Rsocre2_lnc_dis
#####Notably,remove this code comment can write the predicted result in a txt.
  write.table(final_Rscore, file = "final_Rscore.txt", row.names = F, quote = F, sep="\t") 
  return(final_Rscore)
}

repmat <- function(A,M=1,N=1) {
  A = t(matrix(A))
  A = matrix(rep(A,M),M,N*dim(A)[2],byrow = TRUE)
  return (A)
}


