function [final_Rscore]=TPGLDA(A,B,gama)
%The best predictive results achieved by TPGLDA when gama=0.6
%Adjacency matrix A: Konwn binary associations between diseases and lncRNAs
%matrix A is 115*178:115 rows represent 115 lncRNAs and 178 columns represent 178 diseases
%Adjacency matrix B:Konwn binary associations between diseases and genes
%matrix B is 1415*178: 1415 rows represent 1415 genes and 178 columns represent 178 diseases

%lncRNA_disease_interaction: adjacency matrix for the lncRNA_disease associations
%lncRNA_disease_interaction(i,j) means lncRNA i associated with disease j

% ndA:the number of diseases in A
% nlA:the number of lncRNAs in A
[nlA,ndA] = size(A);

% Initialize matrix Rsocre1_lnc_dis
Rscore1_lnc_dis=zeros(nlA,ndA);

%calculate the corresponding weight matrix W in A
for i=1:nlA
        q=bsxfun(@rdivide,repmat(A(i,:),nlA,1).*A,sum(A));
        W_A(i,1:nlA)=1./sum(A,2).*sum(q,2);
end
%calculate the level of consistency between the contribution of resource moved in both directions 
W=W_A';
W=W./(repmat(sum(W),nlA,1));
W_A = (W_A+W);
%obtain the first level of resource score about lncRNA_disease_associations
Rscore1_lnc_dis= W_A*A;

% ndB:the number of diseases in A
% nlB:the number of lncRNAs in A
% ng:the number of genes
[ng,ndB] = size(B);
[nlB,ndB]= size(A);
% Initialize matrix Rsocre2_lnc_dis
Rsocre2_lnc_dis=zeros(nlB,ndB);
%calculate the corresponding weight matrix W in B
for i=1:nlB
        q=bsxfun(@rdivide,repmat(A(i,:),ng,1).*B,sum(B));
        W_B(i,1:ng)=(1./sum(B,2)).*sum(q,2);
end
%obtain the second level of resource score about lncRNA_disease_association by using disease-related genes as collaborative prediction
Rsocre2_lnc_dis=W_B*B;

% Initialize Rscore matrix as zero
final_Rscore=zeros(nlB,ndB);

%construct lncRNA_disease_gene tripartite graph
%calculate the final resource socre Rscore to infer potenial lncRNA-disease associations
final_Rscore=plus(gama*(Rscore1_lnc_dis),(1-gama)*(Rsocre2_lnc_dis));
end


