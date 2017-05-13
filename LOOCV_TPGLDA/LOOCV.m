function [auc] = LOOCV(A,B,gama)
%Leave-one-out cross validation is implied on TPGLDA
lncRNAsimilarity = importdata('lncRNAsimilarity.txt');
diseasesimilarity=importdata('diseasesimilarity.txt');
%lncRNAsimilarity represents lncRNA expression similarity as the Spearman correlation coefficient between the expression profiles of each lncRNA pair.
%diseasesimilarity represents the semantic similarities among all diseases.
A_ori = A;
[score_ori] = TPGLDA(A_ori,B,gama);
index = find(1 == A_ori);
for i = 1:length(index)
    i
    A(index(i)) = 0;
    [lncRNA,disease]=ind2sub(size(A),index(i));
    if sum(A(lncRNA,:))==0
        A(lncRNA,:)=Sim_lnc(A,lncRNAsimilarity,lncRNA);
    end
    if sum(A(:,disease))==0
        A(:,disease)=Sim_dis(A,diseasesimilarity,disease);
    end
    [result]=TPGLDA(A,B,gama);
    score_ori(index(i)) = result(index(i));
    A = A_ori;
end
pre_label_score = score_ori(:);
save score_ori;
label_y = A_ori(:);
auc=roc(pre_label_score(:),label_y,'b')
end
