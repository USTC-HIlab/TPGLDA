function inter_prof_vec=Sim_lnc(A,lncRNAsimilarity,lncRNA)
lncRNA_curr_sim = lncRNAsimilarity(lncRNA,:);
% calculate isolated lncRNA interaction profile vector
inter_prof_vec = lncRNA_curr_sim* A;
end
