function inter_prof_vec=Sim_dis(A,diseasesimilarity,disease)
%A is the lncRNA-disease interaction profile
A=A';
disease_curr_sim = diseasesimilarity(disease,:);
%calculate isolated disease interaction profile vector
inter_prof_vec = disease_curr_sim * A;
inter_prof_vec= inter_prof_vec';
end
