data = readtable('IKAP_input_T84.csv');
data = table2cell(data);

PSP = readtable('phosphositeplus.csv');
PSP = table2cell(PSP);

[data]=SearchMotifs(data,PSP);

[data_red,kin]=MakeKin(data,4);

[a]=CreateTT(data_red,4,kin);

[AP,AP_tab,K,cost,mincost]=FitActivities(a,data_red,kin,1,50);

writematrix(K,'kinase_act_T84.csv')
writecell(kin,'kinase_T84.csv')