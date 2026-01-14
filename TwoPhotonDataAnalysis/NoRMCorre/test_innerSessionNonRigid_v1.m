N = 4;
k1 = 1500;%3000
k2 = 300;%300

temp_options = struct;
temp_options.Nr = Nr;
temp_options.Nc = Nc;

template_nk2 = innerSessionNonRigid_v5(template_in,Y,T,options,N,k1,k2,temp_options);

