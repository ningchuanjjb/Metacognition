function rho = corrHisteq_jjb_v1(I1,I2,alpha_histeq)


I1 = adapthisteq(I1,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',alpha_histeq);
I2 = adapthisteq(I2,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',alpha_histeq);
rho = corr2(I1,I2);






