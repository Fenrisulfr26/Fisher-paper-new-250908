function [mome_m0,mome_m1,mome_m2,mome_m2_cen] = cal_random_moments_new(lambdas,rndnum,binWidth,noise_model)
% generate measured moments under random noise
% input lambdas should be a row vector
% generate random moments value
lambda_rep = repmat(lambdas,[rndnum,1]);

if noise_model == "poisson"
obs = poissrnd(lambda_rep,rndnum,length(lambdas)); 
end

tSeries = binWidth*(0:length(lambdas)-1)+binWidth/2;

mome_m0 = sum(obs,2);
mome_m1 = sum(tSeries.*obs,2)./sum(obs,2);
mome_m2 = sum(tSeries.^2.*obs,2)./sum(obs,2);
mome_m2_cen = mome_m2-mome_m1.^2;
end