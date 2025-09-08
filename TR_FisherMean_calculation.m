function [FiMean_m0,FiMean_m1,FiMean_m2,FiMean_m2c,dm0_dbin,dm1_dbin,dm2_dbin,dm2c_dbin] = TR_FisherMean_calculation(lambdas,rndnum,binWidth)
%TR_FISHER_CALCULATION calculate the Fisher Information of mean value of TR moments
%   Detailed explanation goes here


rng('default')

lambdas = round(lambdas); 
lambda_rep = repmat(lambdas,[rndnum,1]); % generate random sample 

obs = poissrnd(lambda_rep,rndnum,length(lambdas)); % default poisson noise
number_of_timebins = length(lambdas); % number of time bins
tSeries = binWidth*(0:length(lambdas)-1)'+binWidth/2; % center value of each timebin in ps

%----------------calculate moments from random samples-----------------
m0 = sum(obs,2);
Em0 = mean(m0);
dm0_dbin = ones(number_of_timebins,1); % d(EM0)/d(bincounts)

m1 = sum(tSeries'.*obs,2)./sum(obs,2);
Em1 = mean(m1);
dm1_dbin = (tSeries-Em1)./Em0; % d(EM1)/d(bincounts)

m2 = sum(tSeries'.^2.*obs,2)./sum(obs,2);
Em2 = mean(m2);
dm2_dbin = (tSeries.^2-Em2)./Em0; % d(EM2)/d(bincounts)

m2c = m2-m1.^2;
Em2c = mean(m2c);
dm2c_dbin = ((tSeries-Em1).^2-Em2c)./Em0;% d(EM2c)/d(bincounts)
%----------------calculate moments from random samples-----------------

FiMean_m0 = 1/var(m0);
FiMean_m1 = 1/var(m1);
FiMean_m2 = 1/var(m2);
FiMean_m2c = 1/var(m2c);


end