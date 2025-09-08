function [dlnM0,dlnM1,dlnM2,dlnM2c] = calculate_dlnf(lambdas,rndnum,binWidth,num_of_pdf_bins)
%CALCULATE_DLNF Summary of this function goes here
%   Detailed explanation goes here
% input of lambda should be a row vector
% input binwidth is in unit of ps

rng('default')

lambdas = round(lambdas);
lambda_rep = repmat(lambdas,[rndnum,1]);

obs = poissrnd(lambda_rep,rndnum,length(lambdas)); % default poisson noise
number_of_timebins = length(lambdas); % number of time bins

tSeries = binWidth*(0:length(lambdas)-1)'+binWidth/2;

m0 = sum(obs,2);
Em0 = mean(m0);
stdm0 = std(m0);
lambdas_p_m0 = round(1.6*stdm0*ones(number_of_timebins,1)); % perturbation for each timebin 

m1 = sum(tSeries'.*obs,2)./sum(obs,2);
Em1 = mean(m1);
stdm1 = std(m1);
dm1_dbin = (tSeries-Em1)./Em0;
lambdas_p_m1 = round(1.6*stdm1./dm1_dbin); % perturbation for each timebin

m2 = sum(tSeries'.^2.*obs,2)./sum(obs,2);
Em2 = mean(m2);
stdm2 = std(m2);
dm2_dbin = (tSeries.^2-Em2)./Em0;
lambdas_p_m2 = round(1.6*stdm2./dm2_dbin); % perturbation for each timebin

m2c = m2-m1.^2;
Em2c = mean(m2c);
stdm2c = std(m2c);
dm2c_dbin = ((tSeries-Em1).^2-Em2c)./Em0;
lambdas_p_m2c = round(1.6*stdm2c./dm2c_dbin); % perturbation for each timebin


% figure;
% plot(lambdas_p_m1);hold on;
% plot(lambdas_p_m2);hold on;
% plot(lambdas_p_m2c);hold on;

edge_indx = round(linspace(1,rndnum,num_of_pdf_bins+1));

m0_sort = sort(m0);
m0_equal_filling_edge = m0_sort(edge_indx);
m0_equal_filling_edge = [-Inf;m0_equal_filling_edge(2:end-1);  +Inf];  

m1_sort = sort(m1);
m1_equal_filling_edge = m1_sort(edge_indx);
m1_equal_filling_edge = [-Inf; m1_equal_filling_edge(2:end-1);  +Inf];  

m2_sort = sort(m2);
m2_equal_filling_edge = m2_sort(edge_indx);
m2_equal_filling_edge = [-Inf;  m2_equal_filling_edge(2:end-1); +Inf];  

m2c_sort = sort(m2c);
m2c_equal_filling_edge = m2c_sort(edge_indx);
m2c_equal_filling_edge = [-Inf ; m2c_equal_filling_edge(2:end-1); +Inf];  

% % % % def the PDF matrix for each moment
m0_pdf = zeros(num_of_pdf_bins,number_of_timebins);
m0_pdf_p = zeros(num_of_pdf_bins,number_of_timebins);

m1_pdf = zeros(num_of_pdf_bins,number_of_timebins);
m1_pdf_p = zeros(num_of_pdf_bins,number_of_timebins);

m2_pdf = zeros(num_of_pdf_bins,number_of_timebins);
m2_pdf_p = zeros(num_of_pdf_bins,number_of_timebins);

m2c_pdf = zeros(num_of_pdf_bins,number_of_timebins);
m2c_pdf_p = zeros(num_of_pdf_bins,number_of_timebins);

% begin to perturbate on each timebin
for i = 1:number_of_timebins
    f = waitbar(i/number_of_timebins);
    if lambdas(i) == 0 && lamb
        m0_pdf(:,i) = histcounts(m0,m0_equal_filling_edge,'Normalization', 'probability');
        m0_pdf_p(:,i) = m0_pdf(:,i);

        m1_pdf(:,i) = histcounts(m1,m1_equal_filling_edge,'Normalization', 'probability');
        m1_pdf_p(:,i) = m1_pdf(:,i);

        m2_pdf(:,i) = histcounts(m2,m2_equal_filling_edge,'Normalization', 'probability');
        m2_pdf_p(:,i) =  m2_pdf(:,i);

        m2c_pdf(:,i) = histcounts(m2c,m2c_equal_filling_edge,'Normalization', 'probability');
        m2c_pdf_p(:,i) = m2c_pdf(:,i);
        continue
    else
        TR_counts_p = lambdas;
        TR_counts_p(i) = TR_counts_p(i)+lambda_p;

        [m0_p,m1_p,m2_p,m2c_p] = cal_random_moments_new(TR_counts_p,rndnum,binWidth,'poisson');

        m0_pdf(:,i) = histcounts(m0,m0_equal_filling_edge,'Normalization', 'probability');
        m0_pdf_p(:,i) = histcounts(m0_p,m0_equal_filling_edge,'Normalization', 'probability');

        m1_pdf(:,i) = histcounts(m1,m1_equal_filling_edge,'Normalization', 'probability');
        m1_pdf_p(:,i) = histcounts(m1_p,m1_equal_filling_edge,'Normalization', 'probability');

        m2_pdf(:,i) = histcounts(m2,m2_equal_filling_edge,'Normalization', 'probability');
        m2_pdf_p(:,i) = histcounts(m2_p,m2_equal_filling_edge,'Normalization', 'probability');

        m2c_pdf(:,i) = histcounts(m2c,m2c_equal_filling_edge,'Normalization', 'probability');
        m2c_pdf_p(:,i) = histcounts(m2c_p,m2c_equal_filling_edge,'Normalization', 'probability');

    end
end
close(f) % close the prograss bar


% figure,plot(dlnM1)
% figure,plot(dlnM2)
% figure,plot(dlnM2c)

% imagesc(m0_m1_pdf_p-m0_m1_pdf)
% plot((m0_m1_pdf_p(25,:)-m0_m1_pdf(25,:)))

% % % % merge PDF rows to get rid of the zero value
addpath('D:\zhiguan_wang\OneDrive - University of Glasgow\project in UK 2024\tomography241015\Fisher_information_mcx_clone\fisher_cal_mcx');
[m0_pdf_cal,m0_pdf_p_cal] = merge_pdf_row(m0_pdf,m0_pdf_p);
[m1_pdf_cal,m1_pdf_p_cal] = merge_pdf_row(m1_pdf,m1_pdf_p);
[m2_pdf_cal,m2_pdf_p_cal] = merge_pdf_row(m2_pdf,m2_pdf_p);
[m2c_pdf_cal,m2c_pdf_p_cal] = merge_pdf_row(m2c_pdf,m2c_pdf_p);
rmpath('D:\zhiguan_wang\OneDrive - University of Glasgow\project in UK 2024\tomography241015\Fisher_information_mcx_clone\fisher_cal_mcx');
% calculate the H matrix for each moments
dlnM0 = (log(m0_pdf_p_cal)-log(m0_pdf_cal))./lambda_p;
dlnM1 = (log(m1_pdf_p_cal)-log(m1_pdf_cal))./lambda_p;
dlnM2 = (log(m2_pdf_p_cal)-log(m2_pdf_cal))./lambda_p;
dlnM2c = (log(m2c_pdf_p_cal)-log(m2c_pdf_cal))./lambda_p;


end

