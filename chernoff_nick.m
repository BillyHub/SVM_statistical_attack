%numerical results for the Chernoff bounds of adversary test's probability of error
%For separable/non-separable case
%assume f_x1 and f_x2 are symmetric around 0

% refer to: https://www.cs.ubc.ca/~nickhar/W12/Lecture2Notes.pdf

clear all
theta_0 = 5;
m = 300;
Delta = 1;
p_A = 0.3;
p_B = 0.1;

 %distribution = 'uniform';
 %dist_para1 = -5; %for unif_a  
 %dist_para2 = 10; %for unif_b

%distribution = 'exponential';
%dist_para1 = 10; %for exponential distribution mean value mu_e
%dist_para2 = 5; %amount for left-shift, should be positive for non-separable case shift_e

 distribution = 'gaussian'; %only for non-separable case
 dist_para1 = 5; %for mu_g
 dist_para2 = 1; %for sigma

if (strcmp(distribution, 'uniform'))
    %unif_a = -5;
    %unif_b = 10;
    F_X1 = @(x) unifcdf(x, dist_para1, dist_para2);
end

if (strcmp(distribution, 'exponential'))
    %mu_e = 1; %for exponential distribution mean value
    %shift_e = 5; %amount for left-shift, should be positive for non-separable case
    F_X1 = @(x) expcdf(x + dist_para2, dist_para1);
end

if (strcmp(distribution, 'gaussian')) %only for non-separable case
    %mu_g = 1;
    %sigma = 1;
    F_X1 = @(x) normcdf(x, dist_para1, dist_para2);
end

F_X2 = @(x) 1 - F_X1(-x); %assume f_x1 and f_x2 are symmetric
 
%case_study = 'separable';
case_study = 'non-separable';

%chernoff bound
if (strcmp(case_study, 'separable'))
    %need p_A <= 3* p_B
    p_tilde_A = p_A * F_X1(theta_0);
    p_tilde_B =  p_B * F_X1(theta_0);
    bound = 1/2 * exp(-m/(3*p_tilde_A)*((p_tilde_A - p_tilde_B)/2)^2) + 1/2 * exp(-m/(3*p_tilde_B)*((p_tilde_A - p_tilde_B)/2)^2)
end

if (strcmp(case_study, 'non-separable'))
    p_1_A = p_A * F_X1(theta_0);
    p_minus1_A = (1 - p_A) * (1 - F_X2(theta_0));
    p_1_B = p_B * F_X1(theta_0);
    p_minus1_B = (1 - p_B) * (1 - F_X2(theta_0));
    mu_tilde_A = 1/2 * (p_1_A - p_minus1_A) + 1/2;
    mu_tilde_B = 1/2 * (p_1_B - p_minus1_B) + 1/2;
    
    if (mu_tilde_A <= 3 * mu_tilde_B)
        bound = exp(-m * (mu_tilde_A - mu_tilde_B)^2 / (12 * mu_tilde_A))
    else
        "error"
    end
end