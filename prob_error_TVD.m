%numerical results for adversary optimal test's probability of error
%For separable/non-separable case
%assume f_x1 and f_x2 are symmetric around 0

clear all
theta_0 = 5;
m = 50;
Delta = 1;
p_A = 0.7;
p_B = 0.5;

%state vector
theta_s = [];
for s = 1 : 2 * m + 1
    theta_s = [theta_s, theta_0 + (s - m - 1) * Delta];
end

 %distribution = 'uniform';
 %dist_para1 = -5; %for unif_a  
 %dist_para2 = 10; %for unif_b

%distribution = 'exponential';
%dist_para1 = 10; %for exponential distribution mean value mu_e
%dist_para2 = 5; %amount for left-shift, should be positive for non-separable case shift_e

 distribution = 'gaussian'; %only for non-separable case
 dist_para1 = 5; %for mu_g
 dist_para2 = 1; %for sigma
 
%case_study = 'separable';
case_study = 'non-separable';

%transition matrix
P_tran_A = tran_mat(distribution, case_study, m, p_A, theta_s, dist_para1, dist_para2);
P_tran_B = tran_mat(distribution, case_study, m, p_B, theta_s, dist_para1, dist_para2);

%initial row vector
pi_0 = zeros(1, 2 * m + 1);
pi_0(m+1) = 1;

pi_m_A = pi_0 * P_tran_A^m;
pi_m_B = pi_0 * P_tran_B^m;

%total variarion distance
TVD = 1/2 * norm(pi_m_A - pi_m_B);

%optimal test's probability of error
P_ERROR_star = 1/2 - 1/2 * TVD