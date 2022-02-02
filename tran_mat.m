function [P_tran] = tran_mat(distribution, case_study, m, p, theta_s, dist_para1, dist_para2)
%transition matrix

P_tran = zeros(2 * m + 1, 2 * m + 1);

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

if (strcmp(case_study, 'separable'))
    for i = 1 : 2 * m + 1
        for j = 1 : 2 * m + 1
            if (i == 1)
                if (j == 1)
                    P_tran(i, j) = 1;
                else
                    break;
                end
            elseif (i == 2 * m + 1)
                if (j ~= 2 * m + 1)
                    continue;
                else
                    P_tran(i, j) = 1;
                end
            else
                p_left = p * F_X1(theta_s(i)) * (theta_s(i) >= 0);
                p_right = (1 - p) * (1 - F_X2(theta_s(i))) * (theta_s(i) < 0);
                p_stay = 1 - p_left - p_right;
                if (i == j + 1)
                    P_tran(i, j) = p_left;
                end
                if (i == j - 1)
                    P_tran(i, j) = p_right;
                end
                if (i == j)
                    P_tran(i, j) = p_stay;
                end
            end
        end
    end
end

if (strcmp(case_study, 'non-separable'))
    for i = 1 : 2 * m + 1
        for j = 1 : 2 * m + 1
            if (i == 1)
                if (j == 1)
                    P_tran(i, j) = 1;
                else
                    break;
                end
            elseif (i == 2 * m + 1)
                if (j ~= 2 * m + 1)
                    continue;
                else
                    P_tran(i, j) = 1;
                end
            else
                p_left = p * F_X1(theta_s(i));
                p_right = (1 - p) * (1 - F_X2(theta_s(i)));
                p_stay = 1 - p_left - p_right;
                if (i == j + 1)
                    P_tran(i, j) = p_left;
                end
                if (i == j - 1)
                    P_tran(i, j) = p_right;
                end
                if (i == j)
                    P_tran(i, j) = p_stay;
                end
            end
        end
    end
end

end