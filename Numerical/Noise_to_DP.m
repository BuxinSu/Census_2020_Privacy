clear all

% given noise, use local CLT to output epsilon

n = 9;
delta = 10^(-11);
rho = 1.95/100 * 3.325;
sigma = sqrt( 1/(2 * rho/n) );
fprintf('%.20f\n', sigma^2);

epsilon = DP_to_noise(n, sigma, delta);
epsilon_zcdp = rho + 2 * sqrt(- rho * log(delta));
fprintf('%.20f\n', epsilon_zcdp);

fprintf('%.20f\n', epsilon);



% Necessary functions and values
function result = term1(n, epsilon, sigma)
    Bn = sqrt(n) * sigma;
    upper_lim = 50000;
    lower_lim = ceil(Bn * (epsilon* sigma^2 / Bn - n / (2 * Bn)));
    answer = 0;
    for i = lower_lim:upper_lim-1
        answer = answer + normpdf(i/Bn);
    end
    result = answer/Bn;
end

function result = term2(n, epsilon, sigma)
    Bn = sqrt(n) * sigma;
    upper_lim = 50000;
    lower_lim = ceil(Bn * (epsilon* sigma^2 / Bn + n / (2 * Bn)));
    answer = 0;
    for i = lower_lim:upper_lim-1
        answer = answer + normpdf(i/Bn);
    end
    result = answer/Bn;
end

function result = DGM_profile(n, sigma, epsilon)
    t1 = term1(n, epsilon, sigma);
    t2 = term2(n, epsilon, sigma);
    result = t1 - exp(epsilon) * t2;
end


function result = DP_to_noise(n, sigma, delta)
    f = @(x) DGM_profile(n, sigma, x)/delta - 1;
    options = optimset('TolFun', 1e-100);
    result = fzero(f, 10, options);
end



