clear all

% given epsilon, output noise using z-CDP

n = 9;
delta = 10^(-10);
epsilon = 1.25;
result = DP_to_noise(n, delta, epsilon);
disp(result^2);

disp( DGM_profile(n,result,epsilon) )


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

function result = DP_to_noise(n, delta, epsilon)
    f = @(x) DGM_profile(n, x, epsilon)/delta - 1;
    options = optimset('TolFun', 1e-100);
    result = fzero(f, 10, options);
end


