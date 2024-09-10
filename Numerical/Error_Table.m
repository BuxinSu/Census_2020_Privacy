clear all
digits(200);

% List of n and sigma2 values
n_values = [5, 9, 10, 18, 20, 27, 50 ,100];
sigma2_values = [1, 5, 10, 16];

% Loop over the values of n and sigma2
for n = n_values
    for sigma2 = sigma2_values
        disp(['(n,sigma^2): (' num2str(n) ',' num2str(sigma2) ')'])
        B_n = vpa(sqrt( n*sigma2 ));
        % disp(['Period: ' char(vpa(pi * B_n))])
        k = -1000 : 1 : 1000;
        cutoff = floor(double(B_n * pi));
        
        
        
        %%%%%%%%%%%
        %%Omega_1%%
        %%%%%%%%%%%
        % Initialize an array to store the results
        results = zeros(1, cutoff);
        % Loop over t from 1 to cutoff
        for t = 1:cutoff
            char = vpa((sum(exp(-2 .* pi^2 .* k.^2 .* sigma2) .* exp(2 .* k .* pi .* t * sigma2 / B_n)) / sum(exp(-2 .* pi^2 .* k.^2 .* sigma2)))^n);
            result = vpa(char - 1);
            exponential = vpa(exp(- (t-1)^2/2));
            results(t) = vpa( result * exponential / (pi * B_n) );
            % disp(results(t)) 
        end
        
        Omega_1 = vpa( sum(results) );
        scientificNotationStr_1 = sprintf('%.100e\n', Omega_1);
        % disp(['Omega_1:', scientificNotationStr_1])
        
        
        
        %%%%%%%%%%%
        %%Omega_2%%
        %%%%%%%%%%%
        t = cutoff;
        char = vpa((sum(exp(-2 .* pi^2 .* k.^2 .* sigma2) .* exp(2 .* k .* pi .* t * sigma2 / B_n)) / sum(exp(-2 .* pi^2 .* k.^2 .* sigma2)))^n);
        exponential = vpa(exp(- t^2/2) * (pi * B_n - cutoff) / (pi * B_n));
        char_results = vpa(char * exponential);

        Omega_2 = vpa(char_results + exponential);
        scientificNotationStr_2 = sprintf('%.100e\n', Omega_2);
        % disp(['Omega_2:', scientificNotationStr_2])
        
        
        
        %%%%%%%%%%%
        %%Omega_3%%
        %%%%%%%%%%%
        t = pi * B_n;
        exponential = vpa(exp(- t^2/2)/t);
        
        Omega_3 = vpa(exponential / (pi * B_n));
        scientificNotationStr_3 = sprintf('%.100e\n', Omega_3);
        % disp(['Omega_3:', scientificNotationStr_3])
        
        
        
        %%%%%%%%%%%
        %%%Error%%%
        %%%%%%%%%%%
        Omega = Omega_1 + Omega_2 + Omega_3;
        scientificNotationStr = sprintf('%.100e\n', Omega );
        disp(['Error Bound:', scientificNotationStr])
    end
end
