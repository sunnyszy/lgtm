% a = cumsum(rand(100,100)*0.5);
% a = mod(a,2*pi);
% a = unwrap(a,2*pi, 2);
fit_X = [1 2 3 4 5];
fit_Y = [1 8 27 64 125];
result = polyfit(fit_X, fit_Y, 3);