% Sebastian J. Schlecht, Tuesday, 06 December 2022
% https://en.wikipedia.org/wiki/Scale_space_implementation#The_discrete_Gaussian_kernel
% t = variance = sigma.^2
% n = x axis

% discreteGaussian(n,sigma.^2) gives similar result to normpdf(n,0,sigma)
% but has better properties for low sigma

function g = discreteGaussian(n,t)
g = zeros(numel(n),numel(t));
for it = 1:numel(n)
    g(it,:) = besseli(n(it),t,1);
end
end