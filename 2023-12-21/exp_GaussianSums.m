% Sebastian J. Schlecht, Wednesday, 03. January 2024
clear; clc; close all;

winLen = 400;
numSamples = 10000;
level = 0.1;

x = randn(numSamples,1)*level;

y = sum(randn(numSamples,winLen)*level,2) / sqrt(winLen);

figure; hold on;
histogram(x)
histogram(y)