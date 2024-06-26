function [volatility] = findVolatility(time_cor, meas_cor, mask, snr_cor, fb)
% % Fit the volatility 
% % Input:
% % - time_cor = time in seconds of meas_cor
% % - meas_cor = measured correlation between measurements
% % - mask = high energy region
% % - snr_cor = expected correlation based on SNR
% % - fb = center frequencies
% % Output:
% % - volatility in s/sqrt(s)

[numT, numIR, numBands] = size(meas_cor);
for itIR = 1:numIR
    for itBands = 1:numBands
        m = mask(:,itBands);

        T = time_cor(m);
        cor = meas_cor(m,itIR,itBands);
        snr = ones(size(snr_cor(m,itIR,itBands)));
        F = fb(itBands);

        % l1 loss
        loss_fun = @(volatility) sum(abs(correlationModel(F,T,exp(volatility)).*snr - cor));

        % do the fitting in log(volatility) for better scaling
        volatilityMin = -50;
        volatilityMax = -3;
        options = optimset('TolX',0.01);
        %%
        vol = exp(fminbnd(loss_fun,volatilityMin,volatilityMax,options));
        %%
        volatility(itIR,itBands) = vol;    
       
    end
end
end
