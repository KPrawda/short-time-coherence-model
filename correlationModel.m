function [pred_coh,pred_toastd] = correlationModel(F,T,volatility)
%% coherence model
% % Model the coherence for a given volatility
% % Input:
% % - F = fruequencies over which to calculate coherence
% % - T = time over which to calculate coherence
% % - volatility in s/sqrt(s)
% % Output:
% % - pred_coh = modeled coherence
% % - pred_toastd - modeled TOA standard deviation
    magicFactor = 20;
    pred_coh = exp( - magicFactor * (F .* sqrt(T)*volatility).^2 ); % eq. (14) from the paper    
    pred_toastd = T*volatility;
end

