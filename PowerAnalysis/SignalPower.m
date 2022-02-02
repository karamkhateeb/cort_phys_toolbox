function Power = SignalPower(Signal, SamplingFreq)
% this function calculates the signal power

% how many samples are nan? (for troubleshooting)
numBad = sum(isnan(Signal));

time = (length(Signal)-numBad)/SamplingFreq;
Power = sum(Signal.^2,'omitnan')/time;
end