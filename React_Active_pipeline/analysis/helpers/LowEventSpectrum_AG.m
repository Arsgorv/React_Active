function LowEventSpectrum_AG(filename, ch, struc)
% LowEventSpectrum_AG
% Event-optimized low spectrogram (0.5-10 Hz) for behavioural alignment.
%
% Saves: <struc>_LowEvent_Spectrum.mat with Spectro={Sp,t,f}

if filename(end) ~= filesep
    filename = [filename filesep];
end

if isnumeric(ch)
    load(strcat('LFPData/LFP',num2str(ch),'.mat'));
else
    load(strcat('LFPData/',ch,'.mat'));
end

[params,movingwin,~] = SpectrumParametersBM('low_event'); 

disp('... Calculating LOW_EVENT spectrogramm.');
[Sp,t,f] = mtspecgramc(Data(LFP), movingwin, params);

Spectro = {Sp,t,f};
save(strcat([filename,struc,'_LowEvent_Spectrum.mat']), 'Spectro', 'ch', '-v7.3');
end
