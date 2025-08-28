path = 'C:\Users\odunuga\OneDrive - Nexus365\TWPAv2\Experimental Work\0815 - Data from spectrum analyser\';

pumpF = 0.5:0.001:3.999;
pumpP = -4;

data = importdata([path,'ThirdHarmonic6_',num2str(pumpF(1),'%0.6f'),'GHz_',num2str(pumpP,'%0.3f'),'dBm.csv']);

fs = data(1:501*9,1);

% SAdata = zeros(length(fs),length(pumpP));
HarmonicData = zeros(length(pumpF),6);

for i=1:length(pumpF)
    
    data = importdata([path,'ThirdHarmonic6_',num2str(pumpF(i),'%0.6f'),'GHz_',num2str(pumpP,'%0.3f'),'dBm.csv']);
    SAdata = data(:,2);

    for n=1:length(SAdata)/501
        freqs = (1:501) + 501*(n-1);        
        x = transpose(freqs);
        y = SAdata(freqs);
        [p, c] = max(y);
        
        HarmonicData(i,n) = p;
    end
end

%% 

% Measurement system S-parameters for normalization
% This is to account for losses/gains in the measurement setup (cables, amplifiers, etc.)
% and isolate the actual device performance.
% Experimenatally we usually see more loss at higher frequencies, so this correction is important.
% The reference file is measured with a thru connection (no device), so it captures the baseline
input = "C:\Users\odunuga\OneDrive - Nexus365\TWPAv2\Experimental Work\Normalisation files\Input.s2p";
output = "C:\Users\odunuga\OneDrive - Nexus365\TWPAv2\Experimental Work\Normalisation files\Ouput_minus_40.s2p";
cryo = "C:\Users\odunuga\OneDrive - Nexus365\TWPAv2\Experimental Work\Normalisation files\Cryostat.s2p";
ref = "C:\Users\odunuga\OneDrive - Nexus365\TWPAv2\Experimental Work\Normalisation files\reference.s2p";

input_sp = sparameters(input); % read S-parameter data from the file
output_sp = sparameters(output);
cryo_sp = sparameters(cryo);
ref_sp = sparameters(ref);

norm_ref = transpose(z2dB(squeeze(ref_sp.Parameters(2,1,:))));

norm_input = transpose(z2dB(squeeze(input_sp.Parameters(2,1,:)))) - norm_ref;
norm_output = transpose(z2dB(squeeze(output_sp.Parameters(2,1,:)))) + 40 - norm_ref;
norm_cryo = transpose(z2dB(squeeze(cryo_sp.Parameters(2,1,:)))) + 40 - norm_ref;


f = input_sp.Frequencies./1e9;

normalizedHarmonicData = 0.*HarmonicData;

for n=3:2:13

    index = (n-1)/2;

    input_at_wp =  transpose(interp1(f, norm_input, pumpF, "linear"));
    cryro_at_wp =  transpose(interp1(f, norm_cryo, pumpF, "linear"));

    output_at_Nwp =  transpose(interp1(f, norm_output, n*pumpF, "linear"));
    cryro_at_Nwp =  transpose(interp1(f, norm_cryo, n*pumpF, "linear"));

    p_wp = pumpP + input_at_wp + 0.5*cryro_at_wp;
    p_Nwp = HarmonicData(:,index) - output_at_Nwp - 0.5*cryro_at_Nwp;

    normalizedHarmonicData(:,n) =  p_Nwp - p_wp;
    
end


input_at_wp =  transpose(interp1(f, norm_input, pumpF, "linear"));
output_at_wp =  transpose(interp1(f, norm_output, pumpF, "linear"));
cryro_at_wp =  transpose(interp1(f, norm_cryo, pumpF, "linear"));
output_at_3wp =  transpose(interp1(f, norm_output, 3*pumpF, "linear"));
cryro_at_3wp =  transpose(interp1(f, norm_cryo, 3*pumpF, "linear"));
output_at_5wp =  transpose(interp1(f, norm_output, 5*pumpF, "linear"));
cryro_at_5wp =  transpose(interp1(f, norm_cryo, 5*pumpF, "linear"));


p_wp = pumpP + input_at_wp + 0.5*cryro_at_wp;
p_3wp = HarmonicData(:,1) - output_at_3wp - 0.5*cryro_at_3wp;
p_5wp = HarmonicData(:,2) - output_at_5wp - 0.5*cryro_at_5wp;

max_freq = 20; % GHz, spectrum analyzer limit
for n = 3:2:13
    index = (n-1)/2;
    harmonic_freq = n * pumpF;
    out_of_range = harmonic_freq > max_freq;
    normalizedHarmonicData(out_of_range, n) = NaN;
end

%%
fig = figure(1);
hold on

% The raw harmonic frequency curves are quite noisy, so some smooting is helpful for visualization and analysis.
window = 1; % This is the smooting window size.This can be adjusted based on how much smoothing is desired. 1 means no smoothing.
harmonics = [3,5,7,9,11,13];
colors = lines(length(harmonics) + 1);
for idx = 1:length(harmonics)
    n = harmonics(idx);
    curve = normalizedHarmonicData(:,n);
    valid = ~isnan(curve);
    if sum(valid) > window % Only smooth if enough valid points
        smoothed = smooth(curve(valid), window);
        plot(pumpF(valid), smoothed, 'Linewidth',2, 'Color', colors(idx + 1,:))
    end
end


xlabel('Pump Power (dBm)')
ylabel('Harmonic Power (dBm)')
title('Nikita Device - Experimental Measured - Output Power vs Frequency')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
legend([arrayfun(@(n) sprintf('%d Harmonic',n), harmonics, 'UniformOutput', false)])
legend("Location",'best')
grid on
set(gcf,'Position',[1000 300 1500 1000])

%%


function [dB] = z2dB(z)  % function to convert complex impedance to dB
    dB = 10*log10(real(z).^2 + imag(z).^2);
end

saveas(fig, 'C:\Users\odunuga\OneDrive - Nexus365\TWPAv2\Experimental Work\Nikita Device Experimentally measured (Plots)\Nikita Device - Experimental Measured - Output Power vs Frequency (3-13 Harmonics).png')