clc
clear
close all


%Specs

fsamp = 40e6;
fcuts = [9e6 9.5e6 12e6 12.5e6];

rp = 1.5;           
rs = 40;
passBandDev = (10^(rp/20)-1)/(10^(rp/20)+1);
stopBandDev = 10^(-rs/20);
dev = [passBandDev stopBandDev passBandDev];

mags = [1 0 1];



%Kaiser Window Design
[n_kaiser,Wn,beta,ftype] = kaiserord(fcuts,mags,dev,fsamp);
n_kaiser = n_kaiser + rem(n_kaiser,2);
hh = fir1(n_kaiser,Wn,ftype,kaiser(n_kaiser+1,beta),'noscale');

kaisFiltLen = n_kaiser + 1;
disp('Kaiser filter design length');
disp(kaisFiltLen);

zeros_kaiser = hh;
figure;
zplane(zeros_kaiser, 1);
title('Kaiser Window Design')

figure;
stem(hh);
title('Stem plot of Kaiser Filter Coefficients');
xlim([0 n_kaiser]);

[H_kaiser,f_kaiser] = freqz(hh,1,1024,fsamp);
figure;
plot(f_kaiser, 20*log10(abs(H_kaiser)));
title('Kaiser Window Design')
grid;

[pbVar_kaiser, peakSBGain_kaiser] = calcMaxMinPBGain(hh, fsamp, fcuts);
disp('Passband Variation (dB)');
disp(pbVar_kaiser);
disp('Peak Stopband Gain (dB)');
disp(peakSBGain_kaiser);


%Equiripple FIR Filter Design

[n_pm,fo,ao,w] = firpmord(fcuts,mags,dev,fsamp);
b = firpm(n_pm,fo,ao,w);

pmFiltLen = n_pm + 1;
disp('Equiripple FIR Filter Design length');
disp(pmFiltLen);

zeros_pm = b;
figure;
zplane(zeros_pm, 1);
title('Equiripple FIR Filter Design')

figure;
stem(b);
title('Stem plot of Equiripple Filter Coefficients');
xlim([0 n_pm]);

[H_pm,f_pm] = freqz(b,1,1024,fsamp);
figure;
plot(f_pm, 20*log10(abs(H_pm)));
title('Equiripple FIR Filter Design')
grid;

[pbVar_pm, peakSBGain_pm] = calcMaxMinPBGain(b, fsamp, fcuts);
disp('Passband Variation (dB)');
disp(pbVar_pm);
disp('Peak Stopband Gain (dB)');
disp(peakSBGain_pm);

weightRatio = w(1)/w(2);
devRatio = stopBandDev/passBandDev;

disp('Difference between weight and deviation ratios:');
disp(weightRatio - devRatio); %They are the same

disp('Both filters meet the specs, but if they did not, I would manually raise filter order (n).');




function [passBandVar, peakSBGain] = calcMaxMinPBGain(coeff, fs, fcuts)
    
    fpass1 = fcuts(1);
    fstop1 = fcuts(2);
    fstop2 = fcuts(3);
    fpass2 = fcuts(4);

    frequencies = linspace(0, fpass1, 1000);

    peak_gain = 0;
    min_passband_gain = inf;

    for freq = frequencies
        z = exp(1i * 2 * pi * freq / fs); 
        H = polyval(coeff, z);
        gain = abs(H);

        if gain > peak_gain
            peak_gain = gain;
        end
        if gain < min_passband_gain
            min_passband_gain = gain;
        end
    end
    
    frequencies = linspace(fpass2, 20e6, 1000);

    for freq = frequencies
        z = exp(1i * 2 * pi * freq / fs); 
        H = polyval(coeff, z);
        gain = abs(H);

        if gain > peak_gain
            peak_gain = gain;
        end
        if gain < min_passband_gain
            min_passband_gain = gain;
        end
    end
    
    passBandVar = 20*log10(peak_gain)-20*log10(min_passband_gain);
    
    peak_gain = 0;
    
    frequencies = linspace(fstop1, fstop2, 1000);

    for freq = frequencies
        z = exp(1i * 2 * pi * freq / fs); 
        H = polyval(coeff, z);
        gain = abs(H);

        if gain > peak_gain
            peak_gain = gain;
        end
    end
    
    peakSBGain = 20 * log10(peak_gain);
    
end









