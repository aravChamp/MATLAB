clc
clear
close all

%Part (c)

N = 5; % substitute a value between 1 and 45
wname = ['db', int2str(N)];
[h0, h1, f0, f1]= wfilters(wname);

e00 = h0(1:2:end);
e01 = h0(2:2:end);
e10 = h1(1:2:end);
e11 = h1(2:2:end);

E = {e00, e01; e10, e11};

e00n = fliplr(e00);
e01n = fliplr(e10); 
e10n = fliplr(e01);
e11n = fliplr(e11);

E11 = conv(e00n, e00) + conv(e10n, e10);
E12 = conv(e00n, e01) + conv(e10n, e11);
E21 = conv(e01n, e00) + conv(e11n, e10);
E22 = conv(e01n, e01) + conv(e11n, e11);

% The cofficient values (absolute values) lie
% around the center of the vector (z^1 z^0 z^-1)
% this is why it is  proportional

E11(abs(E11) < 0.47) = 0;
E12(abs(E12) < 0.47) = 0;
E21(abs(E21) < 0.47) = 0;
E22(abs(E22) < 0.47) = 0;

EnotE = 2 * [norm(E12), norm(E11); norm(E22), norm(E21)];

%reverse the orders of entries

%polynomial multiplication is convolution of the coefficient vectors

%removed small values

% is proportional to the identity matrix

[gd_h0, ~] = grpdelay(h0);
[gd_h1, ~] = grpdelay(h1);

[gd_f0, ~] = grpdelay(f0);
[gd_f1, ~] = grpdelay(f1);

D_anl = max([gd_h0; gd_h1]);
D_synth = max([gd_f0; gd_f1]);
D_tot = D_anl + D_synth;

% D_tot is the end to end delay

w = linspace(0, pi, 1e4);

H0 = freqz(h0, 1, w);
H1 = freqz(h1, 1, w);

figure;
plot(w, abs(H0));
title('Magnitude Responses |H_0(\omega)| & |H_1(\omega)|');
xlabel('Frequency (\omega)');
ylabel('|H(\omega)|');
hold on;
plot(w, abs(H1));
grid on;
legend('|H_0(\omega)|', '|H_1(\omega)|');
hold off;

H0H1 = H0 .* freqz(h1, 1, w*2); 

H0H0H1 = H0 .* freqz(h0, 1, w*2) .* freqz(h1, 1, w*4);

H0H0H0 = H0 .* freqz(h0, 1, w*2) .* freqz(h0, 1, w*4);


figure;
hold on;
plot(w, abs(H1));
plot(w, abs(H0H1));
plot(w, abs(H0H0H1));
plot(w, abs(H0H0H0));
hold off;

% The passband width increases as frequency increases


P = (1/2) * abs(H1).^2;
P = (1/4) * abs(H0H1).^2 + P;
P = (1/8) * abs(H0H0H1).^2 + P;
P = (1/8) * abs(H0H0H0).^2 + P;

meanP = mean(P);
devP = max(P) - min(P);

%Constant!
