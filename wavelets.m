clc
clear
close all

N = 1; % substitute a value between 1 and 45
wname = ['db', int2str(N)];
[h0, h1, f0, f1]= wfilters(wname);

w = linspace(0, pi, 1e4);

H0_mag = abs(freqz(h0, 1, w));
H1_mag = abs(freqz(h1, 1, w));

figure;
plot(w, H0_mag);
title('Magnitude Responses |H_0(\omega)| & |H_1(\omega)|');
xlabel('Frequency (\omega)');
ylabel('|H(\omega)|');
hold on;
plot(w, H1_mag);
grid on;
legend('|H_0(\omega)|', '|H_1(\omega)|');
hold off;

E = [h0; h1];

% coefficeints are real valued
%reverse the order

B = rot90(E, 2)*E; % = I 2x2

P = H0_mag.^2 + H1_mag.^2;

pwrVariation = max(P) - min(P); % close to zero

avgPwr = mean(P);










