clc
clear
close all

a = [1, -1.5055, 0.12630, -0.3778];
b = [0.1336, -0.0568, 0.0563, 0.1336];

b_bits = 5;

a_quant = fi(a, 1, 6, 3);
b_quant = fi(b, 1, 6, 3);

a1 = [1, -0.4954];
b1 = [-0.4954, 1];

a2 = [1, -1.0101, 0.7626];
b2 = [0.7626, -1.0101, 1];

a1_quant = fi(a1, 1, 6, 3);
b1_quant = fi(b1, 1, 6, 3);

a2_quant = fi(a2, 1, 6, 3);
b2_quant = fi(b2, 1, 6, 3);

a_data = a_quant.data;
b_data = b_quant.data;

a1_data = a1_quant.data;
b1_data = b1_quant.data;

a2_data = a2_quant.data;
b2_data = b2_quant.data;

w_val = [0, pi];

%{
H_og_w = abs(freqz(b, a, exp(1*w_val)));
H_og_quant_w = abs(freqz(b_data, a_data, -exp(1*w_val)));

H_og_w_neg = abs(freqz(b, a, exp(1*w_val)));
H_og_quant_w_neg = abs(freqz(b_data, a_data, -exp(1*w_val)));

error_at_w_0 = abs(H_og_quant_w(1) - H_og_w(1));
error_at_w_pi = abs(H_og_quant_w(2) - H_og_w(2));
%}

%H1
H1_og_w = abs(freqz(b1, a1, exp(1*w_val)));
H1_og_quant_w = abs(freqz(b1_data, a1_data, -exp(1*w_val)));

H1_og_w_neg = abs(freqz(b1, a1, exp(1*w_val)));
H1_og_quant_w_neg = abs(freqz(b1_data, a1_data, -exp(1*w_val)));

error_at_w_0_H1 = abs(H1_og_quant_w(1) - H1_og_w(1));
error_at_w_pi_H1 = abs(H1_og_quant_w(2) - H1_og_w(2));

%H2
H2_og_w = abs(freqz(b2, a2, exp(1*w_val)));
H2_og_quant_w = abs(freqz(b2_data, a2_data, -exp(1*w_val)));

H2_og_w_neg = abs(freqz(b2, a2, exp(1*w_val)));
H2_og_quant_w_neg = abs(freqz(b2_data, a2_data, -exp(1*w_val)));

error_at_w_0_H2 = abs(H2_og_quant_w(1) - H2_og_w(1));
error_at_w_pi_H2 = abs(H2_og_quant_w(2) - H2_og_w(2));

%Ideal gains are 1 for both filters and errors are low

w = linspace(0, pi, 1e4);

H = freqz(b, a, w);

H1 = freqz(b1_data, a1_data, w);

H2 = freqz(b2_data, a2_data, w);

H_Q0 = freqz(b_data, a_data, w);

H_QA = 0.5*(H1 + H2);

disp('Maximum differences:');
disp('|H - H_{Q0}|');
disp(max(abs(H - H_Q0)));
disp('|H - H_{QA}|');
disp(max(abs(H - H_QA)));

figure;
hold on;
plot(w, 20*log10(abs(H)));
plot(w, 20*log10(abs(H_Q0)));
plot(w, 20*log10(abs(H_QA)));
legend('H', 'H_{Q0}', 'H_{QA}');
ylim([-40 0]);
hold off;

disp('Maximum deviations (dB):');
disp('H_{Q0}');
disp(max(20*log10(abs(H)) - 20*log10(abs(H_Q0))));
disp('H_{QA}');
disp(max(20*log10(abs(H)) - 20*log10(abs(H_QA))));

H_pb_max = max(abs(H(1:3000)));
H_Q0_pb_max = max(abs(H_Q0(1:3000)));
H_QA_pb_max = max(abs(H_QA(1:4000)));

H_pb_min = min(abs(H(1:3000)));
H_Q0_pb_min = min(abs(H_Q0(1:3000)));
H_QA_pb_min = min(abs(H_QA(1:3000)));

H_sb_max = max(abs(H(4000:end)));
H_Q0_sb_max = max(abs(H_Q0(4000:end)));
H_QA_sb_max = max(abs(H_QA(4000:end)));

H_sb_min = min(abs(H(3001:end)));
H_Q0_sb_min = min(abs(H_Q0(3001:end)));
H_QA_sb_min = min(abs(H_QA(3001:end)));

disp('Equiripple Characteristics:');

fprintf('H Passband Max: %d \n', H_pb_max);
fprintf('H_Q0 Passband Max: %d \n', H_Q0_pb_max);
fprintf('H_QA Passband Max: %d \n', H_QA_pb_max);

fprintf('H Passband Min: %d \n', H_pb_min);
fprintf('H_Q0 Passband Min: %d \n', H_Q0_pb_min);
fprintf('H_QA Passband Min: %d \n', H_QA_pb_min);

fprintf('H Stopband Max: %d \n', H_sb_max);
fprintf('H_Q0 Stopband Max: %d \n', H_Q0_sb_max);
fprintf('H_QA Stopband Max: %d \n', H_QA_sb_max);

fprintf('H Stopband Min: %d \n', H_sb_min);
fprintf('H_Q0 Stopband Min: %d \n', H_Q0_sb_min);
fprintf('H_QA Stopband Min: %d \n\n', H_QA_sb_min);

H_maxSBgain = 20*log10(H_sb_max);
H_Q0_maxSBgain = 20*log10(H_Q0_sb_max);
H_QA_maxSBgain = 20*log10(H_QA_sb_max);

disp('Maximum Stopband Gain (dB):');
fprintf('H : %d \n', H_maxSBgain);
fprintf('H_Q0 : %d \n', H_Q0_maxSBgain);
fprintf('H_QA : %d \n', H_QA_maxSBgain);

H_pb_groupDelay = (-diff(angle(H(1:4000))) ./ diff(w(1:4000))) / (2 * pi);
H_Q0_pb_groupDelay = (-diff(angle(H_Q0(1:4000))) ./ diff(w(1:4000))) / (2 * pi);
H_QA_pb_groupDelay = (-diff(angle(H_QA(1:4000))) ./ diff(w(1:4000))) / (2 * pi);


figure;
hold on;
plot(w(1:3999), H_pb_groupDelay, '-o');
plot(w(1:3999), H_Q0_pb_groupDelay, '-*');
plot(w(1:3999), H_QA_pb_groupDelay);
legend('H', 'H_{Q0}', 'H_{QA}');
title('Group Delay');
hold off;

%The parallel all-pass realization has better sensitivity
%(it drops off faster)

figure;
%subplot(2, 2, 1);
zplane(b, a);
title('H');

figure;
%subplot(2, 2, 2);
zplane(b_data, a_data);
title('H_{Q0}');

figure;
title('H_{QA}');
hold on;
zplane(b1_data, a1_data);
zplane(b2_data, a2_data);
xlim([-2 2]);
ylim([-2 2]);
hold off;

%The original zero has moved into the unit circle for the approximated
%filters


