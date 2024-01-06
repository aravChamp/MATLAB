clc
clear
close all

rp = 1.5;  %passband deviation
rs = 30;   %stop band deviation

[z_ellip, p_ellip, k_ellip] = ellip(4, rp, rs, [0.3 0.6], 'bandpass');

[b_ellip, a_ellip] = zp2tf(z_ellip, p_ellip, k_ellip);

[H_ellip, ~]= freqz(b_ellip, a_ellip, 1000);

mag_ellip = 20*log10(abs(H_ellip));

phaseAng_ellip = unwrap(angle(H_ellip)) * 180 / pi;

f = linspace(0, 1, 1000);

figure;

%subplot(2, 1, 1);
plot(f, mag_ellip, 'b');
title('Magnitude Response');
xlabel('Frequency');
ylabel('Magnitude');

%{
subplot(2, 1, 2);
plot(f, phaseAng_ellip, 'r');
title('Phase Response');
xlabel('Frequency');
ylabel('Phase');
%}

[H_up, G_up] = zp2sos(z_ellip, p_ellip, k_ellip, 'up', 'inf');
[H_dn, G_dn] = zp2sos(z_ellip, p_ellip, k_ellip, 'down', 'inf');

p_up = zeros(4, 2);
z_up = zeros(4, 2);

p_dn = zeros(4, 2);
z_dn = zeros(4, 2);

b_up = zeros(4, 3);
a_up = zeros(4, 3);

b_dn = zeros(4, 3);
a_dn = zeros(4, 3);

mag_res_up = zeros(1000, 4);
mag_res_dn = zeros(1000, 4);

for k = 1:4 
    
    b_up(k, :) = H_up(k, 1:3);
    a_up(k, :) = H_up(k, 4:6);
    
    b_dn(k, :) = H_dn(k, 1:3);
    a_dn(k, :) = H_dn(k, 4:6);
    
    p_up(k, :) = abs(roots(b_up(k, :)));
    z_up(k, :) = abs(roots(a_up(k, :)));
    
    p_dn(k, :) = abs(roots(b_dn(k, :)));
    z_dn(k, :) = abs(roots(a_dn(k, :)));
    
    %mag_res_up(:, k) = abs(freqz(b_up(k, :), a_up(k, :), 1000));
    %mag_res_dn(:, k) = abs(freqz(b_dn(k, :), a_dn(k, :), 1000)); 
end

% The poles and zeros of the up ordering is descending !

% The poles and zeros of the down ordering is ascending ! 

%denom = a_dn(1, :);
num = 1;
denom = a_up(1, :);

figure;
title('Up Case');
hold on;

for k = 1:3
    
    plot(f, 20*log10(abs(freqz(num, denom, 1000))));
    
    num = conv(num, b_up(k, :));
    
    denom = conv(denom, a_up(k+1, :));  
end

plot(f, 20*log10(abs(freqz(conv(conv(b_up(1, :), b_up(2, :)), b_up(3, :)), conv(conv(conv(a_up(1, :), a_up(2, :)), a_up(3, :)), a_up(4, :)), 1000))));

hold off;

num = 1;
denom = a_dn(1, :);

figure;
title('Down Case');
hold on;

for k = 1:3
    
    plot(f, 20*log10(abs(freqz(num, denom, 1000))));
    
    num = conv(num, b_dn(k, :));
    
    denom = conv(denom, a_dn(k+1, :));  
end

plot(f, 20*log10(abs(freqz(conv(conv(b_dn(1, :), b_dn(2, :)), b_dn(3, :)), conv(conv(conv(a_dn(1, :), a_dn(2, :)), a_dn(3, :)), a_dn(4, :)), 1000))));

hold off;

%{
H2 = abs(freqz(1, a_up(1, :), 1000));
H3 = abs(freqz(b_up(1,:), conv(a_up(1,:), a_up(2,:)), 1000));
H4 = abs(freqz(conv(b_up(1,:),b_up(2,:)), conv(conv(a_up(1,:), a_up(2,:)), a_up(3,:)), 1000));
H5 = abs(freqz(conv(conv(b_up(1,:),b_up(2,:)), b_up(3,:)), conv(conv(conv(a_up(1,:), a_up(2,:)), a_up(3,:)), a_up(4,:)), 1000));
figure;
hold on;
plot(f, 20*log10(H2));
plot(f, 20*log10(H3));
plot(f, 20*log10(H4));
plot(f, 20*log10(H5));
hold off;
%}


