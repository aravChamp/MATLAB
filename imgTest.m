clc
clear
close all

h_lap = (1/6) * [1 4 1; 4 -20 4; 1 4 1];

h_x = (1/4) * [1 0 -1; 2 0 -2; 1 0 -1];

h_y = (1/4) * [-1 -2 -1; 0 0 0; 1 2 1];

% Part A

h_lapsob = conv2(h_x, h_x, 'full') + conv2(h_y, h_y, 'full');

%Part B

[Hf_lap,fx_lap,fy_lap] = freqz2(h_lap);
[Hf_lapsob,fx_lapsob,fy_lapsob] = freqz2(h_lapsob);

%Part C

figure;
colormap(jet);

subplot(2, 2, 1)
contour(fx_lap, fy_lap, abs(Hf_lap), 'Fill', 'on');
title('Contour Plot of H_{Lap}');
xlabel('Frequency');
ylabel('Frequency');

subplot(2, 2, 2)
contour(fx_lapsob, fy_lapsob, abs(Hf_lapsob), 'Fill', 'on');
title('Contour Plot of H_{LapSob}');
xlabel('Frequency');
ylabel('Frequency');

subplot(2, 2, 3)
surf(fx_lap, fy_lap, abs(Hf_lap), 'EdgeColor', 'none');
title('Surface Plot of H_{Lap}');
xlabel('Frequency');
ylabel('Frequency');
zlabel('Magnitude');

subplot(2, 2, 4)
surf(fx_lapsob, fy_lapsob, abs(Hf_lapsob), 'EdgeColor', 'none');
title('Surface Plot of H_{LapSob}');
xlabel('Frequency');
ylabel('Frequency');
zlabel('Magnitude');

%The Laplacian Filter is bandstop
%the sobel laplacian approximation Filter is Bandpass

% You can say they are isotrpoic because they have rotational symmetry

%Part D

load("Rodanimg.mat");
load("LilyImg.mat");

lily = Lilyx;
rodan = Rodanx;

lap_filt_lily = filter2(h_lap, lily);
lap_filt_rodan = filter2(h_lap, rodan);


lapsob_filt_lily = filter2(h_lapsob, lily);
lapsob_filt_rodan = filter2(h_lapsob, rodan);

figure;
colormap(gray);

subplot(1, 3, 1);
image(lily);
title('Original');

subplot(1, 3, 2);
image(lap_filt_lily);
title('H_{Lap}');

subplot(1, 3, 3);
image(lapsob_filt_lily);
title('H_{LapSob}');

figure;
colormap(gray);

subplot(1, 3, 1);
image(rodan);
title('Original');

subplot(1, 3, 2);
image(lap_filt_rodan);
title('H_{Lap}');

subplot(1, 3, 3);
image(lapsob_filt_rodan);
title('H_{LapSob}');




