clc
clear
close all

load("Rodanimg.mat");
load("LilyImg.mat");

%part A

lilyUp = insertZeros(Lilyx);
RodanUp = insertZeros(Rodanx);

figure;
image(Lilyx);
colormap(gray);

figure;
image(lilyUp);
colormap(gray);

% part B

lily_fft = abs(fftshift(fft2(Lilyx)));
figure;
image(lily_fft);
colormap(gray);

lily_up_fft = abs(fftshift(fft2(lilyUp)));
figure;
image(lily_up_fft);
colormap(gray);

%This is showing phase distortion


function resultMatrix = insertZeros(matrix)
    n = length(matrix);
    
    resultMatrix = zeros(2*n);
    
    resultMatrix(1:2:end, 1:2:end) = matrix;
end
