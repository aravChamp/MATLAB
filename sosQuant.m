clc
clear
close all


rp = 1.5;   
rs = 30;    
fpass = [0.3 0.6];  
Fs = 2;     


[n, Wn] = ellipord(fpass, fpass+[0.1 0.1], rp, rs);
[b, a] = ellip(n, rp, rs, Wn, 'bandpass');

[z, p, k] = tf2zpk(b, a);


disp('Zeros:');
disp(z);
disp('Poles:');
disp(p);
disp('Gain:');
disp(k);