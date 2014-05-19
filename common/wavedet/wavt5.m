function wt = wavt5(x,q1,q2,q3,q4,q5)
% Calculates the wavelet transform using quadratic spline wavelet.
% It calculates scales 1 to 4. 
% Juan Pablo Martínez Cortés
% Last update: Rute Almeida  09.07.2007 
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

lx = length(x);
wt = zeros(lx,5);

% Implementazión con comboluzión en o dominio temporal
% Millor filter que comv, ta que sólo bi aiga muestras buenas.
% Feremos una mena d'overlap-save.

% Implementation as a convolution in the temporal domain
% filter rather than conv, so as to ther are only "good" samples

%l1 = length(q1); l2 = length(q2); l3 = length(q3); l4=length(q4);
wt(:,1) = filter(q1,1,x)';   %conv(x,q1)';
wt(:,2) = filter(q2,1,x)';   %conv(x,q2)';
wt(:,3) = filter(q3,1,x)';   %conv(x,q3)';
wt(:,4) = filter(q4,1,x)';   %conv(x,q4)';
wt(:,5) = filter(q5,1,x)';   %conv(x,q5)';
%figure(1)
%hold off
%plot (wt(:,1))
%hold on
%plot (wt(:,2)-10)
%plot (wt(:,3)-20)
%plot (wt(:,4)-30)
%plot (x-40)


% Implementazión en o dominio frecuenzial
% Frequencial domain
%fq1 = fft(q1, 2^(nextpow2(lx+l4-1)));
%fq2 = fft(q2, 2^(nextpow2(lx+l4-1)));
%fq3 = fft(q3, 2^(nextpow2(lx+l4-1)));
%fq4 = fft(q4, 2^(nextpow2(lx+l4-1)));
%fx = conj(fft(x, 2^(nextpow2(lx+l4-1))))';
%wt2(:,1) = real(ifft(fq1.*fx));
%wt2(:,2) = real(ifft(fq2.*fx));
%wt2(:,3) = real(ifft(fq3.*fx));
%wt2(:,4) = real(ifft(fq4.*fx));
%figure(2)
%hold off
%plot (wt2(2:end,1))
%hold on
%plot (wt2(4:end,2)-10)
%plot (wt2(8:end,3)-20)
%plot (wt2(16:end,4)-30)
%plot (x-40)

