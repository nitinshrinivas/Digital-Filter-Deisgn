
% hlpf[n] = sin((n*pi)/6) / n*pi
% hd[n]   = sin((n-nc)*(pi/6)) / (pi*(n-nc))
clear;
%% part a 
N = 51;
n = 1:1:51;
nc = (N-1)/2;
wn = zeros(1,length(n));
for k=1:N
    wn(k) = 1;
end

hdn=zeros(1,length(n));

for k=1:length(n)
    if(k~=nc)
    hdn(k) = (sin((k-nc)*(pi/6))) / (pi*(k-nc));
    else
       hdn(k)=1/6;
    end
end

h1 = hdn.*wn;   

%% part b
Y = fft(h1,1001);
figure(1);
subplot(3,1,1);
stem(linspace(0,2*pi,51),h1);
title('impulse response in rect window');
subplot(3,1,2);
plot(linspace(0,2,1001),mag2db(abs(Y)/max(abs(Y))));
title('mag res in db With Rect Window');
ylim([-100,10]);
subplot(3,1,3);
plot(linspace(0,2*pi,1001),angle(Y));
title('phase response in rect window');

% this is a low pass FIR filter which will have a linear phase.
%% part c

window = blackman(N);
for k=51:length(n)
    window(k) = 0;
end
h2 = hdn.*(window)';

Y2 = fft(h2,1001);
figure(2);
subplot(3,1,1);
stem(h2);
title('impulse response in blackman');
subplot(3,1,2);
plot(mag2db(abs(Y2)/max(abs(Y2))));
ylim([-100,10]);
title('mag resp in db With Blackman');
subplot(3,1,3);
plot(angle(Y2));
title('phase resp with blacman');

%% part d 

% sidelobes levels are more in the rectangular window filter than the blackman
% window filter
% The transition band in retangular window is shorte than in the blackman
% window filter.

%% part e 

n1=1:1:201;
xn=cos(pi*n1/16)+0.25*sin(n1*pi/2);

yn1=conv(xn,h1);
figure(3);
plot(xn);
hold on;
plot(yn1);
hold off;
title('rectangular window 1st sin noise');

yn2=conv(xn,h2);
figure(4);
plot(xn);
hold on;
plot(yn2);
hold off;
title('blackman window 1st sin noise');

x1n=cos(n1*pi/16)+0.25*randn(1,201);

yn11=conv(x1n,h1);
figure(5);
plot(x1n);
hold on;
plot(yn11);
hold off;
title('rectangular window 2nd gaussian noise');

yn12=conv(x1n,h2);
figure(6);
plot(x1n);
hold on;
plot(yn12);
hold off;
title('backman window 2nd gaussian noise');

%% part f
for k=1:length(n)
    h1n(k)=(-1)^k*h1(k);
end
Y = fft(h1n,1001);
figure(7);
subplot(3,1,1)
stem(h1n);
title('impulse response of HPF in rect window');
subplot(3,1,2);
plot(mag2db(abs(Y)./max(abs(Y))));
ylim([-100,10]);
title('mag res in db of HPF in rect window');
subplot(3,1,3);
plot(angle(Y));
title('phase resp of HPF in rect window');

% this is a high pass FIR filter
