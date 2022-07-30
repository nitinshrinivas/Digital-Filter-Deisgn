clear;
%% part a

w0=pi/4;
b0=1/(2-2*cos(w0));     % calculated using H(1)=1
%  H(z)=b0*(1-exp(jw0)z^-1)*(1-exp(-jw0)z^-1);
a1=1;
b1=[b0,-2*b0*cos(w0),b0];
n=2001;
figure(1);
freqz(b1,a1,n,'whole');

%% part b
r0=0.9;
w0=pi/4;
b0=(1+r0^2+2*r0*cos(w0))/(2-2*cos(w0));
a2=[1,-2*r0*cos(w0),r0^2];
b2=[b0,-2*b0*cos(w0),b0];
n=2001;
figure(2);
freqz(b2,a2,n,'whole');

%% part c
% whole z-plane except z=0 is the ROC of the FIR filter and so it is both
% stable and causal.

% there are two possibilities for ROC of IIR filter
% |z|<r0 then it is neither stable nor causal
% |z|>r0 then it is both stable and causal and also we usually choose this
%        ROC bcz it is both causal and stable.

%% part d
% fvtool(b1,a1);

% fvtool(b2,a2);

% r0=0.99;
% w0=pi/4;
% b0=(1+r0^2+2*r0*cos(w0))/(2-2*cos(w0));
% a2=[1,-2*r0*cos(w0),r0^2];
% b2=[b0,-2*b0*cos(w0),b0];
% fvtool(b2,a2);

% r0=0.5;
% w0=pi/4;
% b0=(1+r0^2+2*r0*cos(w0))/(2-2*cos(w0));
% a2=[1,-2*r0*cos(w0),r0^2];
% b2=[b0,-2*b0*cos(w0),b0];
% fvtool(b2,a2);

% when we increase r0 we get a better and better notch filter with much
% better magnitude response and approx 0 freq res

%% part e 
fs=8192;  %Hz
xn=rand(1,2*fs)-0.5;
f0=1024;
n=0:1/fs:2-1/fs;
noise=sin(2*pi*f0*n);
xn1=xn+noise;

% sound(xn,fs);
% pause(6); sound(xn1,fs);

y1=filter(b1,a1,xn1);
y2=filter(b2,a2,xn1);

% pause(6);sound(y1,fs);
% pause(6);sound(y2,fs);

%% part f
figure(3);
subplot(2,2,1);
plot(xn(1:100));
title('original signal');
subplot(2,2,2);
plot(xn1(1:100));
title('input signal');
subplot(2,2,3);
plot(y1(1:100));
title('1st filter o/p');
subplot(2,2,4);
plot(y2(1:100));
title('2nd filter o/p');

% figure(4);plot(abs(fft(xn1,1001)));
% figure(5);plot(abs(fft(y1,1001)));
% figure(6);plot(abs(fft(y2,1001)));

% it can be easily verified by listning and by looking at its DFT that the
% filters are working as they should by eliminating pi/4 frequency.