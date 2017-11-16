close all
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [2 1.5 2];              % Receiver position [x y z] (m)
s = [2 3.5 2];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reverberation time (s)
n = 4096;                   % Number of samples

h = rir_generator(c, fs, r, s, L, beta, n);

y = hilbert(h);
t=(1:length(h))*(1/fs);
% plot(t,h)
% plot(abs(h));
% hold on
% plot(real(y))
% hold on

% plot(t,abs(y))
% e= sqrt(abs(y).^2+h.^2);
[e,ylower] = envelope(a(4,:),3,'rms');
[pks,locs]=findpeaks(-diff(diff(e)));
[pks1,locs1]=findpeaks(e,'MinPeakWidth',1,'Annotate','peaks');
findpeaks(e,'MinPeakWidth',1,'Annotate','peaks');
[pks2,locs2]=findpeaks(e,'MinPeakProminence',0.01);
findpeaks(e,'MinPeakProminence',0.01,'Annotate','extents' );



% subplot(2,1,1)
% plot(e)
plot(t,e)
hold on
plot(t,ylower)
plot(t,h)

plot(locs2,pks2,'*r');
hold off
subplot(2,1,2)
plot(e)
% plot(t,e)
hold on
plot(locs1,pks1,'*r');
hold off
linkaxes
% legend('original','Imaginary Part')
% figure 
% subplot(2,1,1)
% plot(diff(e))
% subplot(2,1,2)
% plot(diff(diff(e)))
% hold on 
% plot(locs,pks,'*r');
