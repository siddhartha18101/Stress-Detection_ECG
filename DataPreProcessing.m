% Author : Alluri L S V Siddhartha Varma

clc;
clear all;
close all;

% Plotting and manuplating a sample of ecg dataset from WESAD dataset
% For any other details refer "README.txt" in wesad dataset

fs = 700; % Sampling Frequency of ECG data

load('test'); % Imported ecg data (10,000 samples) from dataset

% Only 10,000 samples are considered, as a sample of whole data

ecg = ((test/(2^16))-0.5)*(3/1000); % Transfer function for the output of hardware used to get ecg

% For any other details refer "README.txt" in wesad dataset and also
% Datasheet of hardware used 

% Getting the co-ordinates of the time axis
t = (0:(length(ecg)-1))/fs; %(time sample = sample*(1/smapling frequency))
t = t.';

%Plotting the ecg,
subplot(4,1,1);
plot(t,ecg);
title("ECG Signal (Original!)");
xlabel("time");
ylabel("Amplitude");

% Applying the bandpass filter between 1-50Hz inorder to eliminate the noise
% y = bandpass(ecg,[1,50],fs);



% Moving median Filter
k = 150; % width of the filter
M = movmedian(ecg,k);
subplot(4,1,2);
plot(t,ecg-M);
title("ECG Signal (BW removal)");
xlabel("time");
ylabel("Amplitude");

% bw = zeros(size(ecg));
% bw(1) = 
% for i=2:size(ecg)
%     bw(i) = (((-1/2)*M(i-1)) + ((3/2)*M(i)) - ((3/2)*M(i+1)) + ((1/2)*M(i+2)))
%     

% Cutting line for noise elimination
si = zeros(size(ecg));
ni = zeros(size(ecg));
cl = zeros(size(ecg));

necg = zeros(size(ecg));
alpha = 3;

for i=1:length(ecg)-2
    win = [ecg(i);ecg(i+1);ecg(i+2)];
    si(i) = max(win);
    ni(i) = 2*std(win);
    cl(i) = ni(i) + (si(i)-ni(i))*((10-alpha)/9);
    if ecg(i) > cl(i)
        necg(i) = ecg(i)-cl(i);
    else 
        necg(i) = 0;
    end
%      necg(i) = necg(i)/(si(i)-ni(i));
end    
       
% necg = necg/(si-ni);
subplot(4,1,3);
plot(t,necg);
title("ECG Signal (Noise removed)");
xlabel("time");
ylabel("Amplitude");

ecg = necg;




qrsEx = ecg(4500:4800);
[mpdict,~,~,longs] = wmpdictionary(numel(qrsEx),'lstcpt',{{'sym4',3}});
% figure
% plot(qrsEx)
% hold on
% plot(2*circshift(mpdict(:,11),[-2 0]),'r')
% axis tight
% legend('QRS Complex','Sym4 Wavelet')
% title('Comparison of Sym4 Wavelet and QRS Complex')
% 
wt = modwt(ecg,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');

y = abs(y).^2;
[qrspeaks,locs] = findpeaks(y,t);
subplot(4,1,4);
plot(t,y);
hold on;
plot(locs,qrspeaks,'o')
xlabel("time");
ylabel("Amplitude");
title('R Peaks Localized by Wavelet Transform with Automatic Annotations');

% [pxx,f] = plomb(qrspeaks,locs);
% subplot(4,1,4);
% plot(f,pxx);
% xlabel("Frequency");
% ylabel("Power");
% title("Lomb Periodogram");
