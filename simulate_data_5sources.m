function [sens,x,f,randDist,randFreq,distjit] = simulate_data_5sources(SOI,G,SNR,fs,dist,vert,roi,eq_power,freqDist,QG)
% simulate_data_5sources generated synthetic data with one source of
% interest and 5 distractor sources.
%
% Input:
% SOI:  Location of source of interest as indexed by no. vertex
% G:    Lead field matrix (channels x vertices)
% SNR:  Signal-to-noise ratio (dB)
% fs:   Sampling frequency
% dist: distances of the distractor sources from the SOI (m)
% vert: Vertices of the G (No. vertices x 3)
% roi:  Region of interest
% eq_power: Equal power of planted sources
% freqDist: Which frequency range to use for distractor sources
% QG:   Laplacian like smoothing function (No. vertices x No. vertices)
%
% Output:
% sens: Simulated EEG data
% x:    Individual signal components of SOI
%
% Created by Sofie Therese Hansen 2016, edited jan. 2019.
% Ref: "Unmixing oscillatory brain activity by EEG source localization and
% empirical mode decomposition", by ST Hansen et al.

ndists=length(dist);
[M,N] = size(G);
if nargin<8
    eq_power=0;
end
if nargin<9
    freqDist=1;
end

if nargin<10
    QG=eye(N,N);
end

if SOI==0
    SOI=random('Discrete Uniform',size(vert,1));
end

% Find placements of distractors in specified distances from the SOI and
% not in the defined roi
if dist(1)<1
    for d=1:ndists
        redo=1;
        while redo==1
            distjit(d)=dist(d)+randi([-1000,1000],1)/1e6;
            [~,i]=sort(abs(pdist2(vert(SOI,:),vert)-distjit(d)),'ascend');
            if ismember(i(1),roi)==0 % repeat while if distractor location is in roi
                redo=0;
            end
        end
        randDist(d)=i(1);
    end
else
    randDist=dist;
    distjit=0;
end
f=zeros(1,3);
f(3)=27.45;f(1)=10.1;f(2)=18; % frequency of signal components of SOI

% Generate time-varying signal with specifed center frequencies
Tmax=1;
indx=(1:Tmax*fs)/fs;samps=length(indx);
x1=genSim_time_vary_idinput(f(1),fs,Tmax);
x2=genSim_time_vary_idinput(f(2),fs,Tmax);
x3=genSim_time_vary_idinput(f(3),fs,Tmax);
x=[x1;x2;x3];

data = zeros(N,samps);
data(SOI,:) = sum(x,1);
data_SOI=QG*data; % spatially smooth

% Generate distractor sources
data_dist=zeros(N,samps,ndists);
for d=1:ndists
    data_tmp = zeros(N,samps);
    if strcmp(freqDist,'wide')
        randFreq(d)=1;
    elseif strcmp(freqDist,'narrow')
        randFreq(d)=randi([4 30],1)+rand;
    else
        randFreq(d)=freqDist;
    end
    xdist=genSim_time_vary_idinput(randFreq(d),fs,Tmax);
    data_tmp(randDist(d),:)=xdist;
    data_tmp=QG*data_tmp;
    data_dist(:,:,d)=data_tmp;
end

% Project to sensors
if eq_power==0
    sens = G * data;
else
    sens_SOI=G*data_SOI;
    pow_sens_SOI=mean(sens_SOI(:).^2);
    sens_SOI=sens_SOI/sqrt(pow_sens_SOI);
    sens_Distractors=zeros(size(sens_SOI));
    for d=1:ndists
        sens_dist=G*data_dist(:,:,d);
        pow_sens_dist=mean(sens_dist(:).^2);
        sens_dist=sens_dist/sqrt(pow_sens_dist);
        sens_Distractors=sens_Distractors+sens_dist;
    end
    sens=sens_SOI+sens_Distractors;
end
% Add noise
stdnoise = std(reshape(sens,M*samps,1))*10^(-SNR/20);
noise = randn(M,samps) * stdnoise;
sens = sens + noise; 
