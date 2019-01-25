function [Sx]=genSim_time_vary_idinput(f,fs,t,trans)
% genSim_time_vary_idinput generated time-varying signal
%
% Input:
% f:    center frequencies
% fs:   sampling frequencies
% t:    Length of signal (s)
% trans: transition width around center frequency
%
% Output:
% Sx:   Signal with time-varying frequency
%
% Created by Sofie Therese Hansen 2016, edited jan. 2019.
% Ref: "Unmixing oscillatory brain activity by EEG source localization and
% empirical mode decomposition", by ST Hansen et al.

if nargin<4
    trans=0.5;
end
%%
T=t*fs;
P=1e5;
if f==1
    band = [4 30]/fs*2;
else
    band=[f-trans f+trans]/fs*2;
end
indices_t=[];
count=0;
while isempty(indices_t) && count<1e3
    count=count+1;
    u=idinput(P,'sine',band);
    % find segment of u with specified center frequency and transition
    % width
    [~,instFreqAvg,~] = INST_FREQ_local(u');
    F=instFreqAvg*fs;
    if f==1
        indices=find(abs(F-17)<20);
    else
        indices=find(abs(F-f)<0.5);
    end
    indices_t=find(indices(T:end)-indices(1:end-T+1)==T-1);
end
tstart=indices(indices_t(randi(length(indices_t),1)));
Sx=u(tstart:tstart+T-1);
pSx=mean(Sx.^2);
% Scale amplitude to follow the 1/f power law
if f==1
    Sx=Sx'./sqrt(pSx)*sqrt(1/17);
else
    Sx=Sx'./sqrt(pSx)*sqrt(1/f);
end
