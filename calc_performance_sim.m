function out = calc_performance_sim(x,IMFavg,IMFmed,fs,f,nans,wmM,ftrue)
% calc_performance_sim calculates the average deviation between the true
% signal and the extracted IMFs.
% 
% Inputs:
% x:        True signal (channels x samples)
% IMFavg:   mean IMF across noise realizations
% IMFmed:   median IMF across noise realizations
% fs:       sampling frequency (Hz)
% f:        Simulated center frequency of true signal
% nans:     Number of nans in noise realization
% wmM:      Weigh frequency across time by amplitude (default 1)
% ftrue:    Compare frequency difference across time (1) or using the
% center frequency (~1)
%
% Outputs:
% out:      Struct containing phase difference, frequency difference and
% temporal correlation of the true and extracted signals.
%
% Created by Sofie Therese Hansen 2016, edited jan. 2019.
% Ref: "Unmixing oscillatory brain activity by EEG source localization and
% empirical mode decomposition", by ST Hansen et al.

if nargin<7
    wmM=1;
end
if nargin<8
    ftrue=0;
end
samps=size(x,2);
st=75;samps=samps-2*st;
[~,FREQtrue,PHItrue] = INST_FREQ_local(x);
PHItrue=PHItrue(:,st+1:end-st);FREQtrue=FREQtrue(:,st+1:end-st);

[instAmpAvg,instFreqAvg,PHIAvg] = INST_FREQ_local(IMFavg);% calculating the instantaneous amplitude and the instantaneous frequency
instAmpAvg=instAmpAvg(:,st+1:end-st);instFreqAvg=instFreqAvg(:,st+1:end-st);PHIAvg=PHIAvg(:,st+1:end-st);
if wmM==1
    if ftrue==1
          freqsAvg = abs(nanmean(instFreqAvg*fs-FREQtrue*fs,2)); 
    else  
          freqsAvg = abs(nanmean(instFreqAvg*fs,2)-f);
    end
else
    freqsAvg = abs(nansum(instFreqAvg.*(instAmpAvg./repmat(nansum(instAmpAvg,2),[1,samps]))*fs,2)-f);
end

[instAmpMed,instFreqMed,PHIMed] = INST_FREQ_local(IMFmed);% calculating the instantaneous amplitude and the instantaneous frequency
instAmpMed=instAmpMed(:,st+1:end-st);instFreqMed=instFreqMed(:,st+1:end-st);PHIMed=PHIMed(:,st+1:end-st);
if wmM==1
    if ftrue==1
        freqsMed = abs(nanmean(instFreqMed*fs-FREQtrue*fs,2));        
    else
        freqsMed = abs(nanmean(instFreqMed*fs,2)-f);
    end
else
    freqsMed = abs(nansum(instFreqMed.*(instAmpMed./repmat(nansum(instAmpMed,2),[1,samps]))*fs,2)-f);
end

phaseAvg=NaN(1,3);phaseMed=NaN(1,3);corrAvg=NaN(1,3);corrMed=NaN(1,3);
for fr=1:3
phaseAvg(fr) = nanmean(abs(circ_dist(squeeze(PHIAvg(fr,:)),PHItrue(fr,:))));
phaseMed(fr) = nanmean(abs(circ_dist(squeeze(PHIMed(fr,:)),PHItrue(fr,:))));
co=corrcoef(IMFavg(fr,st+1:end-st),x(fr,st+1:end-st));
corrAvg(fr)=abs(co(1,2));
co=corrcoef(IMFmed(fr,st+1:end-st),x(fr,st+1:end-st));
corrMed(fr)=abs(co(1,2));
end

freqsMed(isnan(phaseMed))=NaN;
freqsAvg(isnan(phaseAvg))=NaN;

%%% based on knowing the frequency band, IMFs averaged/median across reps
out.phaseMedIMFs=phaseMed;
out.phaseAvgIMFs=phaseAvg;
out.freqsMedIMFs=freqsMed;
out.freqsAvgIMFs=freqsAvg;
out.corrAvgIMFs=corrAvg;
out.corrMedIMFs=corrMed;
out.nans = nans;
