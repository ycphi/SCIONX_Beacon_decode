function [summary,bits]=scionx_compare_bit_windows(t,y,bounds,lag,rate,polarity)
% NRZ hard decisions of corresponding waveform windows, before destuff/CRC.
% t,bounds,lag are milliseconds. Recover phase/threshold from EARLIER data
% only; never optimize agreement or edit bits. No time stretching.
t=t(:); y=y(:); period=1000/rate;
assert(numel(t)==numel(y) && all(diff(t)>0) && all(isfinite(y)));
sourceBounds=bounds-lag;
assert(sourceBounds(1)>=t(1) && bounds(2)<=t(end));
source=y(t>=sourceBounds(1) & t<=sourceBounds(2)); sorted=sort(source);
tones=sorted(max(1,round([.1 .9]*numel(sorted))));
threshold=mean(tones); deviation=diff(tones)/2;
assert(deviation>300,'Not enough modulation for bit decisions.');
integral=cumtrapz(t,y); phases=(0:63)'/64*period;
quality=nan(size(phases));
for j=1:numel(phases)
    edges=makeEdges(sourceBounds,period,phases(j));
    levels=diff(interp1(t,integral,edges,'linear'))/period;
    quality(j)=mean((levels-threshold).^2);
end
[~,pick]=max(quality); phase=phases(pick);
edges=makeEdges(sourceBounds,period,phase);
earlier=diff(interp1(t,integral,edges,'linear'))/period;
after=diff(interp1(t,integral,edges+lag,'linear'))/period;
assert(all(isfinite(earlier)) && all(isfinite(after)));
earlierBits=polarity*(earlier-threshold)>0;
afterBits=polarity*(after-threshold)>0;
mismatch=earlierBits~=afterBits;
margin=min(abs(earlier-threshold),abs(after-threshold))/deviation;
% Sensitivity uses the SAME phase offset for both copies; never pick the best.
phaseAgreement=zeros(3,1);
for j=1:3
    e=edges+(j-2)*.1*period;
    keep=e>=sourceBounds(1) & e<=sourceBounds(2);
    e=e(keep);
    a=diff(interp1(t,integral,e,'linear'))/period;
    b=diff(interp1(t,integral,e+lag,'linear'))/period;
    phaseAgreement(j)=100*mean((a>threshold)==(b>threshold));
end
centers=(edges(1:end-1)+edges(2:end))/2;
bits=table((1:numel(after))',centers,centers+lag,earlier,after, ...
    double(earlierBits),double(afterBits),mismatch,margin,margin<.25, ...
    'VariableNames',{'BitIndex','SourceCenter_ms','AfterCenter_ms','EarlierMeanHz', ...
    'AfterMeanHz','EarlierBit','AfterBit','Mismatch','MinimumDecisionMargin', ...
    'LowDecisionMargin'});
summary=table(rate,phase,threshold,lag,numel(after),sum(mismatch), ...
    100*mean(~mismatch),sum(margin<.25),sum(mismatch & margin<.25), ...
    min(phaseAgreement),max(phaseAgreement), ...
    'VariableNames',{'SymbolRate','SourcePhase_ms','ThresholdHz','Lag_ms', ...
    'ComparedBits','DifferentBits','Agreement_pct','LowMarginBits', ...
    'LowMarginDifferentBits','PhaseSensitivityMin_pct','PhaseSensitivityMax_pct'});
end

function edges=makeEdges(bounds,period,phase)
first=ceil((bounds(1)-phase)/period); last=floor((bounds(2)-phase)/period);
edges=(first:last)'*period+phase;
assert(numel(edges)>=17,'Need at least 16 complete symbol intervals.');
end
