function [selected,reason,candidates]=select_scionx_postamble_candidate(candidates,signal,rfEnd,options)
% A stable post-group tone supports an anchor, but does NOT certify frame/CRC.
% No preference for the latest time or a fixed offset from the July reference.
arguments
    candidates table
    signal struct
    rfEnd (1,1) double
    options.AmbiguityMargin (1,1) double = .05
    options.FlagCountTolerance (1,1) double = 1
    options.ToneWindow_ms (1,2) double = [1 21]
    options.MaximumToneStd_Hz (1,1) double = 300
    options.ResolveRepeatedPostambles (1,1) logical = false
end
selected=[]; reason="not_found"; n=height(candidates);
candidates.Plausible=false(n,1);
candidates.PostFlagToneStd_Hz=nan(n,1);
candidates.PostFlagToneMean_Hz=nan(n,1);
candidates.PostFlagToneSupported=false(n,1);
candidates.RepeatPartnerAnchor_s=nan(n,1);
candidates.RepeatLag_ms=nan(n,1);
candidates.RepeatCorrelation=nan(n,1);
candidates.RepeatAgreement_pct=nan(n,1);
candidates.RepeatComparedBits=zeros(n,1);
candidates.RepeatTransitions=zeros(n,1);
candidates.RepeatSupported=false(n,1);
if n==0, return; end
% One extra flag must not suppress a similarly strong opposite-polarity group.
candidates.Plausible=candidates.FlagCount>=candidates.FlagCount(1)-options.FlagCountTolerance & ...
    candidates.FitCorrelation>=candidates.FitCorrelation(1)-options.AmbiguityMargin;
t=signal.Time_s(:); smoothed=movmean(signal.DiscriminatorHz(:),max(1,round(.0001*signal.SampleRate)));
mid=signal.ToneMidpointHz;
if isfield(signal,'DeviationHz'), deviation=signal.DeviationHz;
else
    values=sort(signal.DiscriminatorHz(:));
    deviation=diff(values(max(1,round([.1 .9]*numel(values)))))/2;
end
for j=1:n
    bounds=candidates.EndTime_s(j)+options.ToneWindow_ms/1000;
    if bounds(2)>min(rfEnd,t(end)) || bounds(1)<t(1), continue; end
    hit=t>=bounds(1) & t<=bounds(2);
    if sum(hit)<8, continue; end
    spread=std(smoothed(hit)); level=mean(smoothed(hit));
    candidates.PostFlagToneStd_Hz(j)=spread;
    candidates.PostFlagToneMean_Hz(j)=level;
    % A constant near-zero/no-signal interval is not terminal-tone evidence.
    candidates.PostFlagToneSupported(j)=spread<=options.MaximumToneStd_Hz && ...
        abs(level-mid)>=.5*deviation && abs(level-mid)<=1.8*deviation;
end
plausible=find(candidates.Plausible);
if numel(plausible)==1
    selected=plausible; reason="flag_group";
else
    supported=plausible(candidates.PostFlagToneSupported(plausible));
    if numel(supported)==1
        selected=supported; reason="unique_postflag_tone";
    else
        reason="ambiguous";
    end
end
% Keep the legacy decision unless exactly one earlier group has strong
% independently measured repeated data. Enabled by the folder decoder.
if options.ResolveRepeatedPostambles && numel(plausible)>1 && ...
        all(ismember({'BitPolarity','LineCode'},candidates.Properties.VariableNames)) && ...
        isfield(signal,'RFStart_s')
    candidates=repeatEvidence(candidates,signal,rfEnd);
    supported=find(candidates.RepeatSupported);
    if numel(supported)==1 && ~isequal(selected,supported)
        selected=supported; reason="verified_repeated_postamble";
    end
end
end

function c=repeatEvidence(c,s,rfEnd)
% Detected group separation supplies lag; no lag fitting or bit repair.
for j=find(c.Plausible)'
    if c.LineCode(j)~="NRZ", continue; end
    for h=find(c.Plausible)'
        lag=1000*(c.AnchorTime_s(h)-c.AnchorTime_s(j));
        if lag<60 || lag>160 || c.LineCode(h)~="NRZ" || ...
                c.BitPolarity(j)~=c.BitPolarity(h), continue; end
        bounds=[8.5 49.5]; anchor=c.AnchorTime_s(j);
        if anchor+bounds(1)/1000<c.EndTime_s(j) || ...
                anchor+(bounds(1)-lag)/1000<s.RFStart_s || ...
                anchor+bounds(2)/1000>min(rfEnd,s.Time_s(end)), continue; end
        t=(s.Time_s(:)-anchor)*1000; y=s.DiscriminatorHz(:);
        if bounds(1)-lag<t(1), continue; end
        q=(bounds(1):1000/s.SampleRate:bounds(2))';
        earlier=interp1(t,y,q-lag); after=interp1(t,y,q);
        if std(earlier)<300 || std(after)<300, continue; end
        earlier=earlier-mean(earlier); after=after-mean(after);
        correlation=sum(earlier.*after)/sqrt(sum(earlier.^2)*sum(after.^2));
        source=sort(y(t>=bounds(1)-lag & t<=bounds(2)-lag));
        if numel(source)<16 || diff(source(max(1,round([.1 .9]*numel(source)))))/2<=300
            continue
        end
        [bits,detail]=scionx_compare_bit_windows(t,y,bounds,lag,9600,1);
        transitions=sum(diff(detail.EarlierBit)~=0);
        if isfinite(c.RepeatCorrelation(j)) && correlation<=c.RepeatCorrelation(j), continue; end
        c.RepeatPartnerAnchor_s(j)=c.AnchorTime_s(h); c.RepeatLag_ms(j)=lag;
        c.RepeatCorrelation(j)=correlation; c.RepeatAgreement_pct(j)=bits.Agreement_pct;
        c.RepeatComparedBits(j)=bits.ComparedBits; c.RepeatTransitions(j)=transitions;
        c.RepeatSupported(j)=correlation>=.90 && bits.Agreement_pct>=95 && ...
            bits.ComparedBits>=256 && transitions>=20;
    end
end
end
