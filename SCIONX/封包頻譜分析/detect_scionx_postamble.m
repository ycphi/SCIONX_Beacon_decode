function [candidates, signal] = detect_scionx_postamble(row, options)
%DETECT_SCIONX_POSTAMBLE Find exact consecutive flag groups in an RF burst.
% NRZ by default; optional NRZI. No bit repair, destuffing, or CRC claim.
arguments
    row (1,:) table
    options.SymbolRates (1,:) double {mustBePositive,mustBeFinite} = 9590:5:9635
    options.LineCodes (1,:) string {mustBeMember(options.LineCodes,["NRZ","NRZI"])} = "NRZ"
    options.PhaseCount (1,1) double {mustBeInteger,mustBePositive} = 16
    options.MinimumFlags (1,1) double {mustBeInteger,mustBePositive} = 4
    options.MinimumOffset_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 40
    options.MinimumFitCorrelation (1,1) double {mustBeNonnegative,mustBeFinite} = .65
    options.MaximumFragmentGap_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 1.5
    options.ReceiverLowpassHz (1,1) double {mustBePositive,mustBeFinite} = 8000
    options.Padding_ms (1,1) double {mustBePositive,mustBeFinite} = 300
    options.SignalOnly (1,1) logical = false
end
assert(options.MinimumFlags>=2 && options.MinimumFitCorrelation<=1, ...
    'SCIONX:Postamble:Configuration','MinimumFlags must be >=2 and fit correlation <=1.');
fs = row.SampleRate;
assert(options.ReceiverLowpassHz < fs/2,'SCIONX:Postamble:Nyquist','Lowpass exceeds Nyquist.');
if row.InputFormat=="audio"
    info=audioinfo(row.FilePath);
    assert(info.NumChannels==1,'Audio discriminator input must be mono.');
    first=max(0,row.StartSample-round(options.Padding_ms/1000*fs));
    last=min(info.TotalSamples,row.EndSampleExclusive+round(options.Padding_ms/1000*fs));
    samples=audioread(row.FilePath,[first+1 last]);
    time=first/fs+(0:numel(samples)-1)'/fs;
    active=time>=row.Start_s+.005 & time<=row.End_s-.005;
    fm=conv(samples,lowpass(options.ReceiverLowpassHz,fs),'same');
    values=sort(fm(active));
    assert(~isempty(values) && all(isfinite(fm)),'Invalid audio discriminator samples.');
    tones=values(max(1,round([.1 .9]*numel(values))));
    audioScale=2400/max(diff(tones)/2,eps);
    fm=(fm-mean(tones))*audioScale*row.DisplayPolarity;
    fm=conv(fm,[.25;.5;.25],'same'); offset=0;
else
switch row.InputFormat
    case "cs16", bytes=4; primitive='int16=>double'; scale=32768;
    case "cf32", bytes=8; primitive='single=>double'; scale=1;
    otherwise, error('SCIONX:Postamble:Format','Unsupported IQ format.');
end
info=dir(row.FilePath); total=info.bytes/bytes;
assert(mod(info.bytes,bytes)==0,'SCIONX:Postamble:File','Incomplete I/Q pair.');
first=max(0,row.StartSample-round(options.Padding_ms/1000*fs));
last=min(total,row.EndSampleExclusive+round(options.Padding_ms/1000*fs));
fid=fopen(row.FilePath,'rb','ieee-le'); assert(fid>=0,'Cannot open IQ file.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
assert(fseek(fid,bytes*first,'bof')==0,'Cannot seek IQ file.');
raw=fread(fid,[2,last-first],primitive);
assert(all(isfinite(raw(:))),'SCIONX:Postamble:Nonfinite','Nonfinite IQ samples.');
z=complex(raw(1,:).',raw(2,:).')/scale;
if row.DisplayPolarity==-1, z=conj(z); end
% Estimate carrier before a narrower receiver lowpass to tolerate Doppler.
h=lowpass(min(fs*.4,14000),fs);
wide=conv(z,h,'same');
time=first/fs+((0:numel(z)-2)'+.5)/fs;
fm=angle(wide(2:end).*conj(wide(1:end-1)))*fs/(2*pi);
active=time>=row.Start_s+.005 & time<=row.End_s-.005;
values=sort(fm(active));
assert(~isempty(values),'SCIONX:Postamble:ShortBurst','No usable discriminator samples.');
lo=values(max(1,round(.1*numel(values)))); hi=values(max(1,round(.9*numel(values))));
offset=(lo+hi)/2;
z=z.*exp(-2i*pi*offset*(0:numel(z)-1)'/fs);
filtered=conv(z,lowpass(options.ReceiverLowpassHz,fs),'same');
fm=angle(filtered(2:end).*conj(filtered(1:end-1)))*fs/(2*pi);
fm=conv(fm,[.25;.5;.25],'same');
end
values=sort(fm(active));
tones=values(max(1,round([.1 .9]*numel(values))));
midpoint=mean(tones); deviation=max(diff(tones)/2,100);
signal=struct('Time_s',time,'DiscriminatorHz',fm,'CarrierOffsetHz',offset, ...
    'ToneMidpointHz',midpoint,'DeviationHz',deviation,'SampleRate',fs,'RFStart_s',row.Start_s,'RFEnd_s',row.End_s);
record=struct('AnchorTime_s',0,'EndTime_s',0,'FlagCount',0,'SymbolRate',0, ...
    'LineCode',"",'BitPolarity',0,'ThresholdHz',0,'FitCorrelation',0, ...
    'FragmentCount',1,'FragmentAnchors_s',"",'FragmentFlagCounts',"", ...
    'FragmentSymbolRates',"",'GroupDuration_ms',0,'MaxFragmentGap_ms',0);
found=repmat(record,0,1);
if options.SignalOnly
    candidates=struct2table(found);
    return;
end
flag=[0;1;1;1;1;1;1;0];
integral=cumtrapz(time,fm);
for rate=options.SymbolRates
    for phase=0:options.PhaseCount-1
        origin=row.Start_s-2/rate+phase/options.PhaseCount/rate;
        edges=(origin:1/rate:row.End_s+1/rate)';
        edges=edges(edges>=time(1) & edges<=time(end));
        y=diff(interp1(time,integral,edges,'linear'))*rate;
        for threshold=midpoint+[-.15 0 .15]*deviation
            levels=y>threshold;
            for model=1:3
                if model==3 && ~any(options.LineCodes=="NRZI"), continue; end
                if model~=3 && ~any(options.LineCodes=="NRZ"), continue; end
                if model==3
                    bits=levels(2:end)==levels(1:end-1); code="NRZI"; polarity=1; skip=1;
                else
                    bits=levels;
                    if model==2, bits=~bits; end
                    code="NRZ"; polarity=3-2*model; skip=0;
                end
                hits=find(conv(double(bits),flip(2*flag-1),'valid')==sum(flag));
                if isempty(hits), continue; end
                a=[1;find(diff(hits)~=8)+1]; b=[a(2:end)-1;numel(hits)];
                for group=1:numel(a)
                    count=b(group)-a(group)+1;
                    if count<options.MinimumFlags, continue; end
                    idx=hits(a(group))+skip; lengthBits=8*count;
                    anchor=edges(idx);
                    if (anchor-row.Start_s)*1000<options.MinimumOffset_ms, continue; end
                    if idx+lengthBits-1>numel(y), continue; end
                    expected=2*repmat(flag,count,1)-1;
                    if model==3
                        expected=cumprod(2*repmat(flag,count,1)-1);
                    end
                    actual=y(idx:idx+lengthBits-1);
                    coefficient=pearson(actual,expected);
                    if model==3, coefficient=abs(coefficient); else, coefficient=coefficient*polarity; end
                    if coefficient<options.MinimumFitCorrelation, continue; end
                    r=record; r.AnchorTime_s=anchor; r.EndTime_s=edges(idx+lengthBits);
                    r.FlagCount=count; r.SymbolRate=rate; r.LineCode=code; r.BitPolarity=polarity;
                    r.ThresholdHz=threshold; r.FitCorrelation=coefficient;
                    r.FragmentAnchors_s=string(sprintf('%.9f',anchor));
                    r.FragmentFlagCounts=string(count); r.FragmentSymbolRates=string(rate);
                    r.GroupDuration_ms=(r.EndTime_s-r.AnchorTime_s)*1000;
                    found(end+1,1)=r; %#ok<AGROW>
                end
            end
        end
    end
end
candidates=struct2table(found);
if isempty(candidates), return; end
% Collapse phase/clock/threshold alternatives from the SAME flag group.
candidates=sortrows(candidates,{'FlagCount','FitCorrelation'},{'descend','descend'});
keep=false(height(candidates),1);
for j=1:height(candidates)
    old=find(keep & candidates.LineCode==candidates.LineCode(j) & ...
        candidates.BitPolarity==candidates.BitPolarity(j));
    overlap=min(candidates.EndTime_s(j),candidates.EndTime_s(old))- ...
        max(candidates.AnchorTime_s(j),candidates.AnchorTime_s(old));
    if any(overlap>0), continue; end
    keep(j)=true;
end
candidates=candidates(keep,:);
% Weak/clock-shifted flags can form two exact runs separated by a bad flag.
% Associate nearby runs, preserving their anchors and independent clock fits.
% No missing bits are filled in; FlagCount counts only the matched fragments.
candidates=sortrows(candidates,{'LineCode','BitPolarity','AnchorTime_s'});
j=1;
while j<height(candidates)
    gap=(candidates.AnchorTime_s(j+1)-candidates.EndTime_s(j))*1000;
    sameModel=candidates.LineCode(j)==candidates.LineCode(j+1) && ...
        candidates.BitPolarity(j)==candidates.BitPolarity(j+1);
    if sameModel && gap>=0 && gap<=options.MaximumFragmentGap_ms
        candidates.EndTime_s(j)=candidates.EndTime_s(j+1);
        candidates.FlagCount(j)=candidates.FlagCount(j)+candidates.FlagCount(j+1);
        candidates.FitCorrelation(j)=min(candidates.FitCorrelation(j:j+1));
        candidates.FragmentCount(j)=candidates.FragmentCount(j)+candidates.FragmentCount(j+1);
        candidates.FragmentAnchors_s(j)=candidates.FragmentAnchors_s(j)+";"+candidates.FragmentAnchors_s(j+1);
        candidates.FragmentFlagCounts(j)=candidates.FragmentFlagCounts(j)+";"+candidates.FragmentFlagCounts(j+1);
        candidates.FragmentSymbolRates(j)=candidates.FragmentSymbolRates(j)+";"+candidates.FragmentSymbolRates(j+1);
        candidates.MaxFragmentGap_ms(j)=max([candidates.MaxFragmentGap_ms(j:j+1);gap]);
        candidates.GroupDuration_ms(j)=(candidates.EndTime_s(j)-candidates.AnchorTime_s(j))*1000;
        candidates(j+1,:)=[];
    else
        j=j+1;
    end
end
candidates=sortrows(candidates,{'FlagCount','FitCorrelation'},{'descend','descend'});
end

function h=lowpass(cutoff,fs)
% Approximately the same impulse support at 50/100 kS/s; zero phase FIR.
half=max(16,round(64*fs/100000)); m=(-half:half)'; q=2*cutoff/fs*m;
s=ones(size(q)); mask=q~=0; s(mask)=sin(pi*q(mask))./(pi*q(mask));
h=2*cutoff/fs*s.*(.54+.46*cos(pi*m/half)); h=h/sum(h);
end

function value=pearson(a,b)
a=a-mean(a); b=b-mean(b);
value=sum(a.*b)/max(sqrt(sum(a.^2)*sum(b.^2)),eps);
end
