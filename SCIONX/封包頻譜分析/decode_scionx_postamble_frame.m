function result=decode_scionx_postamble_frame(row,signal,options)
% Search receiver timing around the detected first closing flag. No bit repair.
arguments
    row table
    signal struct
    options.SymbolRates (1,:) double = 9590:5:9635
    options.PhaseOffsets (1,:) double = -.5:1/32:.5
    options.ExpectedFrameBytes (1,1) double = 274
end
t=signal.Time_s(:)-row.AnchorTime_s; y=signal.DiscriminatorHz(:);
integral=cumtrapz(t,y); threshold=signal.ToneMidpointHz; deviation=signal.DeviationHz;
if ~isfinite(deviation), error('Invalid tone deviation'); end
rfStart=row.RFStart_s-row.AnchorTime_s;
flag=logical([0;1;1;1;1;1;1;0]); bestKey=[-Inf -Inf -Inf];
result=struct('Found',false,'Candidates',table(),'Frame',[], ...
    'Bits',[],'Edges_s',[],'SymbolRate',NaN,'PhaseOffset',NaN, ...
    'ThresholdHz',NaN,'Polarity',NaN,'ClosingTime_s',NaN,'FlagCount',0);
records=zeros(0,9);
for rate=options.SymbolRates
    for phase=options.PhaseOffsets
        anchor=phase/rate;
        first=ceil((max(rfStart,t(1))-anchor)*rate);
        last=min(48,floor((min(row.RFEnd_s-row.AnchorTime_s,t(end))-anchor)*rate));
        edges=anchor+(first:last)'/rate; closing=1-first;
        if closing<1 || closing+23>numel(edges)-1, continue; end
        levels=diff(interp1(t,integral,edges,'linear'))*rate;
        for offset=[0 -.15 .15]
            thr=threshold+offset*deviation;
            for polarity=[1 -1]
                bits=polarity*(levels-thr)>0;
                count=0;
                for j=closing:8:numel(bits)-7
                    if ~isequal(bits(j:j+7),flag), break; end
                    count=count+1;
                end
                if count<3, continue; end
                frame=scionx_reverse_hdlc(bits(1:closing-1));
                n=numel(frame.Bytes);
                if n<2 || n>options.ExpectedFrameBytes, continue; end
                % CRC-valid complete frames dominate. Otherwise rank by
                % observed flag waveform fit, not maximum recoverable length.
                f=levels(closing:closing+8*count-1)-thr;
                expected=polarity*(2*double(repmat(flag,count,1))-1);
                fit=sum(f.*expected)/max(sqrt(sum(f.^2)*sum(expected.^2)),eps);
                key=[double(frame.CRC=="PASS") fit n];
                records(end+1,:)=[rate phase thr polarity count n frame.Complete frame.CRC=="PASS" fit]; %#ok<AGROW>
                if better(key,bestKey)
                    bestKey=key; result.Found=true; result.Frame=frame;
                    result.Bits=bits(1:closing-1); result.Edges_s=edges(1:closing)+row.AnchorTime_s;
                    result.SymbolRate=rate; result.PhaseOffset=phase; result.ThresholdHz=thr;
                    result.Polarity=polarity; result.ClosingTime_s=row.AnchorTime_s+anchor;
                    result.FlagCount=count;
                end
            end
        end
    end
end
result.Candidates=array2table(records,'VariableNames',{'SymbolRate','PhaseOffset', ...
    'ThresholdHz','Polarity','FlagCount','ByteCount','Complete','CRCPass','FlagFit'});
end

function yes=better(a,b)
j=find(a~=b,1); yes=~isempty(j) && a(j)>b(j);
end
