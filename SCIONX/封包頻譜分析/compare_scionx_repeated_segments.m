function result = compare_scionx_repeated_segments(alignment, options)
%COMPARE_SCIONX_REPEATED_SEGMENTS Overlay post-flag data with earlier IQ FM.
% Input: align_scionx_postambles result or its Postamble_Alignment.mat file.
% Correlations use the saved discriminator BEFORE additional plot smoothing.
% Same-burst time translations only; no polarity/gain fitting or time warping.
arguments
    alignment
    options.OutputDirectory (1,1) string = ""
    options.AfterWindow_ms (1,2) double {mustBeFinite} = [8.5 49.5]
    options.LagRange_ms (1,2) double {mustBePositive,mustBeFinite} = [50 180]
    options.UseSourcePrior (1,1) logical = true
    options.PreferredSourceStart_ms (1,1) double {mustBeFinite} = -100
    options.SourcePriorRadius_ms (1,1) double {mustBePositive,mustBeFinite} = 10
    options.FractionalSteps (1,1) double {mustBeInteger,mustBePositive} = 32
    options.SourceGuard_ms (1,1) double {mustBeNonnegative,mustBeFinite} = .2
    options.ReceiverLowpassHz (:,1) double {mustBePositive,mustBeFinite} = []
    options.MinimumCorrelation (1,1) double {mustBeNonnegative,mustBeFinite} = .8
    options.MinimumModulationStd_Hz (1,1) double {mustBeNonnegative,mustBeFinite} = 300
    options.LocalWindow_ms (1,1) double {mustBePositive,mustBeFinite} = 5
    options.LocalHop_ms (1,1) double {mustBePositive,mustBeFinite} = 2
    options.LocalLagRadius_ms (1,1) double {mustBePositive,mustBeFinite} = 1.2
    options.StableLagTolerance_ms (1,1) double {mustBePositive,mustBeFinite} = .04
    options.MinimumStableWindows (1,1) double {mustBeInteger,mustBePositive} = 3
    options.DisplaySmoothing_ms (1,1) double {mustBeNonnegative,mustBeFinite} = .1
    options.FrequencyLimitHz (1,1) double {mustBePositive,mustBeFinite} = 5000
    options.IncludeReference (1,1) logical = true
    options.RowsPerPage (1,1) double {mustBePositive,mustBeInteger} = 4
    options.ExportPlots (1,1) logical = true
    options.OpenViewer (1,1) logical = true
end
sourceMat="";
if ischar(alignment) || isstring(alignment)
    sourceMat=string(alignment); loaded=load(sourceMat,'result'); alignment=loaded.result;
end
assert(isstruct(alignment) && all(isfield(alignment, ...
    {'AlignmentResults','Signals','ReferenceRow','ComparisonRows'})), ...
    'SCIONX:Repeat:Input','Provide a postamble alignment result or its MAT file.');
assert(options.AfterWindow_ms(1)>=0 && diff(options.AfterWindow_ms)>0 && ...
    diff(options.LagRange_ms)>0 && options.MinimumCorrelation<=1, ...
    'SCIONX:Repeat:Configuration','Invalid target window, lag range or correlation threshold.');
assert(options.LocalWindow_ms<=diff(options.AfterWindow_ms) && ...
    options.LocalHop_ms<=options.LocalWindow_ms, ...
    'SCIONX:Repeat:Window','Local window/hop must fit inside the target window.');
if strlength(options.OutputDirectory)==0
    options.OutputDirectory=fullfile(fileparts(mfilename('fullpath')),'Repeated_Segment_Output');
end
if ~isfolder(options.OutputDirectory), mkdir(options.OutputDirectory); end
indices=alignment.ComparisonRows(:);
if options.IncludeReference, indices=[alignment.ReferenceRow;indices]; end
indices=unique(indices,'stable');
if ~isempty(options.ReceiverLowpassHz)
    assert(isscalar(options.ReceiverLowpassHz) || ...
        numel(options.ReceiverLowpassHz)>=max(alignment.AlignmentResults.FileIndex), ...
        'SCIONX:Repeat:ReceiverSettings','ReceiverLowpassHz must be scalar or indexed by input FileIndex.');
end
assert(all(isfinite(indices) & indices>=1 & indices<=height(alignment.AlignmentResults) & indices==fix(indices)), ...
    'SCIONX:Repeat:Rows','Invalid alignment row indices.');
records=repmat(emptySummary(),numel(indices),1); views=cell(numel(indices),1);
segments=struct2table(repmat(emptySegment(),0,1));
locals=struct2table(repmat(emptyLocal(),0,1));
for k=1:numel(indices)
    index=indices(k); row=alignment.AlignmentResults(index,:);
    signal=alignment.Signals{index}; cutoff=NaN;
    if isfield(alignment,'Configuration') && isfield(alignment.Configuration,'ReceiverLowpassHz')
        cutoff=alignment.Configuration.ReceiverLowpassHz;
    end
    if ~isempty(options.ReceiverLowpassHz)
        cutoff=options.ReceiverLowpassHz(min(row.FileIndex,numel(options.ReceiverLowpassHz)));
        signalRow=row;
        signalRow.StartSample=round(row.RFStart_s*row.SampleRate);
        signalRow.EndSampleExclusive=round(row.RFEnd_s*row.SampleRate);
        signalRow.Start_s=row.RFStart_s; signalRow.End_s=row.RFEnd_s;
        [~,signal]=detect_scionx_postamble(signalRow,ReceiverLowpassHz=cutoff,SignalOnly=true, ...
            Padding_ms=max(options.LagRange_ms)+max(options.AfterWindow_ms)+20);
    end
    [records(k),views{k},newSegments,newLocals]=analyzeRow(row,signal,index,options);
    records(k).ReceiverLowpassHz=cutoff;
    segments=[segments;newSegments]; locals=[locals;newLocals]; %#ok<AGROW>
    fprintf('File %d burst #%d: %s | single lag %.6f ms | r %.4f | %d stable segments\n', ...
        row.FileIndex,row.BurstID,records(k).Status,records(k).GlobalLag_ms, ...
        records(k).GlobalCorrelation,records(k).SegmentCount);
end
summary=struct2table(records);
result=struct('Configuration',options,'SourceAlignmentMat',sourceMat, ...
    'Summary',summary,'Segments',segments,'LocalWindows',locals,'Views',{views}, ...
    'SummaryCsv',fullfile(options.OutputDirectory,'Repeated_Summary.csv'), ...
    'SegmentsCsv',fullfile(options.OutputDirectory,'Repeated_Segments.csv'), ...
    'LocalWindowsCsv',fullfile(options.OutputDirectory,'Repeated_Local_Windows.csv'), ...
    'PlotFiles',strings(0,1),'DiagnosticFiles',strings(0,1), ...
    'MatFile',fullfile(options.OutputDirectory,'Repeated_Overlay.mat'));
writetable(summary,result.SummaryCsv); writetable(segments,result.SegmentsCsv);
writetable(locals,result.LocalWindowsCsv);
if options.ExportPlots && ~isempty(summary)
    [result.PlotFiles,result.DiagnosticFiles]=exportPages(result,options);
end
save(result.MatFile,'result');
if options.OpenViewer && ~isempty(summary), openViewer(result,options); end
fprintf('Saved repetition comparison: %s\n',options.OutputDirectory);
end

function r=emptySummary()
r=struct('BurstRow',0,'FileIndex',0,'BurstID',0,'FileLabel',"",'AnchorTime_s',NaN, ...
    'AnchorStatus',"",'Status',"unresolved_anchor",'SampleRate',NaN,'ReceiverLowpassHz',NaN, ...
    'AfterStart_ms',NaN,'AfterEnd_ms',NaN,'TargetModulationStd_Hz',NaN, ...
    'GlobalLag_ms',NaN,'GlobalSourceStart_ms',NaN,'GlobalSourceEnd_ms',NaN, ...
    'GlobalCorrelation',NaN,'GlobalSimilarity_pct',NaN,'GlobalRMSE_Hz',NaN, ...
    'SegmentCount',0,'SegmentCoverage_pct',0, ...
    'SearchMode',"not_searched",'GlobalSearchCandidates',0,'AnalysisSeconds',NaN);
end

function r=emptySegment()
r=struct('BurstRow',0,'SegmentID',0,'AfterStart_ms',NaN,'AfterEnd_ms',NaN, ...
    'SourceStart_ms',NaN,'SourceEnd_ms',NaN,'Lag_ms',NaN, ...
    'Correlation',NaN,'Similarity_pct',NaN,'RMSE_Hz',NaN,'SampleCount',0);
end

function r=emptyLocal()
r=struct('BurstRow',0,'Start_ms',NaN,'End_ms',NaN,'Center_ms',NaN, ...
    'Lag_ms',NaN,'Correlation',NaN,'Similarity_pct',NaN,'RMSE_Hz',NaN, ...
    'TargetModulationStd_Hz',NaN,'Status',"",'StableRun',0);
end

function [r,v,segments,locals]=analyzeRow(row,signal,index,o)
analysisTimer=tic;
r=emptySummary(); r.BurstRow=index; r.FileIndex=row.FileIndex; r.BurstID=row.BurstID;
r.FileLabel=row.FileLabel; r.AnchorTime_s=row.AnchorTime_s; r.AnchorStatus=row.Status;
r.SampleRate=signal.SampleRate;
v=struct('Time_ms',[],'AfterHz',[],'GlobalEarlierHz',[],'SegmentEarlierHz',[], ...
    'SearchLag_ms',[],'SearchCorrelation',[]);
segments=struct2table(repmat(emptySegment(),0,1)); locals=struct2table(repmat(emptyLocal(),0,1));
if ~isfinite(row.AnchorTime_s), return; end
fs=signal.SampleRate; t=(signal.Time_s(:)-row.AnchorTime_s)*1000;
y=signal.DiscriminatorHz(:);
assert(numel(t)==numel(y) && numel(t)>2 && all(isfinite(y)) && all(diff(t)>0), ...
    'SCIONX:Repeat:Signal','Invalid discriminator signal.');
assert(max(abs(diff(t)-1000/fs))<.001, ...
    'SCIONX:Repeat:SampleGrid','Discriminator sample times must be uniform.');
% Modulation gate uses 0.1-ms smoothing; correlation uses unsmoothed y.
gate=movmean(y,max(1,round(.0001*fs)));
rfStart=(row.RFStart_s-row.AnchorTime_s)*1000;
rfEnd=(row.RFEnd_s-row.AnchorTime_s)*1000;
window=o.AfterWindow_ms;
if ismember('FlagGroupDuration_ms',row.Properties.VariableNames) && isfinite(row.FlagGroupDuration_ms)
    window(1)=max(window(1),row.FlagGroupDuration_ms+.1);
end
r.AfterStart_ms=window(1); r.AfterEnd_ms=window(2);
if window(1)>=window(2) || window(2)>min(rfEnd,t(end)) || window(1)<max(rfStart,t(1))
    r.Status="target_outside_data"; return;
end
target=find(t>=window(1) & t<=window(2));
if numel(target)<8, r.Status="insufficient_samples"; return; end
v.Time_ms=t(target); v.AfterHz=y(target);
v.GlobalEarlierHz=nan(size(target)); v.SegmentEarlierHz=nan(size(target));
m=[]; searchMode="full_range"; candidates=0;
if o.UseSourcePrior
    % A source start near -100 ms means a lag near target start + 100 ms.
    preferred=v.Time_ms(1)-o.PreferredSourceStart_ms;
    priorRange=[max(o.LagRange_ms(1),preferred-o.SourcePriorRadius_ms), ...
        min(o.LagRange_ms(2),preferred+o.SourcePriorRadius_ms)];
    if diff(priorRange)>0
        m=findMatch(t,y,gate,target,rfStart,priorRange,o);
        candidates=numel(m.SearchLag_ms); searchMode="source_prior";
        % A weak or boundary solution must not lock the search onto the prior.
        edgeMargin=max(2*1000/fs,.02);
        if m.Status~="low_modulation" && (~isfinite(m.Correlation) || ...
                m.Correlation<o.MinimumCorrelation || ...
                m.Lag_ms<priorRange(1)+edgeMargin || m.Lag_ms>priorRange(2)-edgeMargin)
            m=[]; searchMode="full_range_fallback";
        end
    end
end
if isempty(m), m=findMatch(t,y,gate,target,rfStart,o.LagRange_ms,o); end
if searchMode~="source_prior", candidates=candidates+numel(m.SearchLag_ms); end
r.SearchMode=searchMode; r.GlobalSearchCandidates=candidates;
r.Status=m.Status; r.TargetModulationStd_Hz=m.TargetModulationStd_Hz;
v.SearchLag_ms=m.SearchLag_ms; v.SearchCorrelation=m.SearchCorrelation;
if ~isfinite(m.Lag_ms), r.AnalysisSeconds=toc(analysisTimer); return; end
r.GlobalLag_ms=m.Lag_ms; r.GlobalCorrelation=m.Correlation;
r.GlobalSimilarity_pct=max(0,m.Correlation)*100; r.GlobalRMSE_Hz=m.RMSE_Hz;
r.GlobalSourceStart_ms=v.Time_ms(1)-m.Lag_ms;
r.GlobalSourceEnd_ms=v.Time_ms(end)-m.Lag_ms;
v.GlobalEarlierHz=interp1(t,y,v.Time_ms-m.Lag_ms,'linear');
if m.Correlation>=o.MinimumCorrelation, r.Status="single_delay_candidate"; else, r.Status="low_correlation"; end
localRange=[max(o.LagRange_ms(1),m.Lag_ms-o.LocalLagRadius_ms), ...
    min(o.LagRange_ms(2),m.Lag_ms+o.LocalLagRadius_ms)];
starts=window(1):o.LocalHop_ms:window(2)-o.LocalWindow_ms;
if isempty(starts) || starts(end)<window(2)-o.LocalWindow_ms-1e-9
    starts(end+1)=window(2)-o.LocalWindow_ms;
end
localRecords=repmat(emptyLocal(),numel(starts),1);
for j=1:numel(starts)
    targetLocal=find(t>=starts(j) & t<=starts(j)+o.LocalWindow_ms);
    match=findMatch(t,y,gate,targetLocal,rfStart,localRange,o);
    q=emptyLocal(); q.BurstRow=index; q.Start_ms=starts(j); q.End_ms=starts(j)+o.LocalWindow_ms;
    q.Center_ms=(q.Start_ms+q.End_ms)/2; q.Lag_ms=match.Lag_ms; q.Correlation=match.Correlation;
    if isfinite(q.Correlation), q.Similarity_pct=max(0,q.Correlation)*100; end
    q.RMSE_Hz=match.RMSE_Hz; q.Status=match.Status;
    q.TargetModulationStd_Hz=match.TargetModulationStd_Hz;
    localRecords(j)=q;
end
locals=struct2table(localRecords);
good=find(locals.Correlation>=o.MinimumCorrelation & isfinite(locals.Lag_ms));
groups=cell(0,1);
for j=good(:)'
    if isempty(groups) || locals.Center_ms(j)-locals.Center_ms(groups{end}(end))>2*o.LocalWindow_ms || ...
            abs(locals.Lag_ms(j)-median(locals.Lag_ms(groups{end})))>o.StableLagTolerance_ms
        groups{end+1,1}=j; %#ok<AGROW>
    else
        groups{end}(end+1)=j;
    end
end
groups=groups(cellfun(@numel,groups)>=o.MinimumStableWindows);
for j=1:numel(groups)
    g=groups{j}; bounds=[locals.Center_ms(g(1)),locals.Center_ms(g(end))];
    if g(1)==1, bounds(1)=window(1); end
    if g(end)==height(locals), bounds(2)=window(2); end
    center=median(locals.Lag_ms(g));
    range=[max(localRange(1),center-o.StableLagTolerance_ms),min(localRange(2),center+o.StableLagTolerance_ms)];
    targetSegment=find(t>=bounds(1) & t<=bounds(2));
    match=findMatch(t,y,gate,targetSegment,rfStart,range,o);
    if ~isfinite(match.Correlation) || match.Correlation<o.MinimumCorrelation, continue; end
    q=emptySegment(); q.BurstRow=index; q.SegmentID=height(segments)+1;
    q.AfterStart_ms=t(targetSegment(1)); q.AfterEnd_ms=t(targetSegment(end));
    q.Lag_ms=match.Lag_ms; q.SourceStart_ms=q.AfterStart_ms-q.Lag_ms; q.SourceEnd_ms=q.AfterEnd_ms-q.Lag_ms;
    q.Correlation=match.Correlation; q.Similarity_pct=max(0,q.Correlation)*100;
    q.RMSE_Hz=match.RMSE_Hz; q.SampleCount=numel(targetSegment);
    segments=[segments;struct2table(q)]; %#ok<AGROW>
    locals.StableRun(g)=q.SegmentID;
    hit=v.Time_ms>=q.AfterStart_ms & v.Time_ms<=q.AfterEnd_ms;
    v.SegmentEarlierHz(hit)=interp1(t,y,v.Time_ms(hit)-q.Lag_ms,'linear');
end
r.SegmentCount=height(segments);
r.SegmentCoverage_pct=100*sum(isfinite(v.SegmentEarlierHz))/numel(v.Time_ms);
if r.SegmentCount>1 && max(segments.Lag_ms)-min(segments.Lag_ms)>o.StableLagTolerance_ms
    r.Status="piecewise_delay_candidate";
elseif r.SegmentCount==1 && r.GlobalCorrelation<o.MinimumCorrelation
    r.Status="partial_match_candidate";
end
r.AnalysisSeconds=toc(analysisTimer);
end

function m=findMatch(t,y,gate,target,rfStart,lagRange,o)
m=struct('Status',"no_valid_source",'Lag_ms',NaN,'Correlation',NaN,'RMSE_Hz',NaN, ...
    'TargetModulationStd_Hz',NaN,'SearchLag_ms',[],'SearchCorrelation',[]);
if numel(target)<8, m.Status="insufficient_samples"; return; end
n=numel(target); step=median(diff(t));
m.TargetModulationStd_Hz=std(gate(target));
if m.TargetModulationStd_Hz<o.MinimumModulationStd_Hz
    m.Status="low_modulation"; return;
end
% Source interval must lie entirely within THIS burst, before first flag.
earliest=max([t(target(1))-lagRange(2),rfStart+.2,t(1)]);
latest=min([t(target(1))-lagRange(1),-o.SourceGuard_ms-(n-1)*step,t(end)-(n-1)*step]);
first=max(1,ceil((earliest-t(1))/step)+1);
last=min(numel(t)-n+1,floor((latest-t(1))/step)+1);
if last<first, return; end
source=y(first:last+n-1); targetY=y(target); centered=targetY-mean(targetY);
fftLength=2^nextpow2(numel(source)+n-1);
convolution=real(ifft(fft(source,fftLength).*fft(flipud(centered),fftLength)));
numerator=convolution(n:numel(source));
sum1=[0;cumsum(source)]; sum2=[0;cumsum(source.^2)];
variance=max(0,sum2(n+1:end)-sum2(1:end-n)-(sum1(n+1:end)-sum1(1:end-n)).^2/n);
denominator=sqrt(sum(centered.^2)*variance);
score=numerator./max(denominator,eps);
g=gate(first:last+n-1); g1=[0;cumsum(g)]; g2=[0;cumsum(g.^2)];
gvar=max(0,g2(n+1:end)-g2(1:end-n)-(g1(n+1:end)-g1(1:end-n)).^2/n);
score(sqrt(gvar/max(1,n-1))<o.MinimumModulationStd_Hz | denominator<=eps)=NaN;
invalid=~isfinite(score);
score=max(-1,min(1,score));
score(invalid)=NaN;
lags=(target(1)-(first:last)')*step;
m.SearchLag_ms=flipud(lags); m.SearchCorrelation=flipud(score);
valid=find(isfinite(score)); if isempty(valid), return; end
[~,pick]=max(score(valid)); pick=valid(pick);
bestLag=lags(pick); bestScore=score(pick);
bestSource=y(first+pick-1:first+pick+n-2);
% Uniform sample grid: interpolate adjacent source samples directly, avoiding
% a full time-vector search for every fractional lag and local window.
for lag=bestLag+(-1:1/o.FractionalSteps:1)*step
    times=t(target)-lag;
    if lag<lagRange(1) || lag>lagRange(2) || times(1)<max(rfStart+.2,t(1)) || ...
            times(end)>-o.SourceGuard_ms || times(end)>t(end), continue; end
    sampleShift=lag/step;
    positions=target-sampleShift;
    left=floor(positions); fraction=positions-left;
    if any(left<1 | left>=numel(y)), continue; end
    candidate=y(left).*(1-fraction)+y(left+1).*fraction;
    correlation=pearson(targetY,candidate);
    if correlation>bestScore
        bestScore=correlation; m.Lag_ms=lag; bestSource=candidate;
    end
end
if ~isfinite(m.Lag_ms), m.Lag_ms=bestLag; end
m.Correlation=bestScore; m.RMSE_Hz=sqrt(mean((targetY-bestSource).^2)); m.Status="searched";
end

function value=pearson(a,b)
a=a-mean(a); b=b-mean(b); denominator=sqrt(sum(a.^2)*sum(b.^2));
if denominator<=eps, value=NaN; else, value=max(-1,min(1,sum(a.*b)/denominator)); end
end

function plotOverlay(ax,r,v,segments,o)
hold(ax,'on'); colorBlue=[.05 .35 .7]; colorOrange=[.9 .42 .05];
span=[r.AfterStart_ms r.AfterEnd_ms];
if any(~isfinite(span)) || diff(span)<=0, span=o.AfterWindow_ms; end
if isempty(v.Time_ms)
    text(ax,.5,.5,char(r.Status),'Units','normalized','HorizontalAlignment','center');
else
    n=max(1,round(o.DisplaySmoothing_ms/1000*r.SampleRate));
    plot(ax,v.Time_ms,movmean(v.AfterHz,n)/1000,'Color',colorBlue,'DisplayName','After flags');
    if r.SegmentCount>0
        missing=~isfinite(v.SegmentEarlierHz); d=diff([false;missing(:);false]);
        starts=find(d==1); ends=find(d==-1)-1;
        for j=1:numel(starts)
            a=v.Time_ms(starts(j)); b=v.Time_ms(ends(j)); lim=[-1 1]*o.FrequencyLimitHz/1000;
            patch(ax,[a b b a],[lim(1) lim(1) lim(2) lim(2)],[.88 .9 .93], ...
                'FaceAlpha',.35,'EdgeColor','none','HandleVisibility','off');
        end
        for j=1:height(segments)
            hit=v.Time_ms>=segments.AfterStart_ms(j) & v.Time_ms<=segments.AfterEnd_ms(j);
            plot(ax,v.Time_ms(hit),movmean(v.SegmentEarlierHz(hit),n)/1000, ...
                'Color',colorOrange,'DisplayName',sprintf('Earlier S%d',j));
            label=sprintf('S%d: %.4f ms | r=%.3f (%.1f%%)', ...
                j,segments.Lag_ms(j),segments.Correlation(j),segments.Similarity_pct(j));
            text(ax,mean([segments.AfterStart_ms(j) segments.AfterEnd_ms(j)]), ...
                o.FrequencyLimitHz/1000*.9,label,'HorizontalAlignment','center','FontSize',8, ...
                'BackgroundColor','w','Margin',1);
        end
    elseif isfinite(r.GlobalLag_ms)
        plot(ax,v.Time_ms,movmean(v.GlobalEarlierHz,n)/1000,'Color',colorOrange, ...
            'DisplayName','Earlier, single lag');
    end
    legend(ax,'Location','southoutside','Orientation','horizontal','FontSize',8);
end
grid(ax,'on'); xlim(ax,span); ylim(ax,[-1 1]*o.FrequencyLimitHz/1000);
xlabel(ax,'Time after first postamble flag (ms)'); ylabel(ax,'Frequency (kHz)');
if isfinite(r.GlobalCorrelation)
    detail=sprintf('Single lag %.5f ms | r=%.4f (%.1f%%) | RMS difference %.3f kHz | stable coverage %.1f%%', ...
        r.GlobalLag_ms,r.GlobalCorrelation,r.GlobalSimilarity_pct,r.GlobalRMSE_Hz/1000,r.SegmentCoverage_pct);
else
    detail=char(r.Status+" | insufficient changing waveform for a reliable overlap score");
end
title(ax,sprintf('%s | burst #%d | %s\n%s',r.FileLabel,r.BurstID,r.Status,detail), ...
    'Interpreter','none','FontSize',9);
end

function plotLocal(ax,r,locals,o)
hold(ax,'on');
if isempty(locals) || ~any(isfinite(locals.Lag_ms))
    text(ax,.5,.5,char(r.Status),'Units','normalized','HorizontalAlignment','center');
else
    scatter(ax,locals.Center_ms,locals.Lag_ms,32,locals.Correlation,'filled');
    colormap(ax,turbo); clim(ax,[0 1]); cb=colorbar(ax); cb.Label.String='Local Pearson r';
    yline(ax,r.GlobalLag_ms,'k--','single lag');
    center=median(locals.Lag_ms,'omitnan'); spread=max(.2,max(abs(locals.Lag_ms-center),[],'omitnan')+.03);
    ylim(ax,center+[-spread spread]);
end
xlim(ax,o.AfterWindow_ms); grid(ax,'on');
xlabel(ax,'Local window center after first flag (ms)'); ylabel(ax,'Best local delay (ms)');
title(ax,sprintf('%.1f ms windows | local search +/- %.2f ms around global lag', ...
    o.LocalWindow_ms,o.LocalLagRadius_ms),'FontSize',9);
end

function [paths,diagnostics]=exportPages(result,o)
paths=strings(0,1); diagnostics=strings(0,1);
for first=1:o.RowsPerPage:height(result.Summary)
    indices=first:min(first+o.RowsPerPage-1,height(result.Summary));
    fig=figure('Visible','off','Color','w','Position',[50 50 1500 340*numel(indices)]);
    cleanup=onCleanup(@() close(fig));
    layout=tiledlayout(fig,numel(indices),1,'TileSpacing','compact','Padding','compact');
    title(layout,'SCION-X | earlier waveform overlaid on postamble-following data','FontSize',13);
    subtitle(layout,'Blue: received after flags. Orange: earlier waveform shifted forward. Percent = positive Pearson r x 100; not bit agreement.');
    for index=indices
        r=result.Summary(index,:); segments=result.Segments(result.Segments.BurstRow==r.BurstRow,:);
        plotOverlay(nexttile(layout),r,result.Views{index},segments,o);
    end
    path=fullfile(o.OutputDirectory,sprintf('Repeated_Segments_%03d.png',ceil(first/o.RowsPerPage)));
    drawnow; set(fig,'PaperPositionMode','auto'); print(fig,path,'-dpng','-r120');
    paths(end+1,1)=path; %#ok<AGROW>
    clear cleanup
    fig=figure('Visible','off','Color','w','Position',[50 50 1500 280*numel(indices)]);
    cleanup=onCleanup(@() close(fig));
    layout=tiledlayout(fig,numel(indices),2,'TileSpacing','compact','Padding','compact');
    title(layout,'SCION-X | delay search and local timing changes','FontSize',13);
    for index=indices
        r=result.Summary(index,:); v=result.Views{index}; ax=nexttile(layout);
        plot(ax,v.SearchLag_ms,v.SearchCorrelation); grid(ax,'on');
        xlim(ax,o.LagRange_ms); ylim(ax,[-1 1]); xlabel(ax,'Earlier-to-later delay (ms)'); ylabel(ax,'Pearson r');
        if isfinite(r.GlobalLag_ms), xline(ax,r.GlobalLag_ms,'r--'); end
        title(ax,sprintf('%s | #%d | %s',r.FileLabel,r.BurstID,r.SearchMode),'Interpreter','none','FontSize',8);
        plotLocal(nexttile(layout),r,result.LocalWindows(result.LocalWindows.BurstRow==r.BurstRow,:),o);
    end
    path=fullfile(o.OutputDirectory,sprintf('Repeated_Delay_%03d.png',ceil(first/o.RowsPerPage)));
    drawnow; set(fig,'PaperPositionMode','auto'); print(fig,path,'-dpng','-r120');
    diagnostics(end+1,1)=path; %#ok<AGROW>
    clear cleanup
end
end

function openViewer(result,o)
fig=figure('Name','SCION-X Repeated Segment Overlay','NumberTitle','off', ...
    'Color','w','Position',[70 70 1450 850]);
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
layout.OuterPosition=[0 .1 1 .9];
labels=strings(height(result.Summary),1);
for j=1:numel(labels)
    r=result.Summary(j,:); labels(j)=sprintf('%d/%d | #%d | %s | %s',j,numel(labels),r.BurstID,r.Status,r.FileLabel);
end
selection=min(2,numel(labels));
menu=uicontrol(fig,'Style','popupmenu','Units','normalized','Position',[.19 .03 .63 .05], ...
    'String',cellstr(labels),'Value',selection,'Callback',@select);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.03 .03 .13 .05], ...
    'String','Previous','Callback',@(~,~) move(-1));
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.85 .03 .12 .05], ...
    'String','Next','Callback',@(~,~) move(1));
draw();
    function select(src,~)
        selection=src.Value; draw();
    end
    function move(delta)
        selection=max(1,min(numel(labels),selection+delta)); menu.Value=selection; draw();
    end
    function draw()
        delete(layout.Children); r=result.Summary(selection,:);
        plotOverlay(nexttile(layout),r,result.Views{selection},result.Segments(result.Segments.BurstRow==r.BurstRow,:),o);
        plotLocal(nexttile(layout),r,result.LocalWindows(result.LocalWindows.BurstRow==r.BurstRow,:),o);
        drawnow;
    end
end
