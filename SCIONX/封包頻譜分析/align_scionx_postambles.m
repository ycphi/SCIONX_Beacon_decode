function result = align_scionx_postambles(scanResult, options)
%ALIGN_SCIONX_POSTAMBLES Align bursts on first closing-flag candidates.
% Input: result from arrange_scionx_bursts, or its Burst_Results.mat path.
% Works with mixed CS16/CF32 sample rates. See README_Postamble_Alignment.
arguments
    scanResult
    options.OutputDirectory (1,1) string = ""
    options.ReferencePostambleTime_s (1,1) double = NaN
    options.ManualPostambleTimes_s (:,1) double = NaN
    options.SymbolRates (1,:) double {mustBePositive,mustBeFinite} = 9590:5:9635
    options.LineCodes (1,:) string {mustBeMember(options.LineCodes,["NRZ","NRZI"])} = "NRZ"
    options.MinimumFlags (1,1) double {mustBeInteger,mustBePositive} = 4
    options.MinimumOffset_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 40
    options.MinimumFitCorrelation (1,1) double {mustBeNonnegative,mustBeFinite} = .65
    options.MaximumFragmentGap_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 1.5
    options.AmbiguityMargin (1,1) double {mustBeNonnegative,mustBeFinite} = .05
    options.CandidateFlagCountTolerance (1,1) double {mustBeNonnegative,mustBeInteger} = 1
    options.ResolveRepeatedPostambles (1,1) logical = false
    options.ReceiverLowpassHz (1,1) double {mustBePositive,mustBeFinite} = 8000
    options.Smoothing_ms (1,1) double {mustBeNonnegative,mustBeFinite} = .5
    options.DetailSmoothing_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 0
    options.Window_ms (1,2) double {mustBeFinite} = [-260 210]
    options.DetailWindow_ms (1,2) double {mustBeFinite} = [-5 15]
    options.FrequencyLimitHz (1,1) double {mustBePositive,mustBeFinite} = 5000
    options.ShadeProjectedPrefix (1,1) logical = true
    options.BurstsPerPage (1,1) double {mustBeInteger,mustBePositive} = 4
    options.ExportPlots (1,1) logical = true
    options.OpenViewer (1,1) logical = true
end
if ischar(scanResult) || isstring(scanResult)
    loaded=load(scanResult,'result'); scanResult=loaded.result;
end
assert(isstruct(scanResult) && isfield(scanResult,'BurstResults') && ...
    isfield(scanResult,'ReferenceRow'),'SCIONX:Postamble:Input','Provide a burst scan result or MAT file.');
bursts=scanResult.BurstResults; refIndex=scanResult.ReferenceRow;
assert(isfinite(refIndex) && refIndex>=1 && refIndex<=height(bursts), ...
    'SCIONX:Postamble:Reference','A detected reference burst is required.');
assert(all(ismember({'InputFormat','SampleRate','DisplayPolarity'},bursts.Properties.VariableNames)), ...
    'SCIONX:Postamble:Metadata','Re-run the current burst scanner to obtain IQ metadata.');
assert(options.Window_ms(1)<0 && options.Window_ms(2)>0 && ...
    options.DetailWindow_ms(1)<0 && options.DetailWindow_ms(2)>0, ...
    'SCIONX:Postamble:Window','Alignment windows must straddle zero.');
manual=options.ManualPostambleTimes_s;
assert(isscalar(manual) || numel(manual)==height(bursts),'SCIONX:Postamble:ManualTimes', ...
    'ManualPostambleTimes_s must be scalar or one time/NaN per BurstResults row.');
if isscalar(manual), manual=repmat(manual,height(bursts),1); end
if ~isnan(options.ReferencePostambleTime_s), manual(refIndex)=options.ReferencePostambleTime_s; end
assert(all(isnan(manual) | (isfinite(manual) & manual>=0)), ...
    'SCIONX:Postamble:ManualTimes','Manual anchors must be finite nonnegative times or NaN.');
if strlength(options.OutputDirectory)==0
    options.OutputDirectory=fullfile(fileparts(mfilename('fullpath')),'Postamble_Alignment_Output');
end
if ~isfolder(options.OutputDirectory), mkdir(options.OutputDirectory); end
records=repmat(emptyRecord(),height(bursts),1); signals=cell(height(bursts),1);
allCandidates=table();
for j=1:height(bursts)
    row=bursts(j,:);
    [c,signals{j}]=detect_scionx_postamble(row,SymbolRates=options.SymbolRates,LineCodes=options.LineCodes, ...
        MinimumFlags=options.MinimumFlags,MinimumOffset_ms=options.MinimumOffset_ms, ...
        MinimumFitCorrelation=options.MinimumFitCorrelation,ReceiverLowpassHz=options.ReceiverLowpassHz, ...
        MaximumFragmentGap_ms=options.MaximumFragmentGap_ms, ...
        Padding_ms=max(abs(options.Window_ms))+20);
    r=emptyRecord(); r.BurstRow=j; r.FileIndex=row.FileIndex; r.BurstID=row.BurstID;
    r.FilePath=row.FilePath; r.FileLabel=row.FileLabel; r.SampleRate=row.SampleRate;
    r.InputFormat=row.InputFormat; r.RFStart_s=row.Start_s;
    r.RFEnd_s=row.End_s; r.DisplayPolarity=row.DisplayPolarity;
    r.CarrierOffsetHz=signals{j}.CarrierOffsetHz;
    [automatic,selectionReason,c]=select_scionx_postamble_candidate(c,signals{j},row.End_s, ...
        AmbiguityMargin=options.AmbiguityMargin,FlagCountTolerance=options.CandidateFlagCountTolerance, ...
        ResolveRepeatedPostambles=options.ResolveRepeatedPostambles);
    selected=[]; r.SelectionReason=selectionReason;
    if ~isnan(manual(j))
        assert(manual(j)>=row.Start_s && manual(j)<row.End_s, ...
            'SCIONX:Postamble:ManualRange','Manual anchor outside RF burst row %d.',j);
        r.AnchorTime_s=manual(j); r.Status="manual"; r.SelectionReason="manual_override";
        if ~isempty(c)
            [distance,selected]=min(abs(c.AnchorTime_s-manual(j)));
            if distance>.010, selected=[]; end
        end
    elseif ~isempty(c)
        if isempty(automatic)
            r.Status="ambiguous";
        else
            selected=automatic; r.AnchorTime_s=c.AnchorTime_s(selected); r.Status="detected_candidate";
            if c.FragmentCount(selected)>1, r.Status="fragmented_candidate"; end
            if selectionReason=="unique_postflag_tone", r.Status="tone_supported_candidate"; end
            if selectionReason=="verified_repeated_postamble", r.Status="repeat_supported_candidate"; end
        end
    end
    if ~isempty(selected)
        r.DetectedAnchorTime_s=c.AnchorTime_s(selected);
        r.FlagCount=c.FlagCount(selected); r.SymbolRate=c.SymbolRate(selected);
        r.LineCode=c.LineCode(selected); r.BitPolarity=c.BitPolarity(selected);
        r.FitCorrelation=c.FitCorrelation(selected);
        r.FragmentCount=c.FragmentCount(selected);
        r.FlagGroupDuration_ms=c.GroupDuration_ms(selected);
        r.FragmentAnchors_s=c.FragmentAnchors_s(selected);
        r.MaxFragmentGap_ms=c.MaxFragmentGap_ms(selected);
    end
    if isfinite(r.AnchorTime_s)
        r.PostambleOffset_ms=(r.AnchorTime_s-row.Start_s)*1000;
        r.AnchorSampleApprox=round(r.AnchorTime_s*row.SampleRate);
        r.RFStartRelative_ms=(row.Start_s-r.AnchorTime_s)*1000;
        r.RFEndRelative_ms=(row.End_s-r.AnchorTime_s)*1000;
    end
    if ~isempty(c)
        c.BurstRow=repmat(j,height(c),1); c.FileIndex=repmat(row.FileIndex,height(c),1);
        c.BurstID=repmat(row.BurstID,height(c),1);
        c.Selected=false(height(c),1);
        if ~isempty(selected), c.Selected(selected)=true; end
        if isempty(allCandidates), allCandidates=c; else, allCandidates=[allCandidates;c]; end %#ok<AGROW>
    end
    records(j)=r;
    fprintf('File %d burst #%d: %s | anchor %.9f s | %d flags | fit %.3f\n', ...
        r.FileIndex,r.BurstID,r.Status,r.AnchorTime_s,r.FlagCount,r.FitCorrelation);
end
rows=struct2table(records);
csv=fullfile(options.OutputDirectory,'Postamble_Alignment.csv');
candidateCsv=fullfile(options.OutputDirectory,'Postamble_Candidates.csv');
if isfinite(rows.AnchorTime_s(refIndex))
    rows.AdvanceRelativeToReference_ms=rows.PostambleOffset_ms(refIndex)-rows.PostambleOffset_ms;
    resolved=isfinite(rows.AnchorTime_s);
    rows.ProjectedPrefix_ms(resolved)=max(0,rows.RFStartRelative_ms(resolved)-rows.RFStartRelative_ms(refIndex));
else
    warning('SCIONX:Postamble:NoReference','Reference flag group was not unambiguously identified. Set ReferencePostambleTime_s after reviewing candidates.');
end
writetable(rows,csv);
if ~isempty(allCandidates), writetable(allCandidates,candidateCsv); else, candidateCsv=""; end
result=struct('Configuration',options,'AlignmentResults',rows,'Candidates',allCandidates, ...
    'ReferenceRow',refIndex,'ComparisonRows',scanResult.ComparisonRows, ...
    'Signals',{signals},'InventoryCsv',csv,'CandidatesCsv',candidateCsv, ...
    'PlotFiles',strings(0,1),'MatFile',fullfile(options.OutputDirectory,'Postamble_Alignment.mat'));
if isfinite(rows.AnchorTime_s(refIndex))
    indices=scanResult.ComparisonRows;
    indices=indices(isfinite(rows.AnchorTime_s(indices)));
    if options.ExportPlots
        result.PlotFiles=exportPages(rows,signals,refIndex,indices,options);
    end
    if options.OpenViewer && ~isempty(indices)
        openViewer(rows,signals,refIndex,indices,options);
    end
end
save(result.MatFile,'result');
if any(~isfinite(rows.AnchorTime_s))
    warning('SCIONX:Postamble:Unresolved','Some bursts are unresolved; inspect CSV and candidate times. They are omitted from aligned plots.');
end
fprintf('Saved postamble alignment: %s\n',options.OutputDirectory);
end

function r=emptyRecord()
r=struct('BurstRow',0,'FileIndex',0,'BurstID',0,'FilePath',"",'FileLabel',"", ...
    'InputFormat',"",'SampleRate',0,'RFStart_s',0,'RFEnd_s',0,'DisplayPolarity',1, ...
    'Status',"not_found",'AnchorTime_s',NaN,'DetectedAnchorTime_s',NaN, ...
    'AnchorSampleApprox',NaN,'PostambleOffset_ms',NaN,'RFStartRelative_ms',NaN, ...
    'RFEndRelative_ms',NaN,'AdvanceRelativeToReference_ms',NaN,'ProjectedPrefix_ms',NaN, ...
    'FlagCount',0,'SymbolRate',NaN,'LineCode',"",'BitPolarity',NaN, ...
    'FitCorrelation',NaN,'CarrierOffsetHz',NaN,'FragmentCount',0, ...
    'FlagGroupDuration_ms',NaN,'FragmentAnchors_s',"",'MaxFragmentGap_ms',NaN,'SelectionReason',"");
end

function plotRow(layout,position,row,signal,ref,o,role,color)
t=(signal.Time_s-row.AnchorTime_s)*1000;
for column=1:2
    ax=nexttile(layout,position+column-1); hold(ax,'on');
    if column==1
        span=o.Window_ms; smoothing=o.Smoothing_ms;
    else
        span=o.DetailWindow_ms; smoothing=o.DetailSmoothing_ms;
    end
    smoothLength=max(1,round(smoothing/1000*signal.SampleRate));
    frequency=movmean(signal.DiscriminatorHz,smoothLength)/1000;
    frequency(signal.Time_s<row.RFStart_s+.001 | signal.Time_s>row.RFEnd_s-.001)=NaN;
    limits=[-1 1]*o.FrequencyLimitHz/1000;
    if o.ShadeProjectedPrefix && row.RFStartRelative_ms-ref.RFStartRelative_ms>5
        a=max(span(1),ref.RFStartRelative_ms); b=min(span(2),row.RFStartRelative_ms);
        if a<b
            patch(ax,[a b b a],[limits(1) limits(1) limits(2) limits(2)], ...
                [.9 .92 .95],'EdgeColor','none','FaceAlpha',.7);
            text(ax,(a+b)/2,limits(2)*.85,sprintf('%.1f ms\nreference RF extent',b-a), ...
                'HorizontalAlignment','center','FontSize',8,'Color',[.3 .35 .4]);
        end
    end
    duration=row.FlagGroupDuration_ms;
    if isfinite(duration)
        patch(ax,[0 duration duration 0],[limits(1) limits(1) limits(2) limits(2)], ...
            [1 .85 .5],'EdgeColor','none','FaceAlpha',.4);
    end
    plot(ax,t,frequency,'Color',color,'LineWidth',.8);
    xline(ax,0,'k-','first flag','LabelVerticalAlignment','bottom');
    grid(ax,'on'); xlim(ax,span); ylim(ax,limits);
    xlabel(ax,'Time relative to first flag (ms)'); ylabel(ax,'Frequency (kHz)');
    if column==1
        title(ax,sprintf('%s #%03d | anchor %.6f s | offset from RF %.2f ms | advance %+.2f ms\n%s', ...
            role,row.BurstID,row.AnchorTime_s,row.PostambleOffset_ms,row.AdvanceRelativeToReference_ms,row.FileLabel), ...
            'Interpreter','none','FontSize',9);
    else
        title(ax,sprintf('Detail | %d matched flags / %d runs | %s | fit %.3f\n%s | %.2f ms display smoothing', ...
            row.FlagCount,row.FragmentCount,row.LineCode,row.FitCorrelation,row.Status,smoothing), ...
            'Interpreter','none','FontSize',9);
    end
end
end

function paths=exportPages(rows,signals,refIndex,indices,o)
paths=strings(0,1);
% Still export the reference alone when there are no resolved comparisons.
pages=max(1,ceil(numel(indices)/o.BurstsPerPage));
colors=[.1 .35 .7;.75 .2 .3;0 .55 .45;.5 .3 .7];
for page=1:pages
    current=indices((page-1)*o.BurstsPerPage+1:min(page*o.BurstsPerPage,numel(indices)));
    fig=figure('Visible','off','Color','w','Position',[50 50 1500 290*(numel(current)+1)]);
    cleanup=onCleanup(@() close(fig));
    layout=tiledlayout(fig,numel(current)+1,2,'TileSpacing','compact','Padding','compact');
    title(layout,'SCION-X | aligned at first closing-flag / postamble candidate','FontSize',13);
    subtitle(layout,'Zero = flag-group start. Orange = detected group. Grey = projected July RF extent; not proven missing data. No CRC claim.');
    plotRow(layout,1,rows(refIndex,:),signals{refIndex},rows(refIndex,:),o,'JULY REFERENCE',colors(1,:));
    for j=1:numel(current)
        idx=current(j);
        plotRow(layout,2*j+1,rows(idx,:),signals{idx},rows(refIndex,:),o,'CANDIDATE',colors(mod(j-1,3)+2,:));
    end
    path=fullfile(o.OutputDirectory,sprintf('Postamble_Aligned_%03d.png',page));
    drawnow; set(fig,'PaperPositionMode','auto'); print(fig,path,'-dpng','-r120');
    paths(end+1,1)=path; %#ok<AGROW>
    clear cleanup
end
end

function openViewer(rows,signals,refIndex,indices,o)
fig=figure('Name','SCION-X Postamble Alignment','NumberTitle','off','Color','w', ...
    'Position',[80 80 1450 820]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
layout.OuterPosition=[0 .1 1 .9];
labels=strings(numel(indices),1);
for j=1:numel(indices)
    r=rows(indices(j),:);
    labels(j)=sprintf('%d/%d | #%03d | advance %+.2f ms | %s', ...
        j,numel(indices),r.BurstID,r.AdvanceRelativeToReference_ms,r.FileLabel);
end
selection=1;
menu=uicontrol(fig,'Style','popupmenu','Units','normalized','Position',[.19 .035 .63 .045], ...
    'String',cellstr(labels),'Callback',@select);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.03 .035 .13 .045], ...
    'String','Previous','Callback',@(~,~) move(-1));
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.85 .035 .12 .045], ...
    'String','Next','Callback',@(~,~) move(1));
draw();
    function select(src,~)
        selection=src.Value; draw();
    end
    function move(delta)
        selection=max(1,min(numel(indices),selection+delta)); menu.Value=selection; draw();
    end
    function draw()
        delete(layout.Children);
        plotRow(layout,1,rows(refIndex,:),signals{refIndex},rows(refIndex,:),o,'JULY REFERENCE',[.1 .35 .7]);
        idx=indices(selection);
        plotRow(layout,3,rows(idx,:),signals{idx},rows(refIndex,:),o,'CANDIDATE',[.75 .2 .3]);
        drawnow;
    end
end
