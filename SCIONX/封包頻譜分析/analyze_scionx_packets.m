function result=analyze_scionx_packets(comparisonFiles,options)
% Integrated postamble-anchored HDLC/TLM with same-burst bit diagnostics.
% Reference is used internally; only comparison inputs produce TLM content.
arguments
    comparisonFiles (:,1) string
    options.ReferenceFile (1,1) string = ""
    options.ReferenceSampleRate (1,1) double = 100000
    options.ReferenceTime_s (1,1) double = 306.6
    options.ReferencePostambleTime_s (1,1) double = NaN
    options.ResolveRepeatedPostambles (1,1) logical = false
    options.InputFormat (:,1) string = "auto"
    options.SampleRate (:,1) double = NaN
    options.OutputDirectory (1,1) string = ""
    options.AlignmentMatFile (1,1) string = ""
    options.MinimumAgreement_pct (1,1) double {mustBeNonnegative,mustBeFinite} = 80
    options.MinimumComparedBits (1,1) double {mustBePositive,mustBeInteger} = 128
    options.MinimumCoverage_pct (1,1) double {mustBePositive,mustBeFinite} = 80
    options.DecodePolicy (1,1) string {mustBeMember(options.DecodePolicy,["postamble","bit_gate"])} = "postamble"
    options.AfterWindow_ms (1,2) double = [8.5 49.5]
    options.DetectionThreshold_dB (1,1) double = 5
    options.AudioFlagCorrelation (1,1) double {mustBePositive,mustBeFinite} = .70
    options.ExpectedFrameBytes (1,1) double {mustBeInteger,mustBePositive} = 274
    options.ExportPlots (1,1) logical = true
    options.TLMDecoderPath (1,1) string = ""
    options.AudioLevelRecoveryModels (:,2) double = zeros(0,2)
end
base=string(fileparts(mfilename('fullpath')));
assert(~isempty(comparisonFiles),'Provide at least one comparison recording.');
assert(options.MinimumAgreement_pct<=100 && options.MinimumCoverage_pct<=100);
assert(options.AudioFlagCorrelation<=1);
if options.ReferenceFile==""
    options.ReferenceFile=fullfile(base,'..','Beacon Data','CS16_Raw','20260723', ...
        '20260723_091639_100000_98266_x.cs16');
end
if options.TLMDecoderPath==""
    options.TLMDecoderPath=fullfile(base,'..','SignalSim','GNU_Matlab_Demo','decode_scionx_tlm.m');
end
assert(isfile(options.TLMDecoderPath),'TLM decoder missing: %s',options.TLMDecoderPath);
files=[options.ReferenceFile;comparisonFiles];
for k=1:numel(files)
    assert(isfile(files(k)),'Recording missing: %s',files(k));
    files(k)=string(java.io.File(char(files(k))).getCanonicalPath());
end
assert(numel(unique(lower(files)))==numel(files),'Reference and comparison inputs must be distinct.');
% Name outputs after the original recording, including converted IQ inputs.
[~,sourceNames,sourceExtensions]=fileparts(files);
sourceNames=sourceNames+sourceExtensions;
% If filenames collide in different directories, keep every input unique.
if numel(unique(lower(sourceNames)))<numel(sourceNames)
    sourceNames=sourceNames+compose('_File_%03d',(1:numel(files))');
end
formats=expand(options.InputFormat,numel(comparisonFiles));
rates=expand(options.SampleRate,numel(comparisonFiles));
for k=1:numel(comparisonFiles)
    [formats(k),rates(k)]=inputSpec(files(k+1),formats(k),rates(k));
end
formats=["cs16";formats]; rates=[options.ReferenceSampleRate;rates];
root=options.OutputDirectory;
if root=="", root=fullfile(base,'Integrated_Packet_Output'); end
root=string(java.io.File(char(root)).getCanonicalPath());
root=fullfile(root,"Run_"+string(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
mkdir(root);
manifest=table((1:numel(files))',files,formats,rates,'VariableNames', ...
    {'FileIndex','OriginalPath','InputFormat','SampleRate'});
writetable(manifest,fullfile(root,'Input_Manifest.csv'));
if options.AlignmentMatFile~=""
    loaded=load(options.AlignmentMatFile,'result'); alignment=loaded.result;
    % Cache is explicit, never silently reuse an unrelated recording set.
    for k=1:numel(files)
        paths=unique(alignment.AlignmentResults.FilePath(alignment.AlignmentResults.FileIndex==k));
        assert(numel(paths)==1 && strcmpi(char(java.io.File(char(paths)).getCanonicalPath()),files(k)), ...
            'Alignment cache does not match input %d.',k);
    end
    assert(max(alignment.AlignmentResults.FileIndex)==numel(files),'Cache has extra inputs.');
else
    scan=scan_scionx_pipeline_inputs(files,formats,rates,fullfile(root,'Scan'), ...
        Threshold_dB=options.DetectionThreshold_dB,ReferenceTime_s=options.ReferenceTime_s, ...
        AudioFlagCorrelation=options.AudioFlagCorrelation);
    alignment=align_scionx_postambles(scan, ...
        ReferencePostambleTime_s=options.ReferencePostambleTime_s, ...
        ResolveRepeatedPostambles=options.ResolveRepeatedPostambles, ...
        OutputDirectory=fullfile(root,'Alignment'),ExportPlots=options.ExportPlots,OpenViewer=false);
end
cutoffs=repmat(5000,numel(files),1); cutoffs(formats~="cs16")=8000;
repetition=compare_scionx_repeated_segments(alignment, ...
    ReceiverLowpassHz=cutoffs,AfterWindow_ms=options.AfterWindow_ms,IncludeReference=false, ...
    OutputDirectory=fullfile(root,'Repetition'),ExportPlots=options.ExportPlots,OpenViewer=false);
% Internal positive control: no reference hex/bitstream/TLM is published.
ref=alignment.AlignmentResults(alignment.ReferenceRow,:);
control=struct('CRC',"NOT_CHECKED",'FrameBytes',0,'Status',"unresolved_anchor");
if isfinite(ref.AnchorTime_s)
    refSignal=readSignal(ref,cutoffs(ref.FileIndex));
    decoded=decode_scionx_postamble_frame(ref,refSignal,ExpectedFrameBytes=options.ExpectedFrameBytes);
    if decoded.Found
        control.CRC=decoded.Frame.CRC; control.FrameBytes=numel(decoded.Frame.Bytes);
        control.Status="checked";
    end
end
fprintf('Internal reference control: CRC %s, %d bytes (no reference TLM exported).\n',control.CRC,control.FrameBytes);
records=repmat(emptyRecord(),height(repetition.Summary),1);
bitSummary=table(); bitDetails=table();
for k=1:height(repetition.Summary)
    r=emptyRecord(); info=repetition.Summary(k,:);
    r.BurstRow=info.BurstRow; r.FileIndex=info.FileIndex; r.BurstID=info.BurstID;
    r.FileLabel=info.FileLabel; r.AnchorTime_s=info.AnchorTime_s;
    r.WaveformCorrelation=info.GlobalCorrelation;
    row=alignment.AlignmentResults(info.BurstRow,:);
    folder=fullfile(root,sourceNames(r.FileIndex)+sprintf('_Burst_%03d',r.BurstID)); mkdir(folder);
    try
        r.RepetitionStatus=info.Status;
        if ~isfinite(row.AnchorTime_s)
            r.Status="SKIPPED_UNRESOLVED_ANCHOR"; records(k)=r; continue
        end
        signal=readSignal(row,cutoffs(row.FileIndex));
        if isfinite(info.GlobalLag_ms)
        t=(signal.Time_s-row.AnchorTime_s)*1000;
        pieces=repetition.Segments(repetition.Segments.BurstRow==r.BurstRow,:);
        bounds=[info.AfterStart_ms info.AfterEnd_ms info.GlobalLag_ms; ...
            pieces.AfterStart_ms pieces.AfterEnd_ms pieces.Lag_ms];
        counts=zeros(size(bounds,1),3);
        for j=1:size(bounds,1)
            [s,b]=scionx_compare_bit_windows(t,signal.DiscriminatorHz,bounds(j,1:2),bounds(j,3),9600,1);
            label="whole_single_lag"; if j>1, label="stable_segment_"+(j-1); end
            s=addvars(s,r.BurstRow,label,'Before',1,'NewVariableNames',{'BurstRow','Window'});
            b=addvars(b,repmat(r.BurstRow,height(b),1),repmat(label,height(b),1), ...
                'Before',1,'NewVariableNames',{'BurstRow','Window'});
            bitSummary=[bitSummary;s]; bitDetails=[bitDetails;b]; %#ok<AGROW>
            counts(j,:)=[s.ComparedBits s.DifferentBits s.Agreement_pct];
        end
        r.GlobalAgreement_pct=counts(1,3);
        r.GateMode="whole_single_lag"; chosen=counts(1,1:2); coverage=100;
        if height(pieces)>1 && info.SegmentCoverage_pct>=options.MinimumCoverage_pct
            chosen=sum(counts(2:end,1:2),1);
            coverage=100*chosen(1)/counts(1,1); r.GateMode="stable_segments";
        end
        r.ComparedBits=chosen(1); r.DifferentBits=chosen(2);
        r.BitAgreement_pct=100*(1-chosen(2)/chosen(1)); r.Coverage_pct=coverage;
        r.BitGatePassed=r.BitAgreement_pct>=options.MinimumAgreement_pct && ...
            chosen(1)>=options.MinimumComparedBits && coverage>=options.MinimumCoverage_pct;
        end
        r.Eligible=r.BitGatePassed || options.DecodePolicy=="postamble";
        if ~r.Eligible
            r.Status="SKIPPED_BIT_GATE";
            if ~isfinite(info.GlobalLag_ms), r.Status="SKIPPED_"+upper(info.Status); end
            records(k)=r; continue
        end
        r.DecodeRoute="bit_gate";
        if options.DecodePolicy=="postamble", r.DecodeRoute="postamble_anchor"; end
        if row.InputFormat=="audio" && ~isempty(options.AudioLevelRecoveryModels)
            d=decode_scionx_postamble_frame(row,signal,ExpectedFrameBytes=options.ExpectedFrameBytes, ...
                LevelRecoveryModels=options.AudioLevelRecoveryModels);
            write_scionx_level_trials(d,folder,options.TLMDecoderPath,options.ExpectedFrameBytes);
            if d.RecoveryApplied
                r.DecodeRoute="postamble_level_restored";
            end
        else
            d=decode_scionx_postamble_frame(row,signal,ExpectedFrameBytes=options.ExpectedFrameBytes);
        end
        writetable(d.Candidates,fullfile(folder,'Timing_Candidates.csv'));
        if ~d.Found, r.Status="NO_HDLC_CANDIDATE"; records(k)=r; continue; end
        r.FrameBytes=numel(d.Frame.Bytes); r.CRC=d.Frame.CRC; r.Complete=d.Frame.Complete;
        r.DecodeSymbolRate=d.SymbolRate; r.StopReason=d.Frame.StopReason;
        save(fullfile(folder,'HDLC_Diagnostic.mat'),'d');
        serial=table(d.Edges_s(1:end-1),double(d.Bits),'VariableNames',{'BitStart_s','Bit'});
        writetable(serial,fullfile(folder,'Received_Raw_Bits.csv'));
        raw=fullfile(folder,'Received_Candidate.hex'); writeHex(raw,d.Frame.Bytes);
        if d.Frame.Complete && r.FrameBytes==options.ExpectedFrameBytes
            rawPath=fullfile(folder,'Raw_bit.txt'); writeHex(rawPath,d.Frame.Bytes);
            statusFile=fullfile(folder,'Matlab_HDLC_Status.txt');
            writeText(statusFile,"Boundary Mode: POSTAMBLE_REVERSE_STRICT");
            decoderFolder=fileparts(options.TLMDecoderPath); previousPath=path;
            cleanup=onCleanup(@() path(previousPath)); addpath(decoderFolder,'-begin');
            clear decode_scionx_tlm
            decode_scionx_tlm(rawPath); clear cleanup
            r.TLMFile=fullfile(folder,'TLM.txt'); r.Status="COMPLETE_"+r.CRC;
        else
            r.AssumedFirstByte=options.ExpectedFrameBytes-r.FrameBytes+1;
            r.TLMFile=fullfile(folder,'TLM_DIAGNOSTIC.txt');
            decode_scionx_partial_tlm(raw,r.AssumedFirstByte,r.TLMFile, ...
                options.TLMDecoderPath,options.ExpectedFrameBytes);
            r.Status="PARTIAL_DIAGNOSTIC";
            % A complete shorter frame does not justify the assumed map.
            if d.Frame.Complete, r.Status="SHORT_FRAME_DIAGNOSTIC"; end
        end
        metadata=sprintf(['Decode route = %s\nRepetition status = %s\n' ...
            'Bit gate passed = %d\nBit Agreement = %.4f%% (%d / %d equal), window %s\n' ...
            'Compared-window coverage = %.3f%%\nAnchor = %.9f s\n' ...
            'NaN agreement means unavailable, not zero or 100%%.\n' ...
            'Anchor and bit agreement do not validate CRC or telemetry.\n\n'], ...
            r.DecodeRoute,r.RepetitionStatus,r.BitGatePassed, ...
            r.BitAgreement_pct,r.ComparedBits-r.DifferentBits,r.ComparedBits,r.GateMode, ...
            r.Coverage_pct,r.AnchorTime_s);
        writeText(r.TLMFile,string(metadata)+string(fileread(r.TLMFile)));
    catch err
        r.Status="ERROR"; r.Error=string(err.identifier)+": "+string(err.message);
        writeText(fullfile(folder,'Error.txt'),r.Error);
    end
    records(k)=r;
    fprintf('Comparison file %d burst %d: %.2f%% bits, %s, CRC %s\n', ...
        r.FileIndex,r.BurstID,r.BitAgreement_pct,r.Status,r.CRC);
end
summary=struct2table(emptyRecord()); summary=summary([],:);
if ~isempty(records), summary=struct2table(records); end
writetable(summary,fullfile(root,'Packet_Summary.csv'));
recordings=manifest(2:end,:);
recordings.Bursts=zeros(height(recordings),1); recordings.Eligible=zeros(height(recordings),1);
recordings.TLMFiles=zeros(height(recordings),1); recordings.Status=strings(height(recordings),1);
for k=1:height(recordings)
    hit=summary.FileIndex==recordings.FileIndex(k);
    recordings.Bursts(k)=sum(hit); recordings.Eligible(k)=sum(summary.Eligible(hit));
    recordings.TLMFiles(k)=sum(summary.TLMFile(hit)~="");
    if ~any(hit), recordings.Status(k)="NO_BURSTS_FOUND";
    elseif ~any(summary.Eligible(hit)), recordings.Status(k)="NO_ELIGIBLE_BURSTS";
    elseif recordings.TLMFiles(k)==0, recordings.Status(k)="NO_TLM_CANDIDATE";
    else, recordings.Status(k)="TLM_WRITTEN"; end
end
writetable(recordings,fullfile(root,'Recording_Summary.csv'));
if ~isempty(bitSummary), writetable(bitSummary,fullfile(root,'Bit_Summary.csv')); end
if ~isempty(bitDetails), writetable(bitDetails,fullfile(root,'Bit_Details.csv')); end
result=struct('Configuration',options,'OutputDirectory',root,'Inputs',manifest, ...
    'Control',control,'Summary',summary,'Recordings',recordings,'BitSummary',bitSummary, ...
    'SummaryCsv',fullfile(root,'Packet_Summary.csv'),'RepetitionMat',repetition.MatFile);
save(fullfile(root,'Integrated_Result.mat'),'result');
writeText(fullfile(root,'Run_Status.txt'), ...
    "Reference CRC: "+control.CRC+newline+"Reference TLM: not exported"+newline+ ...
    "Comparison bursts: "+height(summary)+newline+"TLM files: "+sum(summary.TLMFile~=""));
disp(summary(:,{'FileLabel','BurstID','BitAgreement_pct','Eligible','Status','FrameBytes','CRC'}));
fprintf('Integrated output: %s\n',root);
end

function r=emptyRecord()
r=struct('BurstRow',0,'FileIndex',0,'BurstID',0,'FileLabel',"",'AnchorTime_s',NaN, ...
    'WaveformCorrelation',NaN,'GlobalAgreement_pct',NaN,'GateMode',"", ...
    'ComparedBits',0,'DifferentBits',0,'BitAgreement_pct',NaN,'Coverage_pct',0, ...
    'BitGatePassed',false,'DecodeRoute',"",'RepetitionStatus',"", ...
    'Eligible',false,'Status',"NOT_PROCESSED",'FrameBytes',0,'Complete',false, ...
    'CRC',"NOT_CHECKED",'DecodeSymbolRate',NaN,'StopReason',"", ...
    'AssumedFirstByte',NaN,'TLMFile',"",'Error',"");
end

function signal=readSignal(row,cutoff)
row.StartSample=round(row.RFStart_s*row.SampleRate); row.EndSampleExclusive=round(row.RFEnd_s*row.SampleRate);
row.Start_s=row.RFStart_s; row.End_s=row.RFEnd_s;
[~,signal]=detect_scionx_postamble(row,SignalOnly=true,ReceiverLowpassHz=cutoff,Padding_ms=280);
end

function value=expand(value,n)
assert(isscalar(value) || numel(value)==n,'Input setting must be scalar or one per comparison file.');
if isscalar(value), value=repmat(value,n,1); end
end

function [format,fs]=inputSpec(file,format,fs)
[~,name,ext]=fileparts(file); ext=lower(ext);
if format=="auto"
    switch ext
        case ".cs16", format="cs16";
        case {".iq",".cf32"}, format="cf32";
        case ".ogg"
            info=audioinfo(file);
            assert(info.NumChannels==1,'Stereo OGG is ambiguous; explicitly set InputFormat="ogg_iq" for I/Q.');
            format="audio";
        otherwise, error('Unknown input extension: %s',ext);
    end
end
assert(ismember(format,["cs16","cf32","audio","ogg_iq"]),'Unsupported format.');
if ismember(format,["audio","ogg_iq"])
    info=audioinfo(file); fs=info.SampleRate;
elseif isnan(fs)
    token=regexp(name,'_(\d+)SPS_','tokens','once');
    if isempty(token), token=regexp(name,'^\d{8}_\d{6}_(\d+)_','tokens','once'); end
    assert(~isempty(token),'Provide SampleRate for headerless raw input: %s',file);
    fs=str2double(token{1});
end
assert(isfinite(fs) && fs>24000,'Sample rate must exceed 24 kS/s.');
end

function writeHex(file,bytes)
fid=fopen(file,'wt'); assert(fid>=0); c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%02X ',bytes); fprintf(fid,'\n');
end

function writeText(file,content)
fid=fopen(file,'wt','n','UTF-8'); assert(fid>=0); c=onCleanup(@() fclose(fid)); %#ok<NASGU>
content=replace(string(content),sprintf('\r\n'),newline);
fprintf(fid,'%s\n',content);
end
