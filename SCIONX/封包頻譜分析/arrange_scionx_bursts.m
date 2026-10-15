function result = arrange_scionx_bursts(referenceFile, comparisonFiles, options)
%ARRANGE_SCIONX_BURSTS Scan CS16/CF32 IQ and align bursts for comparison.
% No Signal Processing Toolbox required. MATLAB R2022a or later.
% Views begin Pretrigger_ms before RF onset (not an HDLC frame boundary).
% InputFormat/SampleRate/DisplayPolarity: scalar or per-input column vector.
% See README_Burst_Comparison_zh-TW.md for configuration and interpretation.
arguments
    referenceFile (1,1) string
    comparisonFiles (:,1) string
    options.OutputDirectory (1,1) string = ""
    options.SampleRate (:,1) double {mustBePositive,mustBeFinite} = 100000
    options.InputFormat (:,1) string {mustBeMember(options.InputFormat,["cs16","cf32"])} = "cs16"
    options.DetectionBandHz (1,2) double {mustBeFinite} = [-12000 12000]
    options.DetectionFFT (1,1) double {mustBeInteger,mustBePositive} = 1024
    options.BatchFrames (1,1) double {mustBeInteger,mustBePositive} = 256
    options.NoiseWindow_s (1,1) double {mustBePositive,mustBeFinite} = 5
    options.Threshold_dB (1,1) double {mustBePositive,mustBeFinite} = 5
    options.MinimumBurst_ms (1,1) double {mustBePositive,mustBeFinite} = 80
    options.MaximumGap_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 25
    options.Pretrigger_ms (1,1) double {mustBeNonnegative,mustBeFinite} = 7
    options.ViewDuration_ms (1,1) double {mustBePositive,mustBeFinite} = 450
    options.DetailDuration_ms (1,1) double {mustBePositive,mustBeFinite} = 40
    options.FrequencySpanHz (1,1) double {mustBePositive,mustBeFinite} = 10000
    options.SpectrogramFFT (1,1) double {mustBeInteger,mustBePositive} = 256
    options.SpectrogramHop (1,1) double {mustBeInteger,mustBePositive} = 25
    options.ColorLimits_dB (1,2) double {mustBeFinite} = [-5 25]
    options.CorrectCarrier (1,1) logical = true
    options.DisplayPolarity (:,1) double {mustBeMember(options.DisplayPolarity,[-1 1])} = 1
    options.ReceiverLowpassHz (1,1) double {mustBePositive,mustBeFinite} = 5000
    options.DiscriminatorLimitHz (1,1) double {mustBePositive,mustBeFinite} = 5000
    options.ReferenceTime_s (1,1) double = NaN
    options.SortBy (1,1) string {mustBeMember(options.SortBy,["time","strength","duration"])} = "time"
    options.BurstsPerPage (1,1) double {mustBeInteger,mustBePositive} = 4
    options.ExportPlots (1,1) logical = true
    options.OpenViewer (1,1) logical = true
end
paths = [referenceFile; comparisonFiles];
options.SampleRate = expandSetting(options.SampleRate,numel(paths),'SampleRate');
options.InputFormat = expandSetting(options.InputFormat,numel(paths),'InputFormat');
options.DisplayPolarity = expandSetting(options.DisplayPolarity,numel(paths),'DisplayPolarity');
% All frequency settings must be valid at the lowest input sample rate.
fs = min(options.SampleRate);
assert(options.DetectionFFT >= 16 && options.SpectrogramFFT >= 16, ...
    'SCIONX:Burst:FFT','FFT length must be at least 16.');
assert(options.DetectionBandHz(1) < options.DetectionBandHz(2) && ...
    all(abs(options.DetectionBandHz) < fs/2), ...
    'SCIONX:Burst:Band','Detection band must lie inside Nyquist.');
assert(options.FrequencySpanHz < fs/2 && options.ColorLimits_dB(1) < options.ColorLimits_dB(2), ...
    'SCIONX:Burst:Display','Invalid frequency span or color limits.');
assert(options.ReceiverLowpassHz < fs/2 && options.DiscriminatorLimitHz < fs/2, ...
    'SCIONX:Burst:Receiver','Receiver lowpass and discriminator limits must be below Nyquist.');
assert(options.SpectrogramHop <= options.SpectrogramFFT, ...
    'SCIONX:Burst:Hop','SpectrogramHop must not exceed SpectrogramFFT.');
assert(options.ViewDuration_ms > options.Pretrigger_ms && ...
    options.DetailDuration_ms <= options.ViewDuration_ms, ...
    'SCIONX:Burst:Time','Invalid display duration.');
assert(isnan(options.ReferenceTime_s) || (isfinite(options.ReferenceTime_s) && options.ReferenceTime_s >= 0), ...
    'SCIONX:Burst:ReferenceTime','ReferenceTime_s must be NaN or a nonnegative finite time.');
if strlength(options.OutputDirectory) == 0
    options.OutputDirectory = fullfile(fileparts(mfilename('fullpath')), 'Burst_Comparison_Output');
end
% Validate every input before generating output: 4 bytes CS16, 8 bytes CF32.
for k = 1:numel(paths)
    bytes = formatSpec(options.InputFormat(k));
    info = dir(paths(k));
    assert(isfile(paths(k)) && isscalar(info) && info.bytes >= bytes && mod(info.bytes,bytes) == 0, ...
        'SCIONX:Burst:File','Missing, empty, or incomplete %s IQ file: %s',options.InputFormat(k),paths(k));
end
if ~isfolder(options.OutputDirectory), mkdir(options.OutputDirectory); end
records = repmat(emptyRecord(),0,1);
scans = cell(numel(paths),1);
for k = 1:numel(paths)
    fileOptions = settingsForFile(options,k);
    fprintf('[%d/%d] Scanning %s | %s, %.0f samples/s\n', ...
        k,numel(paths),paths(k),fileOptions.InputFormat,fileOptions.SampleRate);
    [r,scans{k}] = scanFile(paths(k),k,fileOptions);
    records = [records;r]; %#ok<AGROW>
    fprintf('  %d RF candidates, duration %.3f s\n',numel(r),scans{k}.Duration_s);
end
bursts = struct2table(records);
csv = fullfile(options.OutputDirectory,'Burst_Inventory.csv');
writetable(bursts,csv);
refIndices = find(bursts.FileIndex == 1);
comparisonIndices = find(bursts.FileIndex > 1);
switch options.SortBy
    case "strength", [~,order] = sort(bursts.PeakExcess_dB(comparisonIndices),'descend');
    case "duration", [~,order] = sort(bursts.Duration_s(comparisonIndices),'descend');
    otherwise, [~,order] = sortrows([bursts.FileIndex(comparisonIndices),bursts.Start_s(comparisonIndices)]);
end
comparisonIndices = comparisonIndices(order);
result = struct('Configuration',options,'BurstResults',bursts,'Scans',{scans}, ...
    'InventoryCsv',csv,'ReferenceRow',NaN,'ComparisonRows',comparisonIndices, ...
    'PlotFiles',strings(0,1),'MatFile',fullfile(options.OutputDirectory,'Burst_Results.mat'));
if isempty(refIndices)
    warning('SCIONX:Burst:NoReference','No reference burst found. Inventory and scan metrics saved; lower Threshold_dB or adjust DetectionBandHz.');
    save(result.MatFile,'result');
    return
end
if isnan(options.ReferenceTime_s)
    [~,pick] = max(bursts.PeakExcess_dB(refIndices));
else
    [~,pick] = min(abs(bursts.Start_s(refIndices)-options.ReferenceTime_s));
end
refIndex = refIndices(pick);
result.ReferenceRow = refIndex;
fprintf('Reference RF candidate #%d at %.6f s\n',bursts.BurstID(refIndex),bursts.Start_s(refIndex));
reference = loadView(bursts(refIndex,:),options);
if options.ExportPlots
    % Reference atlas includes ALL detected reference bursts in time order.
    result.PlotFiles = exportPages(bursts,refIndices,[],options,'Reference_Atlas');
    result.PlotFiles = [result.PlotFiles;exportPages(bursts,comparisonIndices,reference,options,'Comparison')];
    result.PlotFiles(end+1,1) = exportDetection(scans,bursts,options);
end
save(result.MatFile,'result');
if isempty(comparisonIndices) && ~isempty(comparisonFiles)
    warning('SCIONX:Burst:NoComparison','No comparison burst found. Inspect Detection_Overview.png and adjust detector settings.');
elseif options.OpenViewer && ~isempty(comparisonIndices)
    openViewer(bursts,comparisonIndices,reference,options);
end
fprintf('Saved: %s\n',options.OutputDirectory);
end

function r = emptyRecord()
r = struct('FileIndex',0,'FilePath',"",'FileLabel',"",'BurstID',0, ...
    'InputFormat',"",'SampleRate',0,'BytesPerComplexSample',0,'DisplayPolarity',1, ...
    'StartSample',0,'EndSampleExclusive',0,'Start_s',0,'End_s',0, ...
    'Duration_s',0,'PeakExcess_dB',0,'MedianExcess_dB',0, ...
    'NoiseBandPower_dBFS',0,'TouchesFileEdge',false,'ViewClipped',false);
end

function [records,scan] = scanFile(path,fileIndex,o)
[bytes,primitive,scale] = formatSpec(o.InputFormat);
info = dir(path); sampleCount = info.bytes/bytes; fs = o.SampleRate;
n = o.DetectionFFT; count = ceil(sampleCount/n);
f = (-floor(n/2):ceil(n/2)-1)'*fs/n;
band = f >= o.DetectionBandHz(1) & f <= o.DetectionBandHz(2);
assert(any(band),'SCIONX:Burst:BandBins','No FFT bins in detection band.');
w = .5-.5*cos(2*pi*(0:n-1)'/(n-1));
power = zeros(count,1);
fid = fopen(path,'rb','ieee-le'); assert(fid >= 0,'Cannot open %s',path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
for first = 1:o.BatchFrames:count
    number = min(o.BatchFrames,count-first+1);
    raw = fread(fid,[2,n*number],char(primitive+"=>single"));
    assert(all(isfinite(raw(:))),'SCIONX:Burst:Nonfinite','Nonfinite IQ samples in %s',path);
    % Nonconjugating transpose preserves I+jQ and the frequency direction.
    z = complex(raw(1,:).',raw(2,:).')/scale;
    actual = numel(z);
    z(end+1:n*number) = 0;
    z = reshape(z,n,number);
    ps = abs(fftshift(fft(z.*w,[],1),1)).^2/(n*sum(w.^2));
    p = sum(ps(band,:),1);
    % Last partial frame is normalized by its actual window energy.
    if actual < n*number
        valid = actual-(number-1)*n;
        energy = sum(w(1:valid).^2);
        if energy > 0, p(end) = p(end)*sum(w.^2)/energy; else, p(end) = eps; end
    end
    power(first:first+number-1) = p;
end
db = 10*log10(max(power,realmin));
window = max(3,2*floor(o.NoiseWindow_s*fs/n/2)+1);
noise = movmedian(db,window,'Endpoints','shrink');
excess = db-noise;
mask = excess >= o.Threshold_dB;
[starts,ends] = runs(mask);
% Merge short internal dropouts but retain original candidate start/end.
if ~isempty(starts)
    mergedStart = starts(1); mergedEnd = ends(1);
    for j = 2:numel(starts)
        if (starts(j)-mergedEnd(end)-1)*n/fs <= o.MaximumGap_ms/1000
            mergedEnd(end) = ends(j);
        else
            mergedStart(end+1,1) = starts(j); %#ok<AGROW>
            mergedEnd(end+1,1) = ends(j); %#ok<AGROW>
        end
    end
    starts = mergedStart; ends = mergedEnd;
end
records = repmat(emptyRecord(),0,1);
[~,name,ext] = fileparts(path);
for j = 1:numel(starts)
    a = (starts(j)-1)*n; b = min(ends(j)*n,sampleCount);
    if (b-a)/fs < o.MinimumBurst_ms/1000, continue; end
    % Refine RF onset/offset with 1-ms band-limited block energies.
    lo = max(0,a-n); hi = min(sampleCount,b+n);
    z = readIQ(path,lo,hi-lo,o.InputFormat);
    mid = mean(o.DetectionBandHz);
    z = z.*exp(-2i*pi*mid*(0:numel(z)-1)'/fs);
    taps = lowpassTaps(diff(o.DetectionBandHz)/2,fs);
    z = conv(z,taps,'same');
    step = max(1,round(fs*.001)); blocks = floor(numel(z)/step);
    p = mean(abs(reshape(z(1:blocks*step),step,blocks)).^2,1);
    baseline = median(noise(starts(j):ends(j)));
    hit = find(10*log10(max(p,realmin)) >= baseline+o.Threshold_dB);
    if ~isempty(hit)
        a = lo+(hit(1)-1)*step; b = min(hi,lo+hit(end)*step);
    end
    if (b-a)/fs < o.MinimumBurst_ms/1000, continue; end
    r = emptyRecord(); r.FileIndex = fileIndex; r.FilePath = path;
    r.InputFormat = o.InputFormat; r.SampleRate = fs;
    r.BytesPerComplexSample = bytes; r.DisplayPolarity = o.DisplayPolarity;
    r.FileLabel = string(name)+string(ext); r.BurstID = numel(records)+1;
    r.StartSample = a; r.EndSampleExclusive = b;
    r.Start_s = a/fs; r.End_s = b/fs; r.Duration_s = (b-a)/fs;
    r.PeakExcess_dB = max(excess(starts(j):ends(j)));
    r.MedianExcess_dB = median(excess(starts(j):ends(j)));
    r.NoiseBandPower_dBFS = baseline;
    r.TouchesFileEdge = a <= step || b >= sampleCount-step;
    origin = max(0,a-round(o.Pretrigger_ms/1000*fs));
    r.ViewClipped = b > origin+round(o.ViewDuration_ms/1000*fs);
    records(end+1,1) = r; %#ok<AGROW>
end
scan = struct('FilePath',path,'InputFormat',o.InputFormat,'SampleRate',fs, ...
    'DetectionFFT',o.DetectionFFT,'SpectrogramFFT',o.SpectrogramFFT, ...
    'SpectrogramHop',o.SpectrogramHop,'Duration_s',sampleCount/fs, ...
    'Time_s',((0:count-1)'*n+n/2)/fs,'Power_dBFS',db, ...
    'Noise_dBFS',noise,'Excess_dB',excess);
end

function [a,b] = runs(mask)
d = diff([false;mask(:);false]); a = find(d==1); b = find(d==-1)-1;
end

function value = expandSetting(value,count,name)
assert(isscalar(value) || numel(value)==count,'SCIONX:Burst:InputSettings', ...
    '%s must be scalar or contain one value per input (reference first, then comparisons).',name);
if isscalar(value), value = repmat(value,count,1); end
end

function o = settingsForFile(options,index)
o = options;
o.SampleRate = options.SampleRate(index);
o.InputFormat = options.InputFormat(index);
o.DisplayPolarity = options.DisplayPolarity(index);
% Match FFT time support/bin width without resampling original IQ.
ratio = o.SampleRate/options.SampleRate(1);
o.DetectionFFT = max(16,round(options.DetectionFFT*ratio));
o.SpectrogramFFT = max(16,round(options.SpectrogramFFT*ratio));
o.SpectrogramHop = min(o.SpectrogramFFT,max(1,round(options.SpectrogramHop*ratio)));
end

function [bytes,primitive,scale] = formatSpec(format)
switch format
    case "cs16", bytes = 4; primitive = "int16"; scale = 32768;
    case "cf32", bytes = 8; primitive = "single"; scale = 1;
    otherwise, error('SCIONX:Burst:Format','Unsupported IQ format: %s',format);
end
end

function z = readIQ(path,start,count,format)
if count<=0, z = complex(zeros(0,1)); return; end
[bytes,primitive,scale] = formatSpec(format);
fid = fopen(path,'rb','ieee-le'); assert(fid>=0,'Cannot open %s',path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
assert(fseek(fid,bytes*start,'bof')==0,'Cannot seek in IQ file.');
r = fread(fid,[2,count],char(primitive+"=>double"));
assert(all(isfinite(r(:))),'SCIONX:Burst:Nonfinite','Nonfinite IQ samples in %s',path);
z = complex(r(1,:).',r(2,:).')/scale;
end

function h = lowpassTaps(cutoff,fs)
% Windowed sinc; odd length => conv(...,'same') preserves sample alignment.
m = (-64:64)'; q = 2*cutoff/fs*m;
s = ones(size(q)); hit = q~=0; s(hit) = sin(pi*q(hit))./(pi*q(hit));
h = 2*cutoff/fs*s.*(.54+.46*cos(pi*m/64)); h = h/sum(h);
end

function v = loadView(row,o)
o = settingsForFile(o,row.FileIndex);
fs = o.SampleRate; info = dir(row.FilePath); total = info.bytes/row.BytesPerComplexSample;
origin = max(0,row.StartSample-round(o.Pretrigger_ms/1000*fs));
count = min(total-origin,round(o.ViewDuration_ms/1000*fs));
z = readIQ(row.FilePath,origin,count,row.InputFormat);
if o.DisplayPolarity == -1, z = conj(z); end
% Estimate constant carrier midpoint from trimmed instantaneous frequencies.
% This is a visual aid, not Doppler tracking, symbol recovery, or CRC validation.
filtered = conv(z,lowpassTaps(min(fs*.4,max(abs(o.DetectionBandHz))),fs),'same');
fm = angle(filtered(2:end).*conj(filtered(1:end-1)))*fs/(2*pi);
time = (origin+(1:numel(fm))')/fs;
active = fm(time>=row.Start_s+.005 & time<=row.End_s-.005);
active = sort(active(isfinite(active)));
offset = 0;
if o.CorrectCarrier && ~isempty(active)
    offset = mean(active(max(1,round([.1 .9]*numel(active)))));
end
z = z.*exp(-2i*pi*offset*(0:numel(z)-1)'/fs);
filtered = conv(z,lowpassTaps(o.ReceiverLowpassHz,fs),'same');
fm = angle(filtered(2:end).*conj(filtered(1:end-1)))*fs/(2*pi)/1000;
fm = conv(fm,[.25;.5;.25],'same');
[ps,t,f] = stftPower(z,o);
% Nearby OFF-burst windows on either side, separated by 20 ms guards.
guard = round(.020*fs); lengthNoise = round(.100*fs);
beforeEnd = max(0,row.StartSample-guard);
afterStart = min(total,row.EndSampleExclusive+guard);
before = readIQ(row.FilePath,max(0,beforeEnd-lengthNoise),min(beforeEnd,lengthNoise),row.InputFormat);
after = readIQ(row.FilePath,afterStart,min(lengthNoise,total-afterStart),row.InputFormat);
noiseFrames = [];
for part = {before,after}
    nz = part{1};
    if numel(nz)>=o.SpectrogramFFT
        if o.DisplayPolarity == -1, nz = conj(nz); end
        nz = nz.*exp(-2i*pi*offset*(0:numel(nz)-1)'/fs);
        np = stftPower(nz,o); noiseFrames = [noiseFrames,np]; %#ok<AGROW>
    end
end
if isempty(noiseFrames)
    noisePSD = 10^(row.NoiseBandPower_dBFS/10)/diff(o.DetectionBandHz);
else
    noisePSD = median(noiseFrames(:));
end
relative = 10*log10(max(ps,realmin)/max(noisePSD,realmin));
onsetMs = (row.StartSample-origin)/fs*1000;
v = struct('Row',row,'Origin_s',origin/fs,'OffsetHz',offset, ...
    'Time_ms',t*1000,'Frequency_kHz',f/1000,'Spectrogram_dB',relative, ...
    'DiscriminatorTime_ms',(1:numel(fm))'/fs*1000, ...
    'Discriminator_kHz',fm,'Onset_ms',onsetMs, ...
    'End_ms',(row.EndSampleExclusive-origin)/fs*1000);
end

function [p,t,f] = stftPower(z,o)
n = o.SpectrogramFFT; hop = o.SpectrogramHop;
if numel(z)<n, z(end+1:n)=0; end
starts = 1:hop:numel(z)-n+1;
indices = (0:n-1)'+starts;
w = .5-.5*cos(2*pi*(0:n-1)'/(n-1));
p = abs(fftshift(fft(z(indices).*w,[],1),1)).^2/(o.SampleRate*sum(w.^2));
f = (-floor(n/2):ceil(n/2)-1)'*o.SampleRate/n;
keep = abs(f)<=o.FrequencySpanHz; f = f(keep); p = p(keep,:);
t = (starts-1+(n-1)/2)/o.SampleRate;
end

function plotRow(layout,position,v,o,role)
r = v.Row;
ax = nexttile(layout,position);
imagesc(ax,v.Time_ms,v.Frequency_kHz,v.Spectrogram_dB); axis(ax,'xy');
colormap(ax,turbo); clim(ax,o.ColorLimits_dB);
xlim(ax,[0 o.ViewDuration_ms]); ylim(ax,[-1 1]*o.FrequencySpanHz/1000);
xline(ax,v.Onset_ms,'w:','RF onset'); xline(ax,v.End_ms,'w:');
ylabel(ax,'Frequency (kHz)'); xlabel(ax,'Time from view origin (ms)');
flag = "";
if r.ViewClipped, flag = " | TRUNCATED VIEW"; end
if r.TouchesFileEdge, flag = flag+" | FILE EDGE"; end
title(ax,sprintf('%s #%03d | RF %.6f s | %.1f ms | peak +%.1f dB%s\n%s | %s, %.0f kS/s', ...
    role,r.BurstID,r.Start_s,r.Duration_s*1000,r.PeakExcess_dB,flag,r.FileLabel,r.InputFormat,r.SampleRate/1000), ...
    'Interpreter','none','FontSize',9);
cb = colorbar(ax); cb.Label.String = 'dB / nearby off-burst PSD';
ax = nexttile(layout,position+1);
plot(ax,v.DiscriminatorTime_ms,v.Discriminator_kHz,'Color',[.03 .3 .65]);
xlim(ax,[0 o.DetailDuration_ms]); ylim(ax,[-1 1]*o.DiscriminatorLimitHz/1000);
grid(ax,'on'); xline(ax,v.Onset_ms,'k:');
xlabel(ax,'Same view origin (ms)'); ylabel(ax,'Discriminator (kHz)');
title(ax,sprintf('First %.0f ms | origin %.6f s\ncarrier correction %+.0f Hz | display polarity %+.0f', ...
    o.DetailDuration_ms,v.Origin_s,v.OffsetHz,r.DisplayPolarity),'FontSize',9);
end

function paths = exportPages(bursts,indices,reference,o,prefix)
paths = strings(0,1);
for first = 1:o.BurstsPerPage:numel(indices)
    current = indices(first:min(first+o.BurstsPerPage-1,numel(indices)));
    rows = numel(current)+~isempty(reference);
    fig = figure('Visible','off','Color','w','Position',[50 50 1500 290*rows]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig,rows,2,'TileSpacing','compact','Padding','compact');
    title(layout,'SCION-X RF burst comparison | common time/frequency/color scales','FontSize',13);
    subtitle(layout,'RF candidates only; reference selection does not validate CRC. Dotted lines: RF onset/end.');
    position = 1;
    if ~isempty(reference)
        plotRow(layout,position,reference,o,'REFERENCE'); position = position+2;
    end
    for idx = current(:)'
        v = loadView(bursts(idx,:),o);
        plotRow(layout,position,v,o,'BURST'); position = position+2;
    end
    path = fullfile(o.OutputDirectory,sprintf('%s_%03d.png',prefix,ceil(first/o.BurstsPerPage)));
    drawnow;
    set(fig,'PaperPositionMode','auto');
    print(fig,path,'-dpng','-r120'); paths(end+1,1) = path; %#ok<AGROW>
    clear cleanup
end
end

function path = exportDetection(scans,bursts,o)
fig = figure('Visible','off','Color','w','Position',[50 50 1400 260*numel(scans)]);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig,numel(scans),1,'TileSpacing','compact');
for j = 1:numel(scans)
    ax = nexttile(layout); s = scans{j};
    plot(ax,s.Time_s,s.Excess_dB); hold(ax,'on'); yline(ax,o.Threshold_dB,'r--','threshold');
    rows = find(bursts.FileIndex==j);
    for idx = rows(:)'
        xline(ax,bursts.Start_s(idx),'k:',sprintf('#%d',bursts.BurstID(idx)));
    end
    grid(ax,'on'); xlim(ax,[0 s.Duration_s]);
    xlabel(ax,'Recording time (s)'); ylabel(ax,'Band power / local median (dB)');
    title(ax,s.FilePath,'Interpreter','none','FontSize',9);
end
path = fullfile(o.OutputDirectory,'Detection_Overview.png');
drawnow;
set(fig,'PaperPositionMode','auto');
print(fig,path,'-dpng','-r120');
end

function openViewer(bursts,indices,reference,o)
fig = figure('Name','SCION-X Burst Browser','NumberTitle','off', ...
    'Color','w','Position',[80 80 1450 820]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
layout.OuterPosition = [0 .1 1 .9];
labels = strings(numel(indices),1);
for j = 1:numel(indices)
    r = bursts(indices(j),:);
    labels(j) = sprintf('%d/%d | #%03d | %.6f s | +%.1f dB | %s', ...
        j,numel(indices),r.BurstID,r.Start_s,r.PeakExcess_dB,r.FileLabel);
end
selection = 1;
menu = uicontrol(fig,'Style','popupmenu','Units','normalized', ...
    'Position',[.19 .035 .63 .045],'String',cellstr(labels),'Callback',@select);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.03 .035 .13 .045], ...
    'String','Previous','Callback',@(~,~) move(-1));
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.85 .035 .12 .045], ...
    'String','Next','Callback',@(~,~) move(1));
draw();
    function select(src,~)
        selection = src.Value; draw();
    end
    function move(delta)
        selection = max(1,min(numel(indices),selection+delta));
        menu.Value = selection; draw();
    end
    function draw()
        delete(layout.Children);
        plotRow(layout,1,reference,o,'REFERENCE');
        plotRow(layout,3,loadView(bursts(indices(selection),:),o),o,'COMPARISON');
        drawnow;
    end
end
