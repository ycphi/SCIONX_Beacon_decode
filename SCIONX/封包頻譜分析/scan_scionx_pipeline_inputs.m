function scan=scan_scionx_pipeline_inputs(files,formats,rates,outputDirectory,options)
% First input is the reference. Mono OGG is demodulated baseband, not I/Q.
arguments
    files (:,1) string
    formats (:,1) string
    rates (:,1) double
    outputDirectory (1,1) string
    options.Threshold_dB (1,1) double = 5
    options.ReferenceTime_s (1,1) double = 306.6
    options.AudioFlagCorrelation (1,1) double = .70
end
assert(numel(files)==numel(formats) && numel(files)==numel(rates));
assert(formats(1)=="cs16" || formats(1)=="cf32",'Reference must be raw IQ.');
rows=table(); scans=cell(numel(files),1);
for k=1:numel(files)
    folder=fullfile(outputDirectory,sprintf('Input_%03d',k));
    if ~isfolder(folder), mkdir(folder); end
    file=files(k); format=formats(k); fs=rates(k);
    if format=="ogg_iq"
        info=audioinfo(file); assert(info.NumChannels==2,'OGG IQ needs exactly two channels.');
        [audio,fs]=audioread(file); assert(all(isfinite(audio(:))));
        file=string(fullfile(folder,'Decoded_OGG_IQ.cf32'));
        fid=fopen(file,'wb','ieee-le'); assert(fid>=0); cleaner=onCleanup(@() fclose(fid));
        fwrite(fid,audio','single'); clear cleaner
        format="cf32";
    end
    if format=="audio"
        assert(~isempty(rows),'Scan the reference IQ first.');
        [r,scans{k}]=audioScan(file,k,rows(1,:),options.Threshold_dB,options.AudioFlagCorrelation);
    else
        part=arrange_scionx_bursts(file,strings(0,1),InputFormat=format,SampleRate=fs, ...
            ReferenceTime_s=options.ReferenceTime_s,Threshold_dB=options.Threshold_dB, ...
            OutputDirectory=folder,ExportPlots=false,OpenViewer=false);
        r=part.BurstResults; r.FileIndex(:)=k; scans{k}=part.Scans{1};
        if k==1, assert(~isempty(r),'No reference bursts found.'); end
    end
    rows=[rows;r]; %#ok<AGROW>
end
ref=find(rows.FileIndex==1); [~,pick]=min(abs(rows.Start_s(ref)-options.ReferenceTime_s));
scan=struct('BurstResults',rows,'ReferenceRow',ref(pick), ...
    'ComparisonRows',find(rows.FileIndex>1),'Scans',{scans});
writetable(rows,fullfile(outputDirectory,'Burst_Inventory.csv'));
result=scan; save(fullfile(outputDirectory,'Burst_Results.mat'),'result');
end

function [rows,scan]=audioScan(file,index,template,threshold,flagThreshold)
[audio,fs]=audioread(file); assert(size(audio,2)==1,'Discriminator OGG must be mono.');
assert(all(isfinite(audio)) && numel(audio)>fs*.1,'Invalid or too short audio.');
% Spectral energy removes slow baseline wander while preserving 9600-baud
% burst energy. One 10-ms record per block; median noise over 5 seconds.
n=max(128,round(.01*fs)); count=floor(numel(audio)/n);
window=.5-.5*cos(2*pi*(0:n-1)'/(n-1)); f=(0:n-1)'*fs/n; power=zeros(count,1);
for first=1:256:count
    last=min(count,first+255); blocks=reshape(audio((first-1)*n+1:last*n),n,[]);
    spectrum=abs(fft(blocks.*window)).^2;
    power(first:last)=sum(spectrum(f>=300 & f<=min(12000,fs/2-1),:),1)'/(n*sum(window.^2));
end
db=10*log10(max(power,realmin)); noise=movmedian(db,501,'Endpoints','shrink');
mask=db-noise>=threshold;
% FM discriminator audio can become QUIETER during a packet, so RF-style
% power-rise detection alone is insufficient. Add a normalized four-flag
% waveform search, independent of any callsign or telemetry contents.
flag=[0;1;1;1;1;1;1;0]; pattern=repmat(flag,4,1);
sampleIndex=floor((0:ceil(numel(pattern)*fs/9600)-1)'*9600/fs)+1;
sampleIndex=sampleIndex(sampleIndex<=numel(pattern));
flagTemplate=2*pattern(sampleIndex)-1; flagTemplate=flagTemplate-mean(flagTemplate);
m=numel(flagTemplate); ranked=zeros(0,2); maximumFlagCorrelation=0;
for first=1:262144:numel(audio)-m+1
    last=min(numel(audio)-m+1,first+262143); part=audio(first:last+m-1);
    fftLength=2^nextpow2(numel(part)+m-1);
    c=real(ifft(fft(part,fftLength).*fft(flipud(flagTemplate),fftLength)));
    c=c(m:numel(part)); s1=[0;cumsum(part)]; s2=[0;cumsum(part.^2)];
    v=max(0,s2(m+1:end)-s2(1:end-m)-(s1(m+1:end)-s1(1:end-m)).^2/m);
    score=abs(c)./max(sqrt(sum(flagTemplate.^2)*v),eps); hits=find(score>=flagThreshold);
    maximumFlagCorrelation=max(maximumFlagCorrelation,max(score));
    ranked=[ranked;hits+first-1 score(hits)]; %#ok<AGROW>
end
ranked=sortrows(ranked,-2); selected=[];
for j=1:size(ranked,1)
    position=ranked(j,1);
    if any(abs(selected-position)<.25*fs), continue; end
    selected(end+1)=position; %#ok<AGROW>
    lo=max(1,floor((position-.30*fs)/n)); hi=min(count,ceil((position+.40*fs)/n));
    mask(lo:hi)=true;
end
% Merge gaps <=30 ms, then include 20 ms guards for lossy audio edges.
[a,b]=runs(mask);
for j=1:numel(a)-1
    if a(j+1)-b(j)-1<=3, mask(b(j):a(j+1))=true; end
end
[a,b]=runs(mask); rows=template([],:); [~,name,ext]=fileparts(file);
for j=1:numel(a)
    if (b(j)-a(j)+1)*n/fs<.08, continue; end
    r=template; r.FileIndex=index; r.FilePath=file; r.FileLabel=string(name)+ext;
    r.BurstID=height(rows)+1; r.InputFormat="audio"; r.SampleRate=fs;
    r.BytesPerComplexSample=0; r.DisplayPolarity=1;
    r.StartSample=max(0,(a(j)-1)*n-round(.02*fs));
    r.EndSampleExclusive=min(numel(audio),b(j)*n+round(.02*fs));
    r.Start_s=r.StartSample/fs; r.End_s=r.EndSampleExclusive/fs; r.Duration_s=r.End_s-r.Start_s;
    r.PeakExcess_dB=max(db(a(j):b(j))-noise(a(j):b(j)));
    r.MedianExcess_dB=median(db(a(j):b(j))-noise(a(j):b(j)));
    r.NoiseBandPower_dBFS=median(noise(a(j):b(j)));
    r.TouchesFileEdge=r.StartSample==0 || r.EndSampleExclusive==numel(audio);
    r.ViewClipped=false; rows=[rows;r]; %#ok<AGROW>
end
scan=struct('FilePath',file,'SampleRate',fs,'Duration_s',numel(audio)/fs, ...
    'Time_s',((0:count-1)'+.5)*n/fs,'Power_dBFS',db,'Noise_dBFS',noise, ...
    'MaximumFlagCorrelation',maximumFlagCorrelation);
fprintf('Audio %s: %d burst candidates, %.0f samples/s, peak flag correlation %.3f\n', ...
    file,height(rows),fs,maximumFlagCorrelation);
end

function [a,b]=runs(mask)
d=diff([false;mask(:);false]); a=find(d==1); b=find(d==-1)-1;
end
