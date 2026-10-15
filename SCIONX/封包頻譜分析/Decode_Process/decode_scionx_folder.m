function result=decode_scionx_folder(inputDirectory,options)
%DECODE_SCIONX_FOLDER Decode CS16 recordings using July 23 postamble alignment.
% result=decode_scionx_folder("D:\captures")
% result=decode_scionx_folder("D:\captures",Recursive=true,SampleRate=100000)
% result=decode_scionx_folder("D:\captures\recording.cs16",SampleRate=100000)
% Accept a folder or one CS16 file. Recursive only applies to folders.
% With no input, show a folder picker. Outputs always go under ./Results.
arguments
    inputDirectory (1,1) string = ""
    options.Recursive (1,1) logical = false
    options.SampleRate (1,1) double = NaN
    options.ExportPlots (1,1) logical = true
    options.ResolveRepeatedPostambles (1,1) logical = true
    options.DetectionThreshold_dB (1,1) double {mustBeFinite,mustBeNonnegative} = 5
end
base=string(fileparts(mfilename('fullpath'))); analysisDirectory=fileparts(base);
if inputDirectory==""
    selected=uigetdir(char(fullfile(analysisDirectory,'..','Beacon Data','CS16_Raw')), ...
        'Select a folder containing CS16 recordings');
    if isequal(selected,0)
        error('SCIONX:Folder:Cancelled','Folder selection cancelled; no recordings processed.');
    end
    inputDirectory=string(selected);
end
assert(isfolder(inputDirectory) || isfile(inputDirectory),'SCIONX:Folder:Missing', ...
    'Input folder or file not found: %s',inputDirectory);
inputPath=canonical(inputDirectory);
singleFile=isfile(inputPath);
inputMode="folder";
if singleFile
    [inputDirectory,~,extension]=fileparts(inputPath);
    assert(strcmpi(extension,'.cs16'),'SCIONX:Folder:InputFormat', ...
        'Single-file input must have a .cs16 extension: %s',inputPath);
    inputMode="file";
else
    inputDirectory=inputPath;
end
reference=canonical(fullfile(analysisDirectory,'..','Beacon Data','CS16_Raw','20260723', ...
    '20260723_091639_100000_98266_x.cs16'));
assert(isfile(reference),'SCIONX:Folder:Reference','July 23 reference not found: %s',reference);
assert(isnan(options.SampleRate) || (isfinite(options.SampleRate) && options.SampleRate>24000), ...
    'SCIONX:Folder:SampleRate','SampleRate must exceed 24000 samples/s, or be NaN for filename inference.');

if singleFile
    listing=dir(inputPath);
else
    pattern=fullfile(inputDirectory,'*');
    if options.Recursive, pattern=fullfile(inputDirectory,'**','*'); end
    listing=dir(pattern);
end
listing=listing(~[listing.isdir]);
names=string({listing.name})';
listing=listing(endsWith(names,'.cs16','IgnoreCase',true));
assert(~isempty(listing),'SCIONX:Folder:NoCS16','No .cs16 recordings found in: %s',inputDirectory);
paths=strings(numel(listing),1);
for k=1:numel(listing), paths(k)=canonical(fullfile(listing(k).folder,listing(k).name)); end
[~,order]=sort(lower(paths)); paths=paths(order); listing=listing(order);
isReference=strcmpi(paths,reference);
isMetadata=false(numel(paths),1);
for k=1:numel(paths)
    if startsWith(string(listing(k).name),'._')
        fid=fopen(paths(k),'rb','ieee-be');
        assert(fid>=0,'SCIONX:Folder:Read','Cannot read: %s',paths(k));
        closer=onCleanup(@() fclose(fid));
        magic=fread(fid,1,'uint32=>uint32'); clear closer
        isMetadata(k)=isequal(magic,uint32(hex2dec('00051607')));
    end
end
excluded=isReference | isMetadata;
files=paths(~excluded);
assert(~isempty(files),'SCIONX:Folder:ReferenceOnly', ...
    'No comparison recordings remain after excluding the reference and AppleDouble metadata.');
bytes=reshape([listing.bytes],[],1);
rates=nan(numel(paths),1); rates(isReference)=100000;
for k=find(~excluded)'
    assert(bytes(k)>0 && mod(bytes(k),4)==0,'SCIONX:Folder:InvalidCS16', ...
        'CS16 must contain complete int16 I/Q pairs (4 bytes per sample): %s',paths(k));
    rates(k)=options.SampleRate;
    if isnan(rates(k))
        [~,stem]=fileparts(paths(k));
        token=regexp(stem,'_(\d+)SPS_','tokens','once');
        if isempty(token), token=regexp(stem,'^\d{8}_\d{6}_(\d+)_','tokens','once'); end
        assert(~isempty(token),'SCIONX:Folder:UnknownSampleRate', ...
            'Cannot infer sample rate for %s. Supply SampleRate explicitly.',paths(k));
        rates(k)=str2double(token{1});
    end
    assert(isfinite(rates(k)) && rates(k)>24000,'SCIONX:Folder:SampleRate', ...
        'Unsupported sample rate in %s: %.9g',paths(k),rates(k));
end
roles=repmat("comparison",numel(paths),1);
roles(isReference)="reference_excluded"; roles(isMetadata)="appledouble_excluded";
inventory=table(paths,bytes,rates,roles,'VariableNames', ...
    {'SourcePath','Bytes','SampleRate','Role'});
fprintf('CS16 input (%s): %s\nComparison recordings: %d; reference excluded: %d; metadata excluded: %d\n', ...
    inputMode,inputPath,numel(files),sum(isReference),sum(isMetadata));
fprintf('Reference: %s\nPolicy: decode backwards from resolved postamble anchors.\n',reference);
previousPath=path; cleanup=onCleanup(@() path(previousPath)); %#ok<NASGU>
addpath(analysisDirectory,'-begin');
result=analyze_scionx_packets(files,ReferenceFile=reference, ...
    ReferenceSampleRate=100000,ReferenceTime_s=306.6,InputFormat="cs16", ...
    SampleRate=rates(~excluded),OutputDirectory=fullfile(base,'Results'), ...
    DecodePolicy="postamble",ExportPlots=options.ExportPlots, ...
    ResolveRepeatedPostambles=options.ResolveRepeatedPostambles, ...
    DetectionThreshold_dB=options.DetectionThreshold_dB);
result.SourceDirectory=inputDirectory;
result.SourcePath=inputPath;
result.InputMode=inputMode;
result.FolderOptions=options;
result.FolderInventory=inventory;
writetable(inventory,fullfile(result.OutputDirectory,'Folder_Inventory.csv'));
save(fullfile(result.OutputDirectory,'Folder_Result.mat'),'result');
save(fullfile(base,'Results','Latest_Result.mat'),'result');
fprintf('CS16 decoding complete: %s\n',result.OutputDirectory);
end

function path=canonical(path)
path=string(java.io.File(char(path)).getCanonicalPath());
end
