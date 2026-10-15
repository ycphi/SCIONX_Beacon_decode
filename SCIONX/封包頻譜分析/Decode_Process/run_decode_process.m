% Portable CS16 entry point. See the toolkit README_zh-TW.md.
% Set inputDirectory to a folder or one .cs16 file; empty opens a folder picker.
% Edit inputDirectory below, then press Run on THIS script.
% Example: inputDirectory="E:\captures\recording.cs16";
% Example: inputDirectory="E:\captures";
processDirectory=string(fileparts(mfilename('fullpath')));
scionxDirectory=fileparts(fileparts(processDirectory));
inputDirectory=fullfile(scionxDirectory,'Beacon Data','CS16_Raw','20260720','20260720_220206_100000_98266_x.cs16');
recursive=true;       % Include nested folders; ignored for a single file.
sampleRate=100000;    % Override all comparison rates; NaN infers from filenames.
exportPlots=true;     % Export postamble alignment and repetition diagnostic plots.
resolveRepeatedPostambles=true; % Earlier flags require strong repeated-data evidence.

addpath(processDirectory);
decodeResult=decode_scionx_folder(inputDirectory,Recursive=recursive, ...
    SampleRate=sampleRate,ExportPlots=exportPlots, ...
    ResolveRepeatedPostambles=resolveRepeatedPostambles);
% Outputs: this Decode_Process folder / Results / Run_...
