function result=selftest_decode_single_file(inputFile)
% Verify single-file isolation and persisted provenance with a real recording.
arguments
    inputFile (1,1) string
end
base=string(fileparts(mfilename('fullpath')));
inputFile=string(java.io.File(char(inputFile)).getCanonicalPath());
assert(isfile(inputFile),'Provide an existing comparison CS16 recording.');
previousPath=path; cleanup=onCleanup(@() path(previousPath)); %#ok<NASGU>
addpath(base);

% Reject other file types before any scan or output is created.
try
    decode_scionx_folder(fullfile(base,'run_decode_process.m'));
    error('SCIONX:Test:ExpectedRejection','A MATLAB script was accepted as CS16.');
catch err
    assert(strcmp(err.identifier,'SCIONX:Folder:InputFormat'),'%s',err.message);
end

% Recursive=true must not include siblings when the input is a file.
result=decode_scionx_folder(inputFile,Recursive=true,SampleRate=100000,ExportPlots=false);
assert(height(result.Recordings)==1 && height(result.Inputs)==2);
assert(height(result.FolderInventory)==1 && result.FolderInventory.Role=="comparison");
assert(result.FolderInventory.SourcePath==inputFile && result.Inputs.OriginalPath(2)==inputFile);
assert(result.SourcePath==inputFile && result.InputMode=="file");
assert(result.SourceDirectory==string(fileparts(inputFile)));
assert(all(result.Summary.FileIndex==2));
tlmFiles=result.Summary.TLMFile(result.Summary.TLMFile~="");
assert(all(isfile(tlmFiles)));
saved=load(fullfile(result.OutputDirectory,'Folder_Result.mat'),'result');
assert(saved.result.SourcePath==inputFile && saved.result.InputMode=="file");
inventory=readtable(fullfile(result.OutputDirectory,'Folder_Inventory.csv'),'TextType','string');
assert(height(inventory)==1 && inventory.SourcePath==inputFile);
fprintf('PASS: one comparison recording, no sibling scan, file-type rejection and saved provenance.\n');
fprintf('Bursts: %d; TLM files: %d; results: %s\n',height(result.Summary),numel(tlmFiles),result.OutputDirectory);
end
