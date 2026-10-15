function report=selftest_portable_toolkit()
% Decode both shipped samples and compare exact bytes with source baselines.
root=string(fileparts(mfilename('fullpath')));
integrity=verify_toolkit_integrity();
previousPath=path; cleanup=onCleanup(@() path(previousPath)); %#ok<NASGU>
process=fullfile(root,'SCIONX','封包頻譜分析','Decode_Process');
analysis=fileparts(process);
decoder=fullfile(root,'SCIONX','SignalSim','GNU_Matlab_Demo');
addpath(process,analysis,decoder,'-begin');
required=["decode_scionx_folder","analyze_scionx_packets","scan_scionx_pipeline_inputs", ...
    "arrange_scionx_bursts","align_scionx_postambles","detect_scionx_postamble", ...
    "select_scionx_postamble_candidate","compare_scionx_repeated_segments", ...
    "scionx_compare_bit_windows","decode_scionx_postamble_frame","scionx_reverse_hdlc", ...
    "scionx_crc16_x25","decode_scionx_partial_tlm","decode_scionx_tlm"];
for name=required
    resolved=string(which(name));
    assert(startsWith(resolved,root+filesep),'SCIONX:Portable:Path', ...
        'Function did not resolve inside this toolkit: %s -> %s',name,resolved);
end
expected=jsondecode(fileread(fullfile(root,'EXPECTED_RESULTS.json')));
cases=repmat(struct('Input',"",'OutputDirectory',"",'Bursts',0,'Status',""),numel(expected),1);
for k=1:numel(expected)
    e=expected(k); input=fullfile(root,string(e.InputRelativePath));
    r=decode_scionx_folder(input,SampleRate=100000,Recursive=true,ExportPlots=false);
    assert(height(r.Recordings)==1 && height(r.FolderInventory)==1 && r.InputMode=="file");
    assert(r.Control.CRC=="PASS" && r.Control.FrameBytes==274);
    assert(height(r.Summary)==numel(e.Packets));
    for j=1:numel(e.Packets)
        p=e.Packets(j); row=r.Summary(r.Summary.BurstID==p.BurstID,:);
        assert(height(row)==1 && row.FrameBytes==p.FrameBytes && row.Status==string(p.Status));
        assert(row.CRC==string(p.CRC) && row.DecodeSymbolRate==p.DecodeSymbolRate);
        assert(row.StopReason==string(p.StopReason) && abs(row.AnchorTime_s-p.AnchorTime_s)<1e-6);
        assert(isfile(row.TLMFile),'Missing TLM output.');
        received=fileread(fullfile(fileparts(row.TLMFile),'Received_Candidate.hex'));
        received=regexprep(upper(received),'\s+','');
        assert(strcmp(received,p.Hex),'SCIONX:Portable:Bytes','Decoded bytes differ from source baseline.');
    end
    cases(k)=struct('Input',input,'OutputDirectory',r.OutputDirectory, ...
        'Bursts',height(r.Summary),'Status',"PASS_EXACT_BYTES");
end
report=struct('Status',"PASS",'MatlabVersion',string(version),'Integrity',integrity,'Cases',cases);
out=fullfile(root,'Verification'); if ~isfolder(out), mkdir(out); end
fid=fopen(fullfile(out,'Portable_Test_Report.json'),'wt','n','UTF-8'); assert(fid>=0);
closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report));
fprintf('PASS: both relocated samples match the source bytes and packet statuses.\n');
end
