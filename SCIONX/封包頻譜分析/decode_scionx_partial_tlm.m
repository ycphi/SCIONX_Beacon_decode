function outputFile = decode_scionx_partial_tlm(candidateFile, firstByte, outputFile, decoderPath, frameLength)
%DECODE_SCIONX_PARTIAL_TLM Map a tail-aligned candidate to diagnostic TLM.
% candidateFile: pure ASCII hex bytes, ending with the candidate two-byte FCS.
% firstByte: assumed 1-based position in the full frame (e.g. 47 or 204).
% decoderPath: path to the user's original decode_scionx_tlm.m.
% Missing fields are redacted. Frame positions are ASSUMED, not validated.
% The original decoder is copied into an isolated temporary work directory.
% Internal zero placeholders enable its field conversions but are NEVER
% reported as received data. They are not used for a claimed CRC result.

    if nargin < 5, frameLength = 274; end
    if nargin < 4 || isempty(decoderPath)
        decoderPath = which('decode_scionx_tlm');
    end
    if isempty(decoderPath) || ~isfile(decoderPath)
        error('SCIONX:DecoderMissing','Provide the original decode_scionx_tlm.m path.');
    end
    sourceText = strtrim(fileread(candidateFile));
    if isempty(regexp(sourceText,'^(?:[0-9A-Fa-f]{2}\s*)+$','once'))
        error('SCIONX:HexOnly','Candidate input must contain only two-digit hex bytes and whitespace.');
    end
    tokens = regexp(sourceText,'[0-9A-Fa-f]{2}','match');
    bytes = uint8(hex2dec(tokens(:))).';
    validateattributes(firstByte,{'numeric'},{'scalar','integer','positive'});
    validateattributes(frameLength,{'numeric'},{'scalar','integer','>=',274});
    if numel(bytes)<2 || firstByte + numel(bytes)-1 ~= frameLength
        error('SCIONX:TailPosition','The candidate must end at the assumed final FCS byte.');
    end
    outputDir = fileparts(outputFile);
    if isempty(outputDir), outputDir = pwd; end
    if ~isfolder(outputDir), mkdir(outputDir); end
    workDir = tempname(outputDir);
    mkdir(workDir);
    previousPath = path;
    cleanup = onCleanup(@() cleanupWork(workDir,outputDir,previousPath));
    copyfile(decoderPath,fullfile(workDir,'decode_scionx_tlm.m'));
    addpath(workDir,'-begin');
    clear decode_scionx_tlm;
    mask = false(1,frameLength);
    mask(firstByte:end)=true;
    internal = zeros(1,frameLength,'uint8');
    internal(mask)=bytes;
    rawFile=fullfile(workDir,'Raw_bit.txt');
    fid=fopen(rawFile,'wt');
    if fid<0, error('SCIONX:Write','Cannot create temporary input.'); end
    fprintf(fid,'%02X ',internal); fprintf(fid,'\n'); fclose(fid);
    fid=fopen(fullfile(workDir,'Matlab_HDLC_Status.txt'),'wt');
    fprintf(fid,'Boundary Mode: REVERSE_PARTIAL_DIAGNOSTIC\n'); fclose(fid);
    decode_scionx_tlm(rawFile);
    generated = splitlines(string(fileread(fullfile(workDir,'TLM.txt'))));
    lines = {
        'SCION-X reverse-HDLC partial telemetry - DIAGNOSTIC ONLY'
        ['Input candidate = ' char(candidateFile)]
        ['Reference decoder = ' char(decoderPath)]
        sprintf('Assumed full frame length = %d bytes, including FCS',frameLength)
        sprintf('Received candidate suffix = %d bytes; assumed positions %d-%d',numel(bytes),firstByte,frameLength)
        'Boundary Mode = REVERSE_PARTIAL_DIAGNOSTIC'
        'CRC Check = UNVERIFIABLE_PARTIAL (full frame unavailable)'
        'Frame Validation = NOT VALIDATED'
        'Telemetry Use = MANUAL DIAGNOSTIC REVIEW ONLY'
        'Byte numbering = 1-based AX.25 positions, inferred from the assumed closing boundary'
        'Integer byte order and conversions follow the supplied decoder'
        'MISSING = at least one byte of this field was not received; no value is inferred'
        'Present bytes may still contain bit/stuffing/alignment errors; plausible values do not validate them'
        ''
    };
    inFields=false; suppress=false;
    for k=1:numel(generated)
        line=char(generated(k));
        if startsWith(line,'[AX.25 / HDLC]'), inFields=true; end
        if ~inFields, continue; end
        if startsWith(line,'[Frame FCS]'), break; end
        if startsWith(line,'['), suppress=false; end
        match=regexp(line,'^\s*Bytes?\s+(\d+(?:-\d+)?)\s+(.+?)\s*=','tokens','once');
        if ~isempty(match)
            positions=sscanf(strrep(match{1},'-',' '),'%d');
            lo=positions(1); hi=positions(end);
            suppress = hi>frameLength || any(~mask(lo:hi));
            if suppress
                label=regexprep(line,'\s*=.*$','');
                lines{end+1}=[label ' = MISSING (field incomplete)']; %#ok<AGROW>
                continue;
            end
        end
        if ~suppress, lines{end+1}=line; end %#ok<AGROW>
    end
    lines{end+1}='';
    lines{end+1}='[Candidate Frame FCS]';
    lines{end+1}=sprintf('Assumed bytes %d-%d = %02X %02X',frameLength-1,frameLength,bytes(end-1),bytes(end));
    lines{end+1}='CRC remains UNVERIFIABLE_PARTIAL; no missing prefix was used to validate it.';
    lines{end+1}='';
    lines{end+1}='[Received candidate bytes - no missing-byte padding]';
    for k=1:16:numel(bytes)
        last=min(numel(bytes),k+15);
        lines{end+1}=sprintf('Assumed bytes %03d-%03d: %s',firstByte+k-1,firstByte+last-1,strtrim(sprintf('%02X ',bytes(k:last)))); %#ok<AGROW>
    end
    fid=fopen(outputFile,'wt','n','UTF-8');
    if fid<0, error('SCIONX:Write','Cannot create output file.'); end
    closeOutput=onCleanup(@() fclose(fid));
    for k=1:numel(lines), fprintf(fid,'%s\n',lines{k}); end
    fprintf('Partial diagnostic TLM: %s\n',outputFile);
end

function cleanupWork(workDir,outputDir,previousPath)
    path(previousPath);
    clear decode_scionx_tlm;
    % Verify resolved absolute paths before recursive deletion.
    root=char(java.io.File(outputDir).getCanonicalPath());
    target=char(java.io.File(workDir).getCanonicalPath());
    if startsWith(lower(target),[lower(root) filesep]) && isfolder(target)
        rmdir(target,'s');
    end
end

