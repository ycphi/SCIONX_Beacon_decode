function results = run_scionx_pipeline(inputPath, outputRoot, radioCondaRoot, parallelWorkers)
%RUN_SCIONX_PIPELINE Run CS16 IQ or OGG audio -> Python -> HDLC -> TLM.
%
% results = run_scionx_pipeline()
%   Opens a dialog to select one CS16/OGG recording or a folder.
%
% results = run_scionx_pipeline(inputPath)
%   inputPath may be one .cs16/.ogg file or a folder. Folder processing is
%   recursive and includes both supported formats.
%
% results = run_scionx_pipeline(inputPath, outputRoot, radioCondaRoot, parallelWorkers)
%   outputRoot defaults to <this demo>/Pipeline_Output.
%   The third argument accepts a Python executable or virtual-environment
%   folder. This portable package defaults to <this demo>/runtime/python.exe.
%   A RadioConda installation folder still selects the original GRC frontend.
%   parallelWorkers controls file-level process parallelism:
%     [] or 0 = automatic (up to 2 workers for a folder);
%     1       = serial processing;
%     N > 1   = at most N workers (and never more than the file count).
%   If Parallel Computing Toolbox or a process pool is unavailable, the
%   pipeline warns and automatically falls back to serial processing.
%
% The default portable frontend needs NumPy, SciPy and SoundFile, but no GNU
% Radio/RadioConda. The portable package includes Python and all dependencies.

    demoDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(inputPath)
        inputPath = chooseInputPath();
        if isempty(inputPath)
            results = struct([]);
            fprintf('SCION-X pipeline cancelled by user.\n');
            return;
        end
    end
    if nargin < 2 || isempty(outputRoot)
        outputRoot = fullfile(demoDir, 'Pipeline_Output');
    end
    if nargin < 3 || isempty(radioCondaRoot)
        radioCondaRoot = '';
    end
    if nargin < 4 || isempty(parallelWorkers)
        parallelWorkers = 0;
    end
    validateattributes(parallelWorkers, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, ...
        mfilename, 'parallelWorkers', 4);

    inputPath = char(inputPath);
    outputRoot = char(outputRoot);
    runtime = resolveDecoderRuntime(demoDir, char(radioCondaRoot));
    generatedCs16Python = fullfile(demoDir, 'SCIONX_IQ_BurstScan_CS16.py');
    generatedOggPython = fullfile(demoDir, 'SCIOX_IQ_BurstScan_OGG.py');
    if runtime.portable
        generatedCs16Python = fullfile(demoDir, 'scionx_portable_scan.py');
        generatedOggPython = generatedCs16Python;
        requireFile(generatedCs16Python, 'Portable Python entry point is missing.');
        [checkStatus, checkOutput] = system(strjoin({runtime.prefix, ...
            quoteCommandArg(generatedCs16Python), '--check-environment'}, ' '));
        if checkStatus ~= 0
            error('SCIONX:PythonEnvironment', [ ...
                'Portable Python environment is incomplete. Re-extract the complete ZIP ' ...
                'including runtime and soft_repair, or pass a configured python.exe as argument 3.\n%s'], ...
                checkOutput);
        end
    end
    requireFile(generatedCs16Python, ...
        'Generate SCIONX_IQ_BurstScan_CS16.py from the CS16 GRC flowgraph first.');
    requireFile(generatedOggPython, ...
        'Generate SCIOX_IQ_BurstScan_OGG.py from the OGG GRC flowgraph first.');
    fprintf('SCION-X decoder runtime: %s\n', runtime.description);

    inputFiles = collectRecordingFiles(inputPath);
    if isempty(inputFiles)
        error('SCIONX:PipelineNoInput', ...
            'No .cs16 or .ogg recordings found: %s', inputPath);
    end
    if ~isfolder(outputRoot)
        mkdir(outputRoot);
    end

    resultTemplate = struct( ...
        'input_file', '', ...
        'input_format', '', ...
        'decoder_backend', '', ...
        'output_directory', '', ...
        'gnu_exit_code', -1, ...
        'frame_count', 0, ...
        'crc_pass_count', 0, ...
        'crc_fail_count', 0, ...
        'boundary_recovered_count', 0, ...
        'diagnostic_tlm_count', 0, ...
        'soft_repair_attempted', false, ...
        'soft_repaired_count', 0, ...
        'soft_repair_statuses', '', ...
        'soft_repair_reasons', '', ...
        'tlm_count', 0, ...
        'parallel_workers', 1, ...
        'elapsed_seconds', 0, ...
        'status', '', ...
        'no_tlm_reason', '', ...
        'message', '');
    results = repmat(resultTemplate, 1, numel(inputFiles));

    workerCount = resolveParallelWorkers(parallelWorkers, numel(inputFiles));
    fprintf('SCION-X pipeline: %d recording file(s), workers=%d\n', ...
        numel(inputFiles), workerCount);
    if workerCount > 1
        parfor (fileIndex = 1:numel(inputFiles), workerCount)
            results(fileIndex) = processRecording( ...
                fileIndex, numel(inputFiles), inputFiles{fileIndex}, ...
                outputRoot, runtime, ...
                generatedCs16Python, generatedOggPython, ...
                resultTemplate, workerCount, false);
        end
    else
        for fileIndex = 1:numel(inputFiles)
            results(fileIndex) = processRecording( ...
                fileIndex, numel(inputFiles), inputFiles{fileIndex}, ...
                outputRoot, runtime, ...
                generatedCs16Python, generatedOggPython, ...
                resultTemplate, workerCount, true);
        end
    end

    summaryFile = fullfile(outputRoot, 'Pipeline_Summary.csv');
    resultTable = struct2table(results, 'AsArray', true);
    writetable(resultTable, summaryFile);
    fprintf('\nSCION-X pipeline complete: %s\n', summaryFile);
    disp(resultTable);
end


function result = processRecording(fileIndex, fileCount, inputFile, ...
        outputRoot, runtime, generatedCs16Python, ...
        generatedOggPython, resultTemplate, workerCount, showLiveOutput)
    started = tic;
    result = resultTemplate;
    if runtime.portable
        result.decoder_backend = 'PORTABLE_PYTHON';
    else
        result.decoder_backend = 'GNU_RADIO';
    end
    result.parallel_workers = workerCount;
    [~, recordingName, extension] = fileparts(inputFile);
    if strcmpi(extension, '.cs16')
        inputFormat = 'CS16_IQ';
        generatedPython = generatedCs16Python;
        inputOption = '--iq-file';
    else
        inputFormat = 'OGG_AUDIO';
        generatedPython = generatedOggPython;
        inputOption = '--ogg-file';
    end
    recordingName = regexprep(recordingName, '[^A-Za-z0-9_.-]', '_');
    recordingDir = fullfile(outputRoot, ...
        sprintf('%03d_%s', fileIndex, recordingName));
    framesDir = fullfile(recordingDir, 'Frames');
    gnuHexTemplate = fullfile(framesDir, 'GNU_GFSK_Raw_HEX.txt');
    gnuReport = fullfile(recordingDir, 'GNU_Burst_Scan_Report.csv');
    matlabReport = fullfile(recordingDir, 'Matlab_Frame_Report.csv');
    if ~isfolder(recordingDir)
        mkdir(recordingDir);
    end

    result.input_file = inputFile;
    result.input_format = inputFormat;
    result.output_directory = recordingDir;
    fprintf('\n[%d/%d] %s BurstScan: %s\n', ...
        fileIndex, fileCount, inputFormat, inputFile);

    command = strjoin({ ...
        runtime.prefix, quoteCommandArg(generatedPython), ...
        inputOption, quoteCommandArg(inputFile), ...
        '--output-txt-path', quoteCommandArg(gnuHexTemplate), ...
        '--report-csv-path', quoteCommandArg(gnuReport)}, ' ');
    if showLiveOutput
        [gnuStatus, gnuOutput] = system(command, '-echo');
    else
        [gnuStatus, gnuOutput] = system(command);
    end
    result.gnu_exit_code = gnuStatus;
    if gnuStatus ~= 0
        if runtime.portable
            result.status = 'PYTHON_FAILED';
        else
            result.status = 'GNU_RADIO_FAILED';
        end
        result.message = strtrim(gnuOutput);
        result.no_tlm_reason = [result.status ': ' compactText(gnuOutput)];
        result.elapsed_seconds = toc(started);
        warning('SCIONX:PipelineGnuFailed', ...
            'Python decoder failed for %s (exit %d). Continuing batch.', ...
            inputFile, gnuStatus);
        return;
    end

    [softAttempted, softRepairedCount, softStatuses, softReasons] = ...
        readGnuSoftRepairSummary(gnuReport);
    result.soft_repair_attempted = softAttempted;
    result.soft_repaired_count = softRepairedCount;
    result.soft_repair_statuses = softStatuses;
    result.soft_repair_reasons = softReasons;

    gnuFrameFiles = dir(fullfile(framesDir, 'Frame_*', 'GNU_GFSK_Raw_HEX.txt'));
    if isempty(gnuFrameFiles)
        result.status = 'NO_COMPLETE_FRAME_CANDIDATES';
        result.message = ...
            'Burst scan completed, but no complete HDLC frame candidate was produced.';
        result.no_tlm_reason = summarizeNoTlmReason(gnuReport, '', ...
            'PYTHON_NO_RAW_HDLC_CANDIDATE');
        result.elapsed_seconds = toc(started);
        fprintf('[%d/%d] No complete HDLC frame candidate: %s\n', ...
            fileIndex, fileCount, inputFile);
        return;
    end

    try
        if showLiveOutput
            decode_gnuradio_hex(framesDir, framesDir, matlabReport);
        else
            evalc('decode_gnuradio_hex(framesDir, framesDir, matlabReport);');
        end
        rawFiles = dir(fullfile(framesDir, 'Frame_*', 'Raw_bit.txt'));
        if isempty(rawFiles)
            result.status = 'NO_MATLAB_HDLC_FRAME';
            result.message = ...
                'Python produced raw candidates, but MATLAB could not recover a complete 274-byte HDLC frame.';
            result.no_tlm_reason = summarizeNoTlmReason(gnuReport, matlabReport, ...
                'MATLAB_NO_274_BYTE_HDLC_FRAME');
            result.elapsed_seconds = toc(started);
            fprintf('[%d/%d] No MATLAB HDLC frame: %s\n', ...
                fileIndex, fileCount, inputFile);
            return;
        end
        if showLiveOutput
            decode_scionx_tlm(framesDir);
        else
            evalc('decode_scionx_tlm(framesDir);');
        end
        tlmFiles = dir(fullfile(framesDir, 'Frame_*', 'TLM.txt'));
        result.frame_count = numel(rawFiles);
        result.tlm_count = numel(tlmFiles);
        [crcPassCount, crcFailCount, recoveredCount, diagnosticCount] = ...
            countTlmResults(tlmFiles);
        result.crc_pass_count = crcPassCount;
        result.crc_fail_count = crcFailCount;
        result.boundary_recovered_count = recoveredCount;
        result.diagnostic_tlm_count = diagnosticCount;
        if isempty(tlmFiles)
            result.no_tlm_reason = 'TLM_DECODER_PRODUCED_NO_OUTPUT';
        elseif numel(tlmFiles) < numel(gnuFrameFiles)
            result.no_tlm_reason = summarizeNoTlmReason(gnuReport, matlabReport, ...
                sprintf('%d_OF_%d_RAW_CANDIDATES_DID_NOT_PRODUCE_TLM', ...
                numel(gnuFrameFiles) - numel(tlmFiles), numel(gnuFrameFiles)));
        end
        if softRepairedCount > 0 && crcPassCount == numel(tlmFiles)
            result.status = 'SUCCESS_WITH_SOFT_REPAIR';
        elseif softRepairedCount > 0
            result.status = 'COMPLETED_WITH_SOFT_REPAIR';
        elseif diagnosticCount == numel(tlmFiles) && diagnosticCount > 0
            result.status = 'DIAGNOSTIC_TLM_ONLY';
        elseif diagnosticCount > 0
            result.status = 'COMPLETED_WITH_DIAGNOSTIC_TLM';
        elseif crcFailCount > 0
            result.status = 'COMPLETED_WITH_CRC_FAILURE';
        elseif recoveredCount > 0
            result.status = 'SUCCESS_WITH_BOUNDARY_RECOVERY';
        else
            result.status = 'SUCCESS';
        end
        result.message = sprintf( ...
            ['%d frame(s), CRC PASS=%d, CRC FAIL=%d, ' ...
             'Boundary RECOVERED=%d, Diagnostic TLM=%d, %d TLM file(s), ' ...
             'Soft repair attempted=%d, repaired=%d, status=%s'], ...
            numel(rawFiles), crcPassCount, crcFailCount, ...
            recoveredCount, diagnosticCount, numel(tlmFiles), ...
            softAttempted, softRepairedCount, softStatuses);
    catch decodeError
        result.status = 'MATLAB_DECODE_FAILED';
        result.message = decodeError.message;
        result.no_tlm_reason = ['MATLAB_DECODE_FAILED: ' compactText(decodeError.message)];
        warning('SCIONX:PipelineMatlabFailed', ...
            'MATLAB decode failed for %s: %s. Continuing batch.', ...
            inputFile, decodeError.message);
    end
    result.elapsed_seconds = toc(started);
    fprintf('[%d/%d] Complete in %.1f s: %s (%s)\n', ...
        fileIndex, fileCount, result.elapsed_seconds, inputFile, result.status);
end


function workerCount = resolveParallelWorkers(requested, fileCount)
    workerCount = 1;
    if fileCount < 2 || requested == 1
        return;
    end
    if ~license('test', 'Distrib_Computing_Toolbox')
        if requested > 1
            warning('SCIONX:ParallelUnavailable', ...
                'Parallel Computing Toolbox is unavailable; using serial processing.');
        end
        return;
    end

    if requested == 0
        requested = min(2, fileCount);
    else
        requested = min(double(requested), fileCount);
    end
    if requested < 2
        return;
    end

    try
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('Processes', requested);
        end
        workerCount = min([requested, pool.NumWorkers, fileCount]);
        if workerCount < requested
            warning('SCIONX:ParallelPoolSmaller', ...
                'Requested %d workers, but the active pool has %d; using %d.', ...
                requested, pool.NumWorkers, workerCount);
        end
    catch poolError
        workerCount = 1;
        warning('SCIONX:ParallelPoolFailed', ...
            'Could not start a process pool (%s); using serial processing.', ...
            poolError.message);
    end
end


function inputPath = chooseInputPath()
    choice = questdlg('Select SCION-X input type', 'SCION-X Pipeline', ...
        'Single recording', 'Recording folder', 'Cancel', 'Single recording');
    switch choice
        case 'Single recording'
            [filename, pathname] = uigetfile( ...
                {'*.cs16;*.ogg', 'SCION-X recordings (*.cs16, *.ogg)'; ...
                 '*.cs16', 'CS16 IQ recordings (*.cs16)'; ...
                 '*.ogg', 'OGG audio recordings (*.ogg)'; ...
                 '*.*', 'All files'}, ...
                'Select a SCION-X recording');
            if isequal(filename, 0)
                inputPath = '';
            else
                inputPath = fullfile(pathname, filename);
            end
        case 'Recording folder'
            selected = uigetdir(pwd, ...
                'Select a folder containing CS16 or OGG recordings');
            if isequal(selected, 0)
                inputPath = '';
            else
                inputPath = selected;
            end
        otherwise
            inputPath = '';
    end
end


function files = collectRecordingFiles(inputPath)
    if isfile(inputPath)
        [~, inputName, extension] = fileparts(inputPath);
        if ~ismember(lower(extension), {'.cs16', '.ogg'})
            error('SCIONX:PipelineWrongExtension', ...
                'Input file must have the .cs16 or .ogg extension: %s', inputPath);
        end
        if startsWith(inputName, '._')
            error('SCIONX:PipelineResourceFork', ...
                'macOS resource-fork files (._*) are not recordings: %s', inputPath);
        end
        files = {char(inputPath)};
        return;
    end
    if ~isfolder(inputPath)
        error('SCIONX:PipelineInputMissing', 'Input path does not exist: %s', inputPath);
    end
    listing = [dir(fullfile(inputPath, '**', '*.cs16')); ...
               dir(fullfile(inputPath, '**', '*.ogg'))];
    listing = listing(~startsWith({listing.name}, '._'));
    files = arrayfun(@(item) fullfile(item.folder, item.name), ...
        listing, 'UniformOutput', false);
    files = sort(files);
end


function requireFile(filename, hint)
    if ~isfile(filename)
        error('SCIONX:PipelineRequiredFile', '%s Missing: %s', hint, filename);
    end
end


function runtime = resolveDecoderRuntime(demoDir, requested)
    runtime = struct('prefix', '', 'portable', true, 'description', '');
    if isempty(requested)
        candidates = {fullfile(demoDir, 'runtime', 'python.exe'), ...
            getenv('SCIONX_PYTHON'), ...
            fullfile(demoDir, '.venv', 'Scripts', 'python.exe'), ...
            fullfile(demoDir, 'runtime', 'python.exe')};
        if ~isempty(getenv('VIRTUAL_ENV'))
            candidates{end + 1} = fullfile(getenv('VIRTUAL_ENV'), 'Scripts', 'python.exe');
        end
        for index = 1:numel(candidates)
            if ~isempty(candidates{index}) && isfile(candidates{index})
                requested = candidates{index};
                break;
            end
        end
        if isempty(requested)
            [status, found] = system('where python.exe');
            if status == 0
                paths = regexp(strtrim(found), '\r?\n', 'split');
                for index = 1:numel(paths)
                    if isfile(paths{index}) && ~contains(paths{index}, 'WindowsApps')
                        requested = paths{index};
                        break;
                    end
                end
            end
        end
        if isempty(requested)
            [status, found] = system('py -3.12 -c "import sys; print(sys.executable)"');
            if status == 0 && isfile(strtrim(found))
                requested = strtrim(found);
            end
        end
        if isempty(requested)
            requested = discoverRadioCondaRoot();
        end
    end
    if isempty(requested)
        error('SCIONX:PythonNotFound', [ ...
            'No Python runtime found. Install 64-bit Python 3.12 and run ' ...
            'setup_python.ps1 in the package folder. RadioConda is NOT required. ' ...
            'You can also pass python.exe as argument 3 or set SCIONX_PYTHON.']);
    end
    if isfolder(requested) && isfile(fullfile(requested, 'Scripts', 'conda.exe'))
        runtime.portable = false;
        runtime.prefix = strjoin({quoteCommandArg(fullfile(requested, 'Scripts', 'conda.exe')), ...
            'run', '-p', quoteCommandArg(requested), 'python'}, ' ');
        runtime.description = ['GNU_RADIO: ' requested];
        return;
    end
    if isfolder(requested)
        if isfile(fullfile(requested, 'Scripts', 'python.exe'))
            requested = fullfile(requested, 'Scripts', 'python.exe');
        else
            requested = fullfile(requested, 'python.exe');
        end
    end
    requireFile(requested, 'Pass the full path of python.exe, or run setup_python.ps1.');
    runtime.prefix = quoteCommandArg(requested);
    runtime.description = ['PORTABLE_PYTHON: ' requested];
end


function root = discoverRadioCondaRoot()
    root = '';
    candidates = strings(0, 1);
    configured = string(getenv('RADIOCONDA_ROOT'));
    if strlength(configured) > 0
        candidates(end + 1) = configured; %#ok<AGROW>
    end
    activeConda = string(getenv('CONDA_PREFIX'));
    if strlength(activeConda) > 0
        candidates(end + 1) = activeConda; %#ok<AGROW>
    end
    userProfile = string(getenv('USERPROFILE'));
    if strlength(userProfile) > 0
        candidates(end + 1) = fullfile(userProfile, 'radioconda'); %#ok<AGROW>
    end
    candidates(end + 1) = "C:\radioconda";
    candidates = unique(candidates, 'stable');

    for index = 1:numel(candidates)
        candidate = char(candidates(index));
        if isfile(fullfile(candidate, 'Scripts', 'conda.exe')) && ...
                isfile(fullfile(candidate, 'python.exe'))
            root = candidate;
            fprintf('RadioConda automatically detected: %s\n', root);
            return;
        end
    end
end


function [passCount, failCount, recoveredCount, diagnosticCount] = countTlmResults(tlmFiles)
    passCount = 0;
    failCount = 0;
    recoveredCount = 0;
    diagnosticCount = 0;
    for index = 1:numel(tlmFiles)
        content = fileread(fullfile(tlmFiles(index).folder, tlmFiles(index).name));
        if ~isempty(regexp(content, '(?m)^CRC Check = PASS\s*$', 'once'))
            passCount = passCount + 1;
        elseif ~isempty(regexp(content, '(?m)^CRC Check = FAIL\s*$', 'once'))
            failCount = failCount + 1;
        end
        if ~isempty(regexp(content, ...
                '(?m)^Boundary Mode = (?:RECOVERED|MATLAB_RECOVERED)\s*$', 'once'))
            recoveredCount = recoveredCount + 1;
        end
        if ~isempty(regexp(content, ...
                '(?m)^Boundary Mode = (?:FLAGS_ONLY_DIAGNOSTIC|MATLAB_FLAGS_ONLY_DIAGNOSTIC)\s*$', 'once'))
            diagnosticCount = diagnosticCount + 1;
        end
    end
end


function reason = summarizeNoTlmReason(gnuReport, matlabReport, prefix)
    parts = string(prefix);
    matlabReasons = reportValues(matlabReport, {'status', 'boundary_mode'});
    if ~isempty(matlabReasons)
        parts(end + 1) = "MATLAB=" + strjoin(matlabReasons, '|');
    end
    gnuReasons = reportValues(gnuReport, ...
        {'status', 'soft_repair_status', 'soft_repair_reason'});
    if ~isempty(gnuReasons)
        parts(end + 1) = "PYTHON=" + strjoin(gnuReasons, '|');
    end
    reason = char(strjoin(parts, '; '));
end


function values = reportValues(reportFile, preferredColumns)
    values = strings(0, 1);
    if isempty(reportFile) || ~isfile(reportFile)
        return;
    end
    try
        report = readtable(reportFile, 'TextType', 'string', ...
            'VariableNamingRule', 'preserve');
        names = report.Properties.VariableNames;
        for index = 1:numel(preferredColumns)
            columnIndex = find(strcmp(names, preferredColumns{index}), 1);
            if isempty(columnIndex)
                continue;
            end
            columnValues = strip(string(report.(names{columnIndex})));
            columnValues(ismissing(columnValues)) = "";
            columnValues = columnValues(columnValues ~= "");
            values = [values; columnValues(:)]; %#ok<AGROW>
        end
        values = unique(values, 'stable');
        if numel(values) > 12
            values = [values(1:12); "...TRUNCATED"];
        end
    catch reportError
        values = "REPORT_READ_FAILED:" + compactText(reportError.message);
    end
end


function value = compactText(value)
    value = char(string(value));
    value = regexprep(strtrim(value), '\s+', ' ');
    if numel(value) > 500
        value = [value(1:500) '...'];
    end
end


function [attempted, repairedCount, statuses, reasons] = readGnuSoftRepairSummary(reportFile)
    attempted = false;
    repairedCount = 0;
    statuses = 'UNAVAILABLE';
    reasons = '';
    if ~isfile(reportFile)
        return;
    end

    report = readtable(reportFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    names = report.Properties.VariableNames;
    statusIndex = find(strcmp(names, 'soft_repair_status'), 1);
    repairedIndex = find(strcmp(names, 'soft_repaired_count'), 1);
    reasonIndex = find(strcmp(names, 'soft_repair_reason'), 1);
    if isempty(statusIndex)
        statuses = 'LEGACY_REPORT_NO_SOFT_REPAIR_COLUMNS';
        return;
    end

    statusValues = string(report.(names{statusIndex}));
    statusValues(ismissing(statusValues)) = "";
    statusValues = strip(statusValues);
    statusValues = statusValues(statusValues ~= "");
    if isempty(statusValues)
        statuses = 'NOT_REPORTED';
    else
        statuses = char(strjoin(unique(statusValues, 'stable'), ';'));
        attempted = any(~ismember(statusValues, ["NOT_RUN", "NOT_NEEDED"]));
    end

    if ~isempty(repairedIndex)
        values = report.(names{repairedIndex});
        if isnumeric(values)
            repairedCount = sum(values, 'omitnan');
        else
            repairedCount = sum(str2double(string(values)), 'omitnan');
        end
    end
    if ~isempty(reasonIndex)
        reasonValues = string(report.(names{reasonIndex}));
        reasonValues(ismissing(reasonValues)) = "";
        reasonValues = strip(reasonValues);
        reasonValues = reasonValues(reasonValues ~= "");
        if ~isempty(reasonValues)
            reasons = char(strjoin(unique(reasonValues, 'stable'), ' | '));
        end
    end
end


function quoted = quoteCommandArg(value)
    value = char(value);
    if contains(value, '"')
        error('SCIONX:PipelineQuote', ...
            'Paths containing a double quote are not supported: %s', value);
    end
    quoted = ['"' value '"'];
end
