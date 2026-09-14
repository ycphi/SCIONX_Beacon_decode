function frames = decode_gnuradio_hex(inputRoot, outputRoot, reportCsv)
%DECODE_GNURADIO_HEX Deframe raw SCION-X HDLC candidates from GNU Radio.
%
% GNU Radio/Python supplies candidate raw serial bits only. Each HEX byte is
% packed in chronological LSB-first order and still contains HDLC flags,
% stuffed zeroes, and the two FCS bytes. MATLAB owns both strict and
% diagnostic deframing:
%   1. exact flag-to-flag HDLC deframe;
%   2. if strict deframe cannot produce 274 bytes, fixed-length boundary
%      recovery from the last flag in a >=4 flag preamble;
%   3. CRC-16/X-25 validation and Raw_bit.txt output.
%
% A recovered 274-byte frame is retained even when CRC fails so telemetry can
% be decoded for diagnosis. Matlab_HDLC_Status.txt records its trust level.

    demoDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(inputRoot)
        inputRoot = fullfile(demoDir, 'Frames');
    end
    if nargin < 2 || isempty(outputRoot)
        outputRoot = inputRoot;
    end
    if nargin < 3 || isempty(reportCsv)
        reportCsv = fullfile(demoDir, 'Matlab_Frame_Report.csv');
    end

    if isfolder(inputRoot)
        listing = dir(fullfile(inputRoot, 'Frame_*', 'GNU_GFSK_Raw_HEX.txt'));
        inputFiles = arrayfun(@(item) fullfile(item.folder, item.name), ...
            listing, 'UniformOutput', false);
        inputFiles = sort(inputFiles);
    elseif isfile(inputRoot)
        inputFiles = {inputRoot};
    else
        error('SCIONX:InputMissing', 'Input path does not exist: %s', inputRoot);
    end
    if isempty(inputFiles)
        error('SCIONX:NoFrames', ...
            'No Frame_NNN/GNU_GFSK_Raw_HEX.txt files found under: %s', inputRoot);
    end
    if ~isfolder(outputRoot)
        mkdir(outputRoot);
    end

    flag = uint8([0 1 1 1 1 1 1 0]);
    expectedHeader = uint8(hex2dec(strsplit( ...
        '84 9C 60 86 AA 40 60 84 9C 60 A6 86 B0 E1 03 F0'))).';
    expectedHeaderBits = bytesToLsbBits(expectedHeader);
    frameBitCount = 274 * 8;
    frames = {};
    crcPassWritten = 0;
    crcFailWritten = 0;
    records = emptyRecords();

    for inputIndex = 1:numel(inputFiles)
        inputTxt = inputFiles{inputIndex};
        [inputFrameDir, ~, ~] = fileparts(inputTxt);
        [~, frameDirName] = fileparts(inputFrameDir);
        if ~startsWith(frameDirName, 'Frame_')
            frameDirName = sprintf('Frame_%03d', inputIndex);
        end
        outputFrameDir = fullfile(outputRoot, frameDirName);
        if ~isfolder(outputFrameDir)
            mkdir(outputFrameDir);
        end
        outputTxt = fullfile(outputFrameDir, 'Raw_bit.txt');
        matlabStatusTxt = fullfile(outputFrameDir, 'Matlab_HDLC_Status.txt');
        deleteIfPresent(outputTxt);
        deleteIfPresent(matlabStatusTxt);

        textLines = readTextLines(inputTxt);
        softRepair = readSoftRepairMetadata(inputFrameDir);
        writtenForInput = false;
        bestFallback = struct([]);

        % First pass: standards-compliant exact flag-to-flag deframing.
        for lineIndex = 1:numel(textLines)
            [serialBits, hasHex] = parsePackedBits(textLines{lineIndex});
            if ~hasHex
                continue;
            end
            flagPositions = strfind(serialBits, flag);
            for flagIndex = 1:max(0, numel(flagPositions) - 1)
                opening = flagPositions(flagIndex);
                closing = flagPositions(flagIndex + 1);
                if closing <= opening + numel(flag)
                    continue;
                end
                stuffed = serialBits(opening + numel(flag):closing - 1);
                [destuffed, destuffOk, destuffReason] = removeHdlcStuffing(stuffed);
                record = makeRecord(inputTxt, lineIndex, opening, closing, ...
                    numel(stuffed), numel(destuffed), 'STRICT');
                record = applySoftRepairMetadata(record, softRepair);
                record.status = destuffReason;

                if ~destuffOk
                    records(end + 1) = record; %#ok<AGROW>
                    continue;
                end
                if mod(numel(destuffed), 8) ~= 0
                    record.status = 'DESTUFFED_LENGTH_NOT_MULTIPLE_OF_8';
                    records(end + 1) = record; %#ok<AGROW>
                    continue;
                end
                frameBytes = lsbBitsToBytes(destuffed);
                record.frame_bytes = numel(frameBytes);
                if numel(frameBytes) <= 2
                    record.status = 'FRAME_TOO_SHORT';
                    records(end + 1) = record; %#ok<AGROW>
                    continue;
                end
                if numel(frameBytes) ~= 274
                    record.status = 'FRAME_LENGTH_NOT_274';
                    records(end + 1) = record; %#ok<AGROW>
                    continue;
                end

                [record, crcText] = validateRecord(record, frameBytes);
                if writtenForInput
                    record.status = 'COMPLETE_ADDITIONAL_NOT_WRITTEN';
                    records(end + 1) = record; %#ok<AGROW>
                    continue;
                end
                [record, frames, crcPassWritten, crcFailWritten] = ...
                    publishFrame(record, frameBytes, outputTxt, matlabStatusTxt, ...
                    frames, crcPassWritten, crcFailWritten);
                writtenForInput = true;
                records(end + 1) = record; %#ok<AGROW>
                fprintf(['SCION-X frame %d: %s, STRICT, bytes=%d, FCS=%s, ' ...
                    'CRC Check %s\n'], numel(frames), frameDirName, ...
                    numel(frameBytes), record.received_fcs, crcText);
            end
        end

        if writtenForInput
            continue;
        end

        % Second pass: recover a complete 2192-bit frame from a reliable
        % preamble even when an in-frame false flag or damaged closing flag
        % prevents exact flag-to-flag deframing.
        for lineIndex = 1:numel(textLines)
            [serialBits, hasHex] = parsePackedBits(textLines{lineIndex});
            if ~hasHex
                continue;
            end
            openings = preambleOpeningFlags(serialBits, flag, 4);
            for openingIndex = 1:numel(openings)
                opening = openings(openingIndex);
                firstDataBit = opening + numel(flag);
                [destuffed, closing, violations, complete] = ...
                    destuffFixedLength(serialBits, firstDataBit, frameBitCount);
                if ~complete
                    record = makeRecord(inputTxt, lineIndex, opening, NaN, ...
                        max(0, numel(serialBits) - firstDataBit + 1), ...
                        numel(destuffed), 'UNRECOVERED');
                    record = applySoftRepairMetadata(record, softRepair);
                    record.destuff_violations = violations;
                    record.status = 'RECOVERY_INSUFFICIENT_DESTUFFED_BITS';
                    records(end + 1) = record; %#ok<AGROW>
                    continue;
                end

                frameBytes = lsbBitsToBytes(destuffed);
                headerErrors = sum(destuffed(1:numel(expectedHeaderBits)) ~= ...
                    expectedHeaderBits);
                closingErrors = closingFlagErrors(serialBits, closing, flag);
                trustedBoundary = headerErrors <= 16 && closingErrors <= 2 && ...
                    violations <= 8;
                if trustedBoundary
                    boundaryMode = 'MATLAB_RECOVERED';
                else
                    boundaryMode = 'MATLAB_FLAGS_ONLY_DIAGNOSTIC';
                end
                record = makeRecord(inputTxt, lineIndex, opening, closing, ...
                    closing - firstDataBit, numel(destuffed), boundaryMode);
                record = applySoftRepairMetadata(record, softRepair);
                record.frame_bytes = numel(frameBytes);
                record.header_bit_errors = headerErrors;
                record.closing_flag_bit_errors = closingErrors;
                record.destuff_violations = violations;
                [record, ~] = validateRecord(record, frameBytes);
                record.status = fallbackStatus(boundaryMode, record.crc_ok);

                candidate = struct('record', record, 'frame', frameBytes);
                if isempty(bestFallback) || isBetterFallback(candidate, bestFallback)
                    bestFallback = candidate;
                end
            end
        end

        if ~isempty(bestFallback)
            record = bestFallback.record;
            frameBytes = bestFallback.frame;
            [record, frames, crcPassWritten, crcFailWritten] = ...
                publishFrame(record, frameBytes, outputTxt, matlabStatusTxt, ...
                frames, crcPassWritten, crcFailWritten);
            records(end + 1) = record; %#ok<AGROW>
            fprintf(['SCION-X frame %d: %s, %s, header errors=%d, ' ...
                'closing errors=%d, destuff violations=%d, CRC Check %s\n'], ...
                numel(frames), frameDirName, record.boundary_mode, ...
                record.header_bit_errors, record.closing_flag_bit_errors, ...
                record.destuff_violations, passFail(record.crc_ok));
        elseif isempty(textLines)
            record = makeRecord(inputTxt, 0, NaN, NaN, 0, 0, 'UNRECOVERED');
            record = applySoftRepairMetadata(record, softRepair);
            record.status = 'EMPTY_RAW_HEX_FILE';
            records(end + 1) = record; %#ok<AGROW>
        elseif ~any(cellfun(@hasHexText, textLines))
            record = makeRecord(inputTxt, 0, NaN, NaN, 0, 0, 'UNRECOVERED');
            record = applySoftRepairMetadata(record, softRepair);
            record.status = 'NO_HEX_BYTES';
            records(end + 1) = record; %#ok<AGROW>
        elseif ~hasPreambleInLines(textLines, flag)
            record = makeRecord(inputTxt, 0, NaN, NaN, 0, 0, 'UNRECOVERED');
            record = applySoftRepairMetadata(record, softRepair);
            record.status = 'NO_RELIABLE_HDLC_PREAMBLE';
            records(end + 1) = record; %#ok<AGROW>
        end
    end

    if isempty(records)
        reportTable = table();
    else
        reportTable = struct2table(records, 'AsArray', true);
    end
    writetable(reportTable, reportCsv);
    fprintf(['MATLAB decode complete: %d frame file(s), CRC PASS=%d, ' ...
        'CRC FAIL continued=%d -> %s\n'], ...
        numel(frames), crcPassWritten, crcFailWritten, outputRoot);
end


function records = emptyRecords()
    records = struct( ...
        'input_file', {}, 'input_line', {}, 'opening_flag_bit', {}, ...
        'closing_flag_bit', {}, 'stuffed_bit_count', {}, ...
        'destuffed_bit_count', {}, 'frame_bytes', {}, ...
        'received_fcs', {}, 'calculated_fcs', {}, 'crc_ok', {}, ...
        'boundary_mode', {}, 'header_bit_errors', {}, ...
        'closing_flag_bit_errors', {}, 'destuff_violations', {}, ...
        'soft_repair_applied', {}, 'soft_repair_status', {}, ...
        'soft_repair_bit_flip_count', {}, 'soft_repair_method', {}, ...
        'crc_pass_after_soft_repair', {}, ...
        'written_to_raw_bit', {}, 'output_file', {}, 'status', {});
end


function record = makeRecord(inputFile, inputLine, opening, closing, ...
        stuffedCount, destuffedCount, boundaryMode)
    record = struct( ...
        'input_file', inputFile, 'input_line', inputLine, ...
        'opening_flag_bit', opening, 'closing_flag_bit', closing, ...
        'stuffed_bit_count', stuffedCount, ...
        'destuffed_bit_count', destuffedCount, 'frame_bytes', 0, ...
        'received_fcs', '', 'calculated_fcs', '', 'crc_ok', false, ...
        'boundary_mode', boundaryMode, 'header_bit_errors', NaN, ...
        'closing_flag_bit_errors', NaN, 'destuff_violations', 0, ...
        'soft_repair_applied', false, 'soft_repair_status', '', ...
        'soft_repair_bit_flip_count', 0, 'soft_repair_method', '', ...
        'crc_pass_after_soft_repair', false, ...
        'written_to_raw_bit', false, 'output_file', '', 'status', '');
end


function [record, crcText] = validateRecord(record, frameBytes)
    calculated = crc16X25(frameBytes(1:end - 2));
    received = bitor(uint16(frameBytes(end - 1)), ...
        bitshift(uint16(frameBytes(end)), 8));
    record.received_fcs = sprintf('%04X', received);
    record.calculated_fcs = sprintf('%04X', calculated);
    record.crc_ok = calculated == received;
    record.crc_pass_after_soft_repair = ...
        record.soft_repair_applied && record.crc_ok;
    crcText = passFail(record.crc_ok);
end


function [record, frames, passCount, failCount] = publishFrame( ...
        record, frameBytes, outputTxt, statusTxt, frames, passCount, failCount)
    fileId = fopen(outputTxt, 'wt');
    if fileId < 0
        error('SCIONX:OutputOpen', 'Cannot create output file: %s', outputTxt);
    end
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '%s\n', strtrim(sprintf('%02X ', frameBytes)));
    clear cleanup;
    record.written_to_raw_bit = true;
    record.output_file = outputTxt;
    if strcmp(record.boundary_mode, 'STRICT')
        if record.crc_ok
            record.status = 'CRC_VALID_WRITTEN';
        else
            record.status = 'CRC_FAILED_WRITTEN_CONTINUED';
        end
    end
    if record.crc_ok
        passCount = passCount + 1;
    else
        failCount = failCount + 1;
    end
    frames{end + 1} = frameBytes;
    writeMatlabStatus(statusTxt, record);
end


function writeMatlabStatus(filename, record)
    fileId = fopen(filename, 'wt');
    if fileId < 0
        error('SCIONX:StatusOutputOpen', 'Cannot create status file: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, 'Boundary Mode: %s\n', record.boundary_mode);
    fprintf(fileId, 'Opening Flag Bit Index: %g\n', record.opening_flag_bit);
    fprintf(fileId, 'Expected Closing Flag Bit Index: %g\n', record.closing_flag_bit);
    fprintf(fileId, 'Header Bit Errors: %g\n', record.header_bit_errors);
    fprintf(fileId, 'Closing Flag Bit Errors: %g\n', record.closing_flag_bit_errors);
    fprintf(fileId, 'Destuff Violations: %g\n', record.destuff_violations);
    fprintf(fileId, 'CRC Check: %s\n', passFail(record.crc_ok));
    fprintf(fileId, 'Soft repair applied: %s\n', ...
        yesNo(record.soft_repair_applied));
    fprintf(fileId, 'Soft repair status: %s\n', record.soft_repair_status);
    fprintf(fileId, 'Soft repair bit flips: %d\n', ...
        record.soft_repair_bit_flip_count);
    fprintf(fileId, 'Soft repair method: %s\n', record.soft_repair_method);
    fprintf(fileId, 'CRC PASS after soft repair: %s\n', ...
        yesNo(record.crc_pass_after_soft_repair));
    if record.crc_ok
        fprintf(fileId, 'Frame Validation: PASS\n');
        fprintf(fileId, 'Telemetry Use: NORMAL\n');
    else
        fprintf(fileId, 'Frame Validation: FAIL - DIAGNOSTIC ONLY\n');
        fprintf(fileId, 'Telemetry Use: REFERENCE ONLY - NOT VALIDATED\n');
    end
end


function status = fallbackStatus(boundaryMode, crcOk)
    if strcmp(boundaryMode, 'MATLAB_RECOVERED')
        if crcOk
            status = 'MATLAB_RECOVERED_CRC_VALID_WRITTEN';
        else
            status = 'MATLAB_RECOVERED_CRC_FAILED_WRITTEN_CONTINUED';
        end
    elseif crcOk
        status = 'MATLAB_FLAGS_ONLY_DIAGNOSTIC_CRC_MATCH_UNTRUSTED';
    else
        status = 'MATLAB_FLAGS_ONLY_DIAGNOSTIC_CRC_FAILED_WRITTEN';
    end
end


function better = isBetterFallback(candidate, current)
    a = fallbackKey(candidate.record);
    b = fallbackKey(current.record);
    difference = find(a ~= b, 1);
    better = ~isempty(difference) && a(difference) > b(difference);
end


function key = fallbackKey(record)
    trusted = strcmp(record.boundary_mode, 'MATLAB_RECOVERED');
    key = [double(trusted), double(record.crc_ok), ...
        -double(record.header_bit_errors), ...
        -double(record.closing_flag_bit_errors), ...
        -double(record.destuff_violations)];
end


function [bits, hasHex] = parsePackedBits(line)
    tokens = regexp(line, ...
        '(?i)(?<![0-9A-F])[0-9A-F]{2}(?![0-9A-F])', 'match');
    hasHex = ~isempty(tokens);
    if hasHex
        bits = bytesToLsbBits(uint8(hex2dec(tokens)));
    else
        bits = uint8([]);
    end
end


function lines = readTextLines(filename)
    fileId = fopen(filename, 'rt');
    if fileId < 0
        error('SCIONX:InputOpen', 'Cannot open GNU Radio HEX file: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fileId));
    text = textscan(fileId, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = text{1};
end


function openings = preambleOpeningFlags(bits, flag, minimumFlags)
    positions = strfind(bits, flag);
    openings = [];
    if isempty(positions)
        return;
    end
    runStart = 1;
    for index = 2:(numel(positions) + 1)
        runEnded = index > numel(positions) || ...
            positions(index) - positions(index - 1) ~= numel(flag);
        if runEnded
            if index - runStart >= minimumFlags
                openings(end + 1) = positions(index - 1); %#ok<AGROW>
            end
            runStart = index;
        end
    end
end


function [output, rawIndex, violations, complete] = ...
        destuffFixedLength(input, firstDataBit, targetCount)
    input = uint8(input(:).');
    output = zeros(1, targetCount, 'uint8');
    outputCount = 0;
    onesCount = 0;
    rawIndex = firstDataBit;
    violations = 0;
    while outputCount < targetCount && rawIndex <= numel(input)
        bit = input(rawIndex);
        rawIndex = rawIndex + 1;
        outputCount = outputCount + 1;
        output(outputCount) = bit;
        if bit == 1
            onesCount = onesCount + 1;
        else
            onesCount = 0;
        end
        if onesCount == 5
            if rawIndex > numel(input)
                break;
            end
            if input(rawIndex) ~= 0
                violations = violations + 1;
            end
            rawIndex = rawIndex + 1;
            onesCount = 0;
        end
    end
    complete = outputCount == targetCount;
    output = output(1:outputCount);
end


function errors = closingFlagErrors(bits, closing, flag)
    if closing > numel(bits)
        errors = numel(flag);
        return;
    end
    available = bits(closing:min(closing + numel(flag) - 1, numel(bits)));
    errors = numel(flag) - numel(available) + ...
        sum(available ~= flag(1:numel(available)));
end


function present = hasPreambleInLines(lines, flag)
    present = false;
    for index = 1:numel(lines)
        [bits, hasHex] = parsePackedBits(lines{index});
        if hasHex && ~isempty(preambleOpeningFlags(bits, flag, 4))
            present = true;
            return;
        end
    end
end


function present = hasHexText(line)
    present = ~isempty(regexp(line, ...
        '(?i)(?<![0-9A-F])[0-9A-F]{2}(?![0-9A-F])', 'once'));
end


function deleteIfPresent(filename)
    if isfile(filename)
        delete(filename);
    end
end


function text = passFail(value)
    if value
        text = 'PASS';
    else
        text = 'FAIL';
    end
end


function text = yesNo(value)
    if value
        text = 'YES';
    else
        text = 'NO';
    end
end


function metadata = readSoftRepairMetadata(frameDir)
    metadata = struct('applied', false, 'status', '', ...
        'bit_flip_count', 0, 'method', '');
    filename = fullfile(frameDir, 'Soft_Repair_Status.json');
    if ~isfile(filename)
        return;
    end
    try
        document = jsondecode(fileread(filename));
        if isfield(document, 'soft_repair_status')
            metadata.status = char(string(document.soft_repair_status));
        end
        if isfield(document, 'repair_bit_flip_count')
            metadata.bit_flip_count = double(document.repair_bit_flip_count);
        elseif isfield(document, 'algorithm_result') && ...
                isstruct(document.algorithm_result) && ...
                isfield(document.algorithm_result, 'global_flip_indices')
            metadata.bit_flip_count = ...
                numel(document.algorithm_result.global_flip_indices);
        end
        metadata.applied = metadata.bit_flip_count > 0;
        if isfield(document, 'repair_applied')
            metadata.applied = logical(document.repair_applied) && ...
                metadata.bit_flip_count > 0;
        end
        if isfield(document, 'algorithm_result') && ...
                isstruct(document.algorithm_result) && ...
                isfield(document.algorithm_result, 'repair_method')
            metadata.method = char(string( ...
                document.algorithm_result.repair_method));
        end
    catch metadataError
        metadata.status = ['METADATA_READ_FAILED: ' metadataError.message];
    end
end


function record = applySoftRepairMetadata(record, metadata)
    record.soft_repair_applied = metadata.applied;
    record.soft_repair_status = metadata.status;
    record.soft_repair_bit_flip_count = metadata.bit_flip_count;
    record.soft_repair_method = metadata.method;
end


function bits = bytesToLsbBits(bytes)
    bytes = uint8(bytes(:));
    bitMatrix = zeros(numel(bytes), 8, 'uint8');
    for bitPosition = 1:8
        bitMatrix(:, bitPosition) = uint8(bitget(bytes, bitPosition));
    end
    bits = reshape(bitMatrix.', 1, []);
end


function bytes = lsbBitsToBytes(bits)
    bits = uint8(bits(:).');
    if isempty(bits)
        bytes = uint8([]);
        return;
    end
    bitMatrix = reshape(bits, 8, []).';
    weights = double(2 .^ (0:7)).';
    bytes = uint8(double(bitMatrix) * weights).';
end


function [output, valid, reason] = removeHdlcStuffing(input)
    input = uint8(input(:).');
    output = zeros(1, numel(input), 'uint8');
    outputCount = 0;
    onesCount = 0;
    index = 1;
    valid = true;
    reason = 'DESTUFF_OK';
    while index <= numel(input)
        bit = input(index);
        outputCount = outputCount + 1;
        output(outputCount) = bit;
        if bit == 1
            onesCount = onesCount + 1;
            if onesCount > 6
                valid = false;
                reason = 'HDLC_ABORT_OR_INVALID_ONES_RUN';
                output = output(1:outputCount);
                return;
            end
        else
            onesCount = 0;
        end
        index = index + 1;
        if onesCount == 5
            if index > numel(input) || input(index) ~= 0
                valid = false;
                reason = 'MISSING_STUFFED_ZERO';
                output = output(1:outputCount);
                return;
            end
            index = index + 1;
            onesCount = 0;
        end
    end
    output = output(1:outputCount);
end


function crc = crc16X25(bytes)
    crc = uint16(hex2dec('FFFF'));
    polynomial = uint16(hex2dec('8408'));
    bytes = uint8(bytes(:).');
    for value = bytes
        crc = bitxor(crc, uint16(value));
        for bitIndex = 1:8
            if bitand(crc, uint16(1)) ~= 0
                crc = bitxor(bitshift(crc, -1), polynomial);
            else
                crc = bitshift(crc, -1);
            end
        end
    end
    crc = bitcmp(crc);
end
