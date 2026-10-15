function telemetry = decode_scionx_tlm(inputPath)
%DECODE_SCIONX_TLM Convert per-frame Raw_bit.txt files to SCION-X TLM.txt.
%
% telemetry = decode_scionx_tlm(inputPath)
%
% inputPath may be either:
%   * the Frames directory containing Frame_NNN/Raw_bit.txt; or
%   * one individual Raw_bit.txt file.
%
% Raw_bit.txt is ASCII hexadecimal.  Byte positions in this decoder are
% one-based positions in the AX.25 frame, exactly as defined in the SCION-X
% telemetry map.  Each TLM.txt is written next to its Raw_bit.txt.

    if nargin < 1 || isempty(inputPath)
        inputPath = fullfile(fileparts(mfilename('fullpath')), 'Frames');
    end

    if isfolder(inputPath)
        listing = dir(fullfile(inputPath, 'Frame_*', 'Raw_bit.txt'));
        inputFiles = arrayfun(@(item) fullfile(item.folder, item.name), ...
            listing, 'UniformOutput', false);
        inputFiles = sort(inputFiles);
    elseif isfile(inputPath)
        inputFiles = {inputPath};
    else
        error('SCIONX:TlmInputMissing', 'Input path does not exist: %s', inputPath);
    end
    if isempty(inputFiles)
        error('SCIONX:TlmNoFrames', ...
            'No Frame_NNN/Raw_bit.txt files found under: %s', inputPath);
    end

    telemetry = repmat(struct('input_file', '', 'output_file', '', ...
        'frame_bytes', 0, 'destination', '', 'source', '', 'apid', 0, ...
        'obc_unix_time', uint32(0), 'boundary_mode', 'STRICT'), 1, numel(inputFiles));

    for frameIndex = 1:numel(inputFiles)
        rawFile = inputFiles{frameIndex};
        frameBytes = readHexBytes(rawFile);
        if numel(frameBytes) < 271
            error('SCIONX:TlmFrameTooShort', ...
                'Frame has %d bytes; at least 271 are required: %s', ...
                numel(frameBytes), rawFile);
        end

        outputFile = fullfile(fileparts(rawFile), 'TLM.txt');
        writeTelemetryFile(outputFile, rawFile, frameBytes);

        [destination, ~, ~] = decodeAx25Address(frameBytes(1:7));
        [source, ~, ~] = decodeAx25Address(frameBytes(8:14));
        telemetry(frameIndex).input_file = rawFile;
        telemetry(frameIndex).output_file = outputFile;
        telemetry(frameIndex).frame_bytes = numel(frameBytes);
        telemetry(frameIndex).destination = destination;
        telemetry(frameIndex).source = source;
        telemetry(frameIndex).apid = double(frameBytes(17));
        telemetry(frameIndex).obc_unix_time = leU32(frameBytes(20:23));
        telemetry(frameIndex).boundary_mode = readBoundaryMode(fileparts(rawFile));

        fprintf('SCION-X TLM %d: %s -> %s\n', ...
            frameIndex, rawFile, outputFile);
    end
    fprintf('SCION-X telemetry decode complete: %d TLM.txt file(s)\n', ...
        numel(inputFiles));
end


function bytes = readHexBytes(filename)
    text = fileread(filename);
    tokens = regexp(text, '(?i)(?<![0-9A-F])[0-9A-F]{2}(?![0-9A-F])', 'match');
    if isempty(tokens)
        error('SCIONX:TlmNoHex', 'No hexadecimal bytes found in: %s', filename);
    end
    bytes = uint8(hex2dec(tokens(:))).';
end


function writeTelemetryFile(outputFile, inputFile, b)
    fileId = fopen(outputFile, 'wt', 'n', 'UTF-8');
    if fileId < 0
        error('SCIONX:TlmOutputOpen', 'Cannot create output file: %s', outputFile);
    end
    cleanup = onCleanup(@() fclose(fileId));

    [destination, destinationCall, destinationSsid] = decodeAx25Address(b(1:7));
    [source, sourceCall, sourceSsid] = decodeAx25Address(b(8:14));
    control = b(15);
    pid = b(16);
    apid = b(17);
    obcTime = leU32(b(20:23));
    receivedFcs = bitor(uint16(b(end - 1)), bitshift(uint16(b(end)), 8));
    calculatedFcs = crc16X25(b(1:end - 2));
    if receivedFcs == calculatedFcs
        fcsStatus = 'PASS';
    else
        fcsStatus = 'FAIL';
    end
    boundaryMode = readBoundaryMode(fileparts(inputFile));
    if strcmp(fcsStatus, 'FAIL') || contains(boundaryMode, 'DIAGNOSTIC')
        frameValidation = 'FAIL - DIAGNOSTIC ONLY';
        telemetryUse = 'REFERENCE ONLY - NOT VALIDATED';
    else
        frameValidation = 'PASS';
        telemetryUse = 'NORMAL';
    end

    fprintf(fileId, 'SCION-X Telemetry Decode\n');
    fprintf(fileId, 'Input file = %s\n', inputFile);
    fprintf(fileId, 'Frame length = %d bytes\n', numel(b));
    fprintf(fileId, 'Boundary Mode = %s\n', boundaryMode);
    fprintf(fileId, 'CRC Check = %s\n', fcsStatus);
    fprintf(fileId, 'Frame Validation = %s\n', frameValidation);
    fprintf(fileId, 'Telemetry Use = %s\n', telemetryUse);
    fprintf(fileId, 'Byte numbering = 1-based AX.25 frame positions\n');
    fprintf(fileId, 'Integer byte order = little-endian unless noted\n\n');

    fprintf(fileId, '[AX.25 / HDLC]\n');
    fprintf(fileId, 'Bytes 1-7   Destination = %s (call=%s, SSID=%d)\n', ...
        destination, destinationCall, destinationSsid);
    fprintf(fileId, 'Bytes 8-14  Source      = %s (call=%s, SSID=%d)\n', ...
        source, sourceCall, sourceSsid);
    fprintf(fileId, 'Byte 15     Control     = 0x%02X (%s)\n', ...
        control, controlName(control));
    fprintf(fileId, 'Byte 16     PID         = 0x%02X (%s)\n\n', ...
        pid, pidName(pid));

    fprintf(fileId, '[SCION-X Packet Header]\n');
    fprintf(fileId, 'Byte 17     APID = 0x%02X (%s)\n', apid, packetApidName(apid));
    writeBitGroup(fileId, 18, 'LEOP_and_deployment', b(18), ...
        {'is_first_boot', 'is_going_to_adcs_commission', 'is_adcs_commissioned', ...
         'uhf', 'vhf', 'sp2', 'sp1', 'mag'}, ...
        repmat({'false'}, 1, 8), repmat({'true'}, 1, 8));
    writeBitGroup(fileId, 19, 'system_state_main_SD_nominal', b(19), ...
        {'is_going_to_factory_reset', 'is_going_to_reset', 'main_SD', ...
         'system_state', 'uhf_tx_power', 'reserved_5', 'reserved_6', 'reserved_7'}, ...
        {'false', 'false', 'sd_0', 'contingency (safe/phx)', 'low', '0', '0', '0'}, ...
        {'true', 'true', 'sd_1', 'nominal (charging/COMM/SCI)', 'high', '1', '1', '1'});
    fprintf(fileId, 'Bytes 20-23 OBC_time = %u (0x%08X, %s)\n\n', ...
        obcTime, obcTime, unixTimeText(obcTime));

    fprintf(fileId, '[FSW]\n');
    fprintf(fileId, 'Byte 24     Current_Mode = 0x%02X (%s)\n', b(24), fswModeName(b(24)));
    writeBitGroup(fileId, 25, 'Subsystems_on_off_state', b(25), ...
        {'EPS', 'ADCS', 'UHF', 'XBand', 'CIP', 'HyperSCAN', 'SEUV_Probe', 'APRS'}, ...
        repmat({'off'}, 1, 8), repmat({'on'}, 1, 8));
    writeU16(fileId, 'Bytes 26-27 FSW_Reboot_Count', b(26:27), 'count');
    writeU16(fileId, 'Bytes 28-29 CMD_Accept_Count', b(28:29), 'count');
    writeU16(fileId, 'Bytes 30-31 CMD_Reject_Count', b(30:31), 'count');
    writeU32(fileId, 'Bytes 32-35 last_packet_send', b(32:35), 'raw uint32');
    writeU32(fileId, 'Bytes 36-39 last_packet_send_time', b(36:39), 'raw uint32');
    writeU32(fileId, 'Bytes 40-43 last_Accept_cmd_time', b(40:43), 'raw uint32');
    fprintf(fileId, 'Byte 44     last_Accept_CMD_Type = 0x%02X, APID=0x%X, OPCODE=0x%X\n', ...
        b(44), bitshift(b(44), -4), bitand(b(44), uint8(15)));
    writeU32(fileId, 'Bytes 45-48 FSW_static_data_address', b(45:48), 'address');
    writeU32(fileId, 'Bytes 49-52 FSW_dynamic_data_address', b(49:52), 'address');
    writeU32(fileId, 'Bytes 53-56 SD_Write_Pointer_Beacon', b(53:56), 'raw uint32');
    writeU32(fileId, 'Bytes 57-60 SD_Write_Pointer_plCIP', b(57:60), 'raw uint32');
    writeU32(fileId, 'Bytes 61-64 SD_Write_Pointer_plSEUV', b(61:64), 'raw uint32');
    writeU32(fileId, 'Bytes 65-68 SD_Write_Pointer_plHyperSCAN', b(65:68), 'raw uint32');
    writeU32(fileId, 'Bytes 69-72 SD_Write_Pointer_AdcsDetailData', b(69:72), 'raw uint32');
    writeU32(fileId, 'Bytes 73-76 CMD_Loss_Timer', b(73:76), 'seconds');
    writeU16(fileId, 'Bytes 77-78 HyperSCAN_Xband_Write_file_number', b(77:78), 'file number');
    writeU16(fileId, 'Bytes 79-80 SEUV_Xband_Write_file_number', b(79:80), 'file number');
    writeU16(fileId, 'Bytes 81-82 CIP_Xband_Write_file_number', b(81:82), 'file number');
    fprintf(fileId, '\n');

    fprintf(fileId, '[ADCS]\n');
    adcsTime = leU32(b(83:86));
    fprintf(fileId, 'Bytes 83-86 ADCS_time = %u (0x%08X, %s)\n', ...
        adcsTime, adcsTime, unixTimeText(adcsTime));
    writeU16(fileId, 'Bytes 87-88 ADCS_time_msec', b(87:88), 'milliseconds');
    writeU16(fileId, 'Bytes 89-90 ADCS_3V3_I', b(89:90), 'raw uint16; scale not supplied');
    writeU16(fileId, 'Bytes 91-92 ADCS_5V_I', b(91:92), 'raw uint16; scale not supplied');
    writeU16(fileId, 'Bytes 93-94 ADCS_Vbat_I', b(93:94), 'raw uint16; scale not supplied');
    writeI16(fileId, 'Bytes 95-96 ADCS_MCU_T', b(95:96), 'raw int16; scale not supplied');
    writeI16(fileId, 'Bytes 97-98 ADCS_MAG_T', b(97:98), 'raw int16; scale not supplied');

    estimateMode = bitand(b(99), uint8(15));
    controlMode = bitshift(b(99), -4);
    fprintf(fileId, 'Byte 99     ADCS_AE_Control_Mode = 0x%02X, bits[7:0]=%s\n', ...
        b(99), bitsText(b(99)));
    fprintf(fileId, '  bits[3:0] Attitude_Estimate_Mode = 0x%X (%s)\n', ...
        estimateMode, attitudeEstimateName(estimateMode));
    fprintf(fileId, '  bits[7:4] Control_Mode = 0x%X (%s)\n', ...
        controlMode, adcsControlName(controlMode));

    runMode = bitand(b(100), uint8(3));
    fprintf(fileId, 'Byte 100    ADCS_Run_Mode = 0x%02X, bits[7:0]=%s\n', ...
        b(100), bitsText(b(100)));
    fprintf(fileId, '  bits[1:0] Run_Mode = %u (%s)\n', runMode, adcsRunModeName(runMode));
    fprintf(fileId, '  bits[7:2] Reserved = 0x%02X\n', bitshift(b(100), -2));

    writeBitGroup(fileId, 101, 'ADCS_Sensor_Actuator_state', b(101), ...
        {'CubeControl_Signal', 'CubeControl_Motor', 'CubeSense1', 'CubeSense2', ...
         'CubeWheel1', 'CubeWheel2', 'CubeWheel3', 'Cube_Star'}, ...
        repmat({'off'}, 1, 8), repmat({'enable'}, 1, 8));
    writeBitGroup(fileId, 102, 'ADCS_Error_Flag0', b(102), ...
        {'CubeSense1_communications', 'CubeSense2_communications', ...
         'CubeControl_Signal_communications', 'CubeControl_Motor_communications', ...
         'CubeWheel1_communications', 'CubeWheel2_communications', ...
         'CubeWheel3_communications', 'Cube_Star_communications'}, ...
        repmat({'no error'}, 1, 8), repmat({'error'}, 1, 8));
    writeBitGroup(fileId, 103, 'ADCS_Error_Flag1', b(103), ...
        {'Cam1_SRAM_Overcurrent_Detected', 'Cam1_3V3_Overcurrent_Detected', ...
         'Cam1_Sensor_Busy', 'Cam1_Sensor_Detection', ...
         'Cam2_SRAM_Overcurrent_Detected', 'Cam2_3V3_Overcurrent_Detected', ...
         'Cam2_Sensor_Busy', 'Cam2_Sensor_Detection'}, ...
        repmat({'no error'}, 1, 8), repmat({'error'}, 1, 8));
    writeBitGroup(fileId, 104, 'ADCS_Error_Flag2', b(104), ...
        {'Magnetometer_Range', 'Sun_Sensor_Range', 'Nadir_Sensor_Range', ...
         'Rate_Sensor_Range', 'Wheel_Speed_Range', 'Coarse_Sun_Sensor', ...
         'StarTracker_Match', 'Star_Tracker_Overcurrent_Detected'}, ...
        repmat({'no error'}, 1, 8), repmat({'error'}, 1, 8));

    fprintf(fileId, '\n[ADCS Navigation / Commissioning / Telecommand]\n');
    writeScaledI16(fileId, 'Bytes 105-106 ADCS_ECI_Ref_X', b(105:106), 0.25, 'km');
    writeScaledI16(fileId, 'Bytes 107-108 ADCS_ECI_Ref_Y', b(107:108), 0.25, 'km');
    writeScaledI16(fileId, 'Bytes 109-110 ADCS_ECI_Ref_Z', b(109:110), 0.25, 'km');
    writeScaledI16(fileId, 'Bytes 111-112 ADCS_X_Angular_Rate', b(111:112), 0.01, 'deg/s');
    writeScaledI16(fileId, 'Bytes 113-114 ADCS_Y_Angular_Rate', b(113:114), 0.01, 'deg/s');
    writeScaledI16(fileId, 'Bytes 115-116 ADCS_Z_Angular_Rate', b(115:116), 0.01, 'deg/s');
    writeScaledI16(fileId, 'Bytes 117-118 ADCS_Est_X_Angular_Rate', b(117:118), 0.001, 'deg/s');
    writeScaledI16(fileId, 'Bytes 119-120 ADCS_Est_Y_Angular_Rate', b(119:120), 0.001, 'deg/s');
    writeScaledI16(fileId, 'Bytes 121-122 ADCS_Est_Z_Angular_Rate', b(121:122), 0.001, 'deg/s');

    commissioningStep = bitand(b(123), uint8(15));
    commissioningError = bitshift(b(123), -4);
    fprintf(fileId, 'Byte 123    ADCS_Commissioning_State = 0x%02X, bits[7:0]=%s\n', ...
        b(123), bitsText(b(123)));
    fprintf(fileId, '  bits[3:0] Current_steps = 0x%X (%s)\n', ...
        commissioningStep, commissioningStepName(commissioningStep));
    fprintf(fileId, '  bits[7:4] Error_flag = 0x%X (%s)\n', ...
        commissioningError, commissioningErrorName(commissioningError));
    writeScaledI16(fileId, 'Bytes 124-125 ADCS_SunX', b(124:125), 1 / 10000, 'normalized');
    writeScaledI16(fileId, 'Bytes 126-127 ADCS_SunY', b(126:127), 1 / 10000, 'normalized');
    writeScaledI16(fileId, 'Bytes 128-129 ADCS_SunZ', b(128:129), 1 / 10000, 'normalized');
    writeU16(fileId, 'Bytes 130-131 ADCS_TC_Counter', b(130:131), 'count');
    writeU16(fileId, 'Bytes 132-133 ADCS_TLM_Request_Count', b(132:133), 'count');
    fprintf(fileId, 'Byte 134    ADCS_Last_TC_ID = 0x%02X\n', b(134));
    processed = bitget(b(135), 1);
    tcError = bitand(bitshift(b(135), -1), uint8(15));
    fprintf(fileId, 'Byte 135    ADCS_Processed_Flag_and_TC_Error_Status = 0x%02X, bits[7:0]=%s\n', ...
        b(135), bitsText(b(135)));
    fprintf(fileId, '  bit[0] Processed_flag = %d (%s)\n', processed, ...
        chooseText(processed, 'has NOT been processed', 'has been processed'));
    fprintf(fileId, '  bits[4:1] TC_error_status = 0x%X (%s)\n', ...
        tcError, tcErrorName(tcError));
    fprintf(fileId, '  bits[7:5] Reserved = 0x%X\n', bitshift(b(135), -5));
    fprintf(fileId, 'Byte 136    ADCS_TC_Param_Error_Index = %u (0x%02X)\n', b(136), b(136));
    writeI16(fileId, 'Bytes 137-138 ADCS_X_Wheel_Speed', b(137:138), 'rpm');
    writeI16(fileId, 'Bytes 139-140 ADCS_Y_Wheel_Speed', b(139:140), 'rpm');
    writeI16(fileId, 'Bytes 141-142 ADCS_Z_Wheel_Speed', b(141:142), 'rpm');
    fprintf(fileId, 'Byte 143    ADCS_Sun_Above_Horizon = %u (0x%02X)\n', b(143), b(143));

    fprintf(fileId, '\n[UHF]\n');
    writeI8(fileId, 'Byte 144 UHF_Chassis_Temperature', b(144), 'degC');
    fprintf(fileId, 'Byte 145    UHF_Configuration = 0x%02X, bits[7:0]=%s\n', ...
        b(145), bitsText(b(145)));
    fprintf(fileId, '  bit[0] Transceiver_Unlocked_Condition = %d (%s)\n', ...
        bitget(b(145), 1), chooseText(bitget(b(145), 1), 'lock', 'unlock'));
    fprintf(fileId, '  bit[0] SWD_Function = %d (%s) [source map also assigns bit[0]]\n', ...
        bitget(b(145), 1), chooseText(bitget(b(145), 1), 'off', 'on'));
    fprintf(fileId, '  bit[1] Echo_Mode = %d (%s)\n', bitget(b(145), 2), ...
        chooseText(bitget(b(145), 2), 'off', 'on'));
    fprintf(fileId, '  bit[2] AFC = %d (%s)\n', bitget(b(145), 3), ...
        chooseText(bitget(b(145), 3), 'off', 'on'));
    fprintf(fileId, '  bit[3] SWD_Annunciate = %d (%s)\n', bitget(b(145), 4), ...
        chooseText(bitget(b(145), 4), 'off', 'on'));
    uhfReadback = bitand(bitshift(b(145), -4), uint8(3));
    uhfChannel = bitshift(b(145), -6);
    fprintf(fileId, '  bits[5:4] Readback_Mode = %u (%s)\n', ...
        uhfReadback, readbackModeName(uhfReadback));
    fprintf(fileId, '  bits[7:6] UHF_Channel = %u (%s)\n', ...
        uhfChannel, uhfChannelName(uhfChannel));

    fprintf(fileId, '\n[X-band]\n');
    fprintf(fileId, 'Byte 146    Xband_Symbol_Rate = 0x%02X\n', b(146));
    fprintf(fileId, 'Byte 147    Xband_Tx_Rate = 0x%02X\n', b(147));
    fprintf(fileId, 'Byte 148    Xband_MODCOD = 0x%02X\n', b(148));
    xbandRolloff = bitand(b(149), uint8(3));
    xbandReadback = bitand(bitshift(b(149), -4), uint8(3));
    fprintf(fileId, 'Byte 149    Xband_Configuration = 0x%02X, bits[7:0]=%s\n', ...
        b(149), bitsText(b(149)));
    fprintf(fileId, '  bits[1:0] Rolloff = %u\n', xbandRolloff);
    fprintf(fileId, '  bit[2] AFC = %d (%s)\n', bitget(b(149), 3), ...
        chooseText(bitget(b(149), 3), 'off', 'on'));
    fprintf(fileId, '  bit[3] SWD_Annunciate = %d (%s)\n', bitget(b(149), 4), ...
        chooseText(bitget(b(149), 4), 'off', 'on'));
    fprintf(fileId, '  bits[5:4] Readback_Mode = %u (%s)\n', ...
        xbandReadback, readbackModeName(xbandReadback));
    fprintf(fileId, '  bits[7:6] Reserved = 0x%X\n', bitshift(b(149), -6));
    xbandFrequency = leFloat32(b(150:153));
    fprintf(fileId, 'Bytes 150-153 Xband_Central_Freq = %.9g (IEEE-754 float32; raw %s)\n', ...
        xbandFrequency, hexBytes(b(150:153)));
    writeBitGroup(fileId, 154, 'Xband_Status_Flags', b(154), ...
        {'BUC_PLL_Lock_Unlock', 'BUC_Over_Temperature', 'SDR_Not_Initialized', ...
         'BUC_Communication_Error', 'System_Over_Temperature_Stop', ...
         'Warning_Predistortion_Not_OK', 'reserved_6', 'reserved_7'}, ...
        {'0', 'no error', 'initialized', 'no error', 'normal', 'normal', '0', '0'}, ...
        {'1', 'flag set', 'not initialized', 'error', 'stop', 'warning', '1', '1'});
    writeScaledI16(fileId, 'Bytes 155-156 Xband_PA_Temp', b(155:156), 0.01, 'degC');
    writeScaledI16(fileId, 'Bytes 157-158 Xband_Measured_TxPower', b(157:158), 0.01, 'raw scaled power');

    fprintf(fileId, '\n[EPS Solar Panels / Battery]\n');
    writeAdcVoltage(fileId, 'Bytes 159-160 EPS_SP0_V', b(159:160));
    writeAdcVoltage(fileId, 'Bytes 161-162 EPS_SP1_V', b(161:162));
    writeAdcVoltage(fileId, 'Bytes 163-164 EPS_SP2_V', b(163:164));
    writeAdcVoltage(fileId, 'Bytes 165-166 EPS_SP3_V', b(165:166));
    writeAdcVoltage(fileId, 'Bytes 167-168 EPS_SP4_V', b(167:168));
    writeAdcCurrent(fileId, 'Bytes 169-170 EPS_SP0_I', b(169:170));
    writeAdcCurrent(fileId, 'Bytes 171-172 EPS_SP1_I', b(171:172));
    writeAdcCurrent(fileId, 'Bytes 173-174 EPS_SP2_I', b(173:174));
    writeAdcCurrent(fileId, 'Bytes 175-176 EPS_SP3_I', b(175:176));
    writeAdcCurrent(fileId, 'Bytes 177-178 EPS_SP4_I', b(177:178));
    writeAdcVoltage(fileId, 'Bytes 179-180 EPS_Battery_V', b(179:180));
    writeAdcCurrent(fileId, 'Bytes 181-182 EPS_Battery_I', b(181:182));
    writeScaledU16(fileId, 'Bytes 183-184 EPS_mainBattery_V', b(183:184), 0.001, 'V');
    writeScaledI16(fileId, 'Bytes 185-186 EPS_mainBattery_I', b(185:186), 0.001, 'A');
    writeScaledU16(fileId, 'Bytes 187-188 EPS_main_SOC', b(187:188), 0.01, '%');
    writeI8(fileId, 'Byte 189 EPS_Battery_Temp1', b(189), 'degC');
    writeI8(fileId, 'Byte 190 EPS_Battery_Temp2', b(190), 'degC');
    writeI8(fileId, 'Byte 191 EPS_Battery_Temp3', b(191), 'degC');
    writeI8(fileId, 'Byte 192 EPS_SP0_Temp', b(192), 'degC');
    writeI8(fileId, 'Byte 193 EPS_SP1_Temp', b(193), 'degC');
    writeI8(fileId, 'Byte 194 EPS_SP2_Temp', b(194), 'degC');
    writeI8(fileId, 'Byte 195 EPS_SP3_Temp', b(195), 'degC');
    writeI8(fileId, 'Byte 196 EPS_SP4_Temp', b(196), 'degC');
    writeI8(fileId, 'Byte 197 EPS_main_Temp', b(197), 'degC');
    writeI8(fileId, 'Byte 198 EPS_Distribution_Temp', b(198), 'degC');
    writeI8(fileId, 'Byte 199 EPS_CDG3V3_Temp', b(199), 'degC');

    fprintf(fileId, '\n[EPS Power Rails]\n');
    writeAdcVoltage(fileId, 'Bytes 200-201 EPS_CDH3V3', b(200:201));
    writeAdcCurrent(fileId, 'Bytes 202-203 EPS_CDH3V3_I', b(202:203));
    writeAdcVoltage(fileId, 'Bytes 204-205 EPS_UHF7V', b(204:205));
    writeAdcCurrent(fileId, 'Bytes 206-207 EPS_UHF7V_I', b(206:207));
    writeAdcVoltage(fileId, 'Bytes 208-209 EPS_HSCAN5V', b(208:209));
    writeAdcCurrent(fileId, 'Bytes 210-211 EPS_HSCAN5V_I', b(210:211));
    writeAdcVoltage(fileId, 'Bytes 212-213 EPS_ADCS5V', b(212:213));
    writeAdcCurrent(fileId, 'Bytes 214-215 EPS_ADCS5V_I', b(214:215));
    writeAdcVoltage(fileId, 'Bytes 216-217 EPS_ADCS12V', b(216:217));
    writeAdcCurrent(fileId, 'Bytes 218-219 EPS_ADCS12V_I', b(218:219));
    writeAdcVoltage(fileId, 'Bytes 220-221 EPS_ADCS3V3', b(220:221));
    writeAdcCurrent(fileId, 'Bytes 222-223 EPS_ADCS3V3_I', b(222:223));
    writeAdcVoltage(fileId, 'Bytes 224-225 EPS_CIP12V', b(224:225));
    writeAdcCurrent(fileId, 'Bytes 226-227 EPS_CIP12V_I', b(226:227));
    writeAdcVoltage(fileId, 'Bytes 228-229 EPS_Xband12V', b(228:229));
    writeAdcCurrent(fileId, 'Bytes 230-231 EPS_Xband12V_I', b(230:231));
    writeAdcVoltage(fileId, 'Bytes 232-233 EPS_SEUV5V', b(232:233));
    writeAdcCurrent(fileId, 'Bytes 234-235 EPS_SEUV5V_I', b(234:235));
    writeAdcVoltage(fileId, 'Bytes 236-237 EPS_SEUV3V3', b(236:237));
    writeAdcCurrent(fileId, 'Bytes 238-239 EPS_SEUV3V3_I', b(238:239));
    writeAdcVoltage(fileId, 'Bytes 240-241 EPS_APRS5V', b(240:241));
    writeAdcCurrent(fileId, 'Bytes 242-243 EPS_APRS5V_I', b(242:243));

    fprintf(fileId, '\n[CIP]\n');
    fprintf(fileId, 'Byte 244    CIP_Mode = %u (0x%02X, %s)\n', ...
        b(244), b(244), cipModeName(b(244)));
    writeCipTemperature(fileId, 'Bytes 245-246 CIP_Temp0', b(245:246));
    writeCipTemperature(fileId, 'Bytes 247-248 CIP_Temp1', b(247:248));
    writeCipTemperature(fileId, 'Bytes 249-250 CIP_Temp2', b(249:250));

    fprintf(fileId, '\n[HyperSCAN raw values - conversion not supplied]\n');
    writeU16(fileId, 'Bytes 251-252 HyperSCAN_SF2_I', b(251:252), 'raw uint16');
    writeU16(fileId, 'Bytes 253-254 HyperSCAN_SF2_T', b(253:254), 'raw uint16');
    writeU16(fileId, 'Bytes 255-256 HyperSCAN_CMOSV_1', b(255:256), 'raw uint16');
    writeU16(fileId, 'Bytes 257-258 HyperSCAN_CMOSV_2', b(257:258), 'raw uint16');
    writeU16(fileId, 'Bytes 259-260 HyperSCAN_CMOSV_3', b(259:260), 'raw uint16');
    writeU16(fileId, 'Bytes 261-262 HyperSCAN_CMOSV_4', b(261:262), 'raw uint16');
    writeU16(fileId, 'Bytes 263-264 HyperSCAN_CMOS_I', b(263:264), 'raw uint16');
    writeU16(fileId, 'Bytes 265-266 HyperSCAN_CMOS_T', b(265:266), 'raw uint16');

    fprintf(fileId, '\n[APRS]\n');
    fprintf(fileId, 'Byte 267    APRS_RX_V = %.1f V (raw=%u)\n', double(b(267)) / 10, b(267));
    fprintf(fileId, 'Byte 268    APRS_TX_V = %.1f V (raw=%u)\n', double(b(268)) / 10, b(268));
    aprsRadio = bitand(b(269), uint8(3));
    aprsDigipeater = bitand(bitshift(b(269), -2), uint8(3));
    fprintf(fileId, 'Byte 269    APRS_Configuration = 0x%02X, bits[7:0]=%s\n', ...
        b(269), bitsText(b(269)));
    fprintf(fileId, '  bits[1:0] Radio_on_off = %u (%s)\n', ...
        aprsRadio, aprsRadioName(aprsRadio));
    fprintf(fileId, '  bits[3:2] Digipeater_on_off = %u (%s)\n', ...
        aprsDigipeater, aprsDigipeaterName(aprsDigipeater));
    fprintf(fileId, '  bits[7:4] Reserved = 0x%X\n', bitshift(b(269), -4));
    writeScaledU16(fileId, 'Bytes 270-271 APRS_T', b(270:271), 0.1, 'scaled value');

    if numel(b) >= 272
        fprintf(fileId, '\n[Unmapped]\n');
        fprintf(fileId, 'Byte 272    Unmapped = 0x%02X\n', b(272));
    end

    if numel(b) >= 2
        fprintf(fileId, '\n[Frame FCS]\n');
        fprintf(fileId, 'Final two bytes = %02X %02X (little-endian FCS 0x%04X)\n', ...
            b(end - 1), b(end), receivedFcs);
        fprintf(fileId, 'Calculated CRC-16/X-25 = 0x%04X (%s)\n', ...
            calculatedFcs, fcsStatus);
        fprintf(fileId, 'FCS bytes are part of Raw_bit.txt but not part of the CRC input.\n');
    end
end


function mode = readBoundaryMode(frameDir)
    mode = 'STRICT';
    statusFile = fullfile(frameDir, 'Matlab_HDLC_Status.txt');
    if ~isfile(statusFile)
        % Backward-compatible fallback for outputs created before MATLAB
        % became the sole owner of strict/diagnostic HDLC deframing.
        statusFile = fullfile(frameDir, 'Boundary_Recovery_Status.txt');
    end
    if ~isfile(statusFile)
        return;
    end
    match = regexp(fileread(statusFile), ...
        '(?mi)^Boundary Mode:\s*([A-Z_]+)\s*$', 'tokens', 'once');
    if ~isempty(match)
        mode = match{1};
    end
end


function [withSsid, call, ssid] = decodeAx25Address(addressBytes)
    call = strtrim(char(bitshift(addressBytes(1:6), -1)));
    ssid = double(bitand(bitshift(addressBytes(7), -1), uint8(15)));
    withSsid = sprintf('%s-%d', call, ssid);
end


function writeBitGroup(fileId, byteNumber, title, value, names, zeroMeanings, oneMeanings)
    fprintf(fileId, 'Byte %-3d    %s = 0x%02X, bits[7:0]=%s\n', ...
        byteNumber, title, value, bitsText(value));
    for bitIndex = 0:7
        state = bitget(value, bitIndex + 1);
        if state == 0
            meaning = zeroMeanings{bitIndex + 1};
        else
            meaning = oneMeanings{bitIndex + 1};
        end
        fprintf(fileId, '  bit[%d] %-39s = %d (%s)\n', ...
            bitIndex, names{bitIndex + 1}, state, meaning);
    end
end


function text = bitsText(value)
    text = sprintf('%d%d%d%d%d%d%d%d', bitget(value, 8), bitget(value, 7), ...
        bitget(value, 6), bitget(value, 5), bitget(value, 4), bitget(value, 3), ...
        bitget(value, 2), bitget(value, 1));
end


function value = leU16(bytes)
    value = bitor(uint16(bytes(1)), bitshift(uint16(bytes(2)), 8));
end


function value = leI16(bytes)
    unsigned = double(leU16(bytes));
    if unsigned >= 32768
        value = int16(unsigned - 65536);
    else
        value = int16(unsigned);
    end
end


function value = leU32(bytes)
    value = bitor(uint32(bytes(1)), bitshift(uint32(bytes(2)), 8));
    value = bitor(value, bitshift(uint32(bytes(3)), 16));
    value = bitor(value, bitshift(uint32(bytes(4)), 24));
end


function value = leFloat32(bytes)
    value = double(typecast(leU32(bytes), 'single'));
end


function text = hexBytes(bytes)
    text = strtrim(sprintf('%02X ', bytes));
end


function text = chooseText(state, zeroText, oneText)
    if state == 0
        text = zeroText;
    else
        text = oneText;
    end
end


function writeScaledU16(fileId, label, bytes, scale, unitText)
    raw = leU16(bytes);
    value = double(raw) * scale;
    fprintf(fileId, '%s = %.9g %s (raw uint16=%u, bytes=%s)\n', ...
        label, value, unitText, raw, hexBytes(bytes));
end


function writeScaledI16(fileId, label, bytes, scale, unitText)
    raw = leI16(bytes);
    value = double(raw) * scale;
    fprintf(fileId, '%s = %.9g %s (raw int16=%d, bytes=%s)\n', ...
        label, value, unitText, raw, hexBytes(bytes));
end


function writeI8(fileId, label, byte, unitText)
    unsigned = double(byte);
    if unsigned >= 128
        value = unsigned - 256;
    else
        value = unsigned;
    end
    fprintf(fileId, '%s = %d %s (raw=0x%02X)\n', label, value, unitText, byte);
end


function writeAdcVoltage(fileId, label, bytes)
    rawWord = leU16(bytes);
    shifted = bitshift(rawWord, -3);
    volts = double(shifted) * 0.008;
    fprintf(fileId, '%s = %.6f V (raw=0x%04X, raw>>3=%u)\n', ...
        label, volts, rawWord, shifted);
end


function writeAdcCurrent(fileId, label, bytes)
    rawWord = leU16(bytes);
    shifted = bitshift(rawWord, -3);
    signedCode = double(shifted);
    if bitget(shifted, 13) ~= 0
        signedCode = signedCode - 8192;
    end
    amperes = signedCode * 40e-6 / 0.01;
    fprintf(fileId, '%s = %.6f A (raw=0x%04X, raw>>3=%u, signed13=%d)\n', ...
        label, amperes, rawWord, shifted, signedCode);
end


function writeCipTemperature(fileId, label, bytes)
    value = leU16(bytes);
    fprintf(fileId, ...
        '%s = %u (little-endian uint16, bytes=%s; source map says Uint8 for a 2-byte field)\n', ...
        label, value, hexBytes(bytes));
end


function crc = crc16X25(bytes)
    crc = uint16(hex2dec('FFFF'));
    polynomial = uint16(hex2dec('8408'));
    for value = uint8(bytes(:).')
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


function writeU16(fileId, label, bytes, unitText)
    value = leU16(bytes);
    fprintf(fileId, '%s = %u (0x%04X, %s)\n', label, value, value, unitText);
end


function writeI16(fileId, label, bytes, unitText)
    value = leI16(bytes);
    fprintf(fileId, '%s = %d (raw bytes %02X %02X, %s)\n', ...
        label, value, bytes(1), bytes(2), unitText);
end


function writeU32(fileId, label, bytes, unitText)
    value = leU32(bytes);
    fprintf(fileId, '%s = %u (0x%08X, %s)\n', label, value, value, unitText);
end


function text = unixTimeText(value)
    try
        time = datetime(double(value), 'ConvertFrom', 'posixtime', ...
            'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss ''UTC''');
        text = char(time);
    catch
        text = 'UTC conversion unavailable';
    end
end


function text = controlName(value)
    if value == hex2dec('03')
        text = 'UI frame';
    else
        text = 'unmapped control';
    end
end


function text = pidName(value)
    if value == hex2dec('F0')
        text = 'no layer-3 protocol';
    else
        text = 'unmapped PID';
    end
end


function text = packetApidName(value)
    switch double(value)
        case 0, text = 'Beacon';
        case 1, text = 'ADCS HK';
        otherwise, text = 'Unknown';
    end
end


function text = fswModeName(value)
    switch double(value)
        case hex2dec('00'), text = 'Phoenix';
        case hex2dec('01'), text = 'Safe';
        case hex2dec('02'), text = 'Charging';
        case hex2dec('03'), text = 'COMM_XBand';
        case hex2dec('0A'), text = 'SCI_HyperSCAN';
        case hex2dec('0B'), text = 'SCI_CIP';
        case hex2dec('0C'), text = 'Reentry';
        case hex2dec('FF'), text = 'Undefined';
        otherwise, text = 'Unknown';
    end
end


function text = attitudeEstimateName(value)
    names = {'No Attitude Estimation', 'MEMS rate sensing', 'Mag rate filter', ...
        'Mag rate filter with pitch estimation', 'Mag and Fine-sun TRIAD algorithm', ...
        'Full-state EKF', 'MEMS gyro EKF', 'User Coded Estimation Mode'};
    index = double(value) + 1;
    if index <= numel(names), text = names{index}; else, text = 'Reserved'; end
end


function text = adcsControlName(value)
    names = {'No control', 'Detumbling control', 'Y-Thomson spin', ...
        'Y-Wheel momentum stabilized - Initial Pitch Acquisition', ...
        'Y-Wheel momentum stabilized - Steady State', 'XYZ-Wheel control', ...
        'Rwheel sun tracking control', 'Rwheel target tracking control', ...
        'Very Fast-spin Detumbling control', 'Fast-spin Detumbling control', ...
        'User Specific Control Mode 1', 'User Specific Control Mode 2', ...
        'Stop R-wheels', 'User Coded Control Mode', ...
        'Yaw-only or roll-only sun tracking mode', ...
        'Target-tracking yaw-only wheel control Mode'};
    text = names{double(value) + 1};
end


function text = adcsRunModeName(value)
    names = {'off', 'Enable', 'Triggered', 'Simulation'};
    text = names{double(value) + 1};
end


function text = commissioningStepName(value)
    names = {'Initial Angular Rate Estimation', 'Detumbling', ...
        'Magnetometer Deployment', 'Magnetometer Calibration (update TLE)', ...
        'XYZ reaction wheel control', 'Sun Tracking 3-Axis Control', ...
        'None, all finished'};
    index = double(value) + 1;
    if index <= numel(names), text = names{index}; else, text = 'Reserved'; end
end


function text = commissioningErrorName(value)
    names = {'no error', 'IARE_faulty_main_magnetometer_error_flag', ...
        'IARE_unexpected_error_flag', 'IARE_angular_rate_too_fast_error_flag', ...
        'detumbling_error_flag', 'mag_deploy_error_flag', ...
        'xyzWheel_polarity_error_flag', 'cubeWheel_performance_error_flag', ...
        'fine_tuning_yWheel_bias_flag', ...
        'sunTracking_pointing_performance_error_flag'};
    index = double(value) + 1;
    if index <= numel(names), text = names{index}; else, text = 'Reserved'; end
end


function text = tcErrorName(value)
    names = {'No Error', 'Invalid telecommand ID', ...
        'Incorrect TC parameter length', 'Incorrect TC parameter value', ...
        'CRC check failed'};
    index = double(value) + 1;
    if index <= numel(names), text = names{index}; else, text = 'Reserved'; end
end


function text = readbackModeName(value)
    names = {'Decimal', 'Binary', 'Hex', 'Reserved'};
    text = names{double(value) + 1};
end


function text = uhfChannelName(value)
    names = {'437.5 MHz', '437.345 MHz', '436.5 MHz', '435.5 MHz'};
    text = names{double(value) + 1};
end


function text = cipModeName(value)
    switch double(value)
        case 7, text = 'CIP Idle Mode';
        case 8, text = 'CIP Normal Mode';
        case 9, text = 'CIP Fast Mode';
        otherwise, text = 'Unknown';
    end
end


function text = aprsRadioName(value)
    switch double(value)
        case 0, text = 'APRS Radio Module sleep';
        case 1, text = 'APRS Radio Module active';
        otherwise, text = 'Reserved';
    end
end


function text = aprsDigipeaterName(value)
    switch double(value)
        case 2, text = 'Digipeater function Enable';
        case 3, text = 'Digipeater function Disable';
        otherwise, text = 'Reserved';
    end
end
