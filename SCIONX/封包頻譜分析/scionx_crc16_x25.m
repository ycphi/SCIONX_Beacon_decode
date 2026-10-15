function crc=scionx_crc16_x25(bytes)
% Reflected CRC-16/X-25, init FFFF, poly 8408, xorout FFFF.
crc=uint16(65535);
for value=uint8(bytes(:))'
    crc=bitxor(crc,uint16(value));
    for k=1:8
        if bitand(crc,1), crc=bitxor(bitshift(crc,-1),uint16(hex2dec('8408')));
        else, crc=bitshift(crc,-1); end
    end
end
crc=bitxor(crc,uint16(65535));
end
