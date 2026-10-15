function result=scionx_reverse_hdlc(raw)
% raw is chronological NRZ bits ending immediately BEFORE the closing flag.
% Prefer a complete frame. Otherwise retain only a conservatively destuffed
% suffix; never reverse the transmission bit order or repair invalid bits.
raw=logical(raw(:)); flag=[0;1;1;1;1;1;1;0];
r=struct('Bytes',uint8([]),'Complete',false,'CRC',"UNVERIFIABLE_PARTIAL", ...
    'StopReason',"left_boundary_unknown",'StuffedZeros',0,'LeadingPartialBits',0, ...
    'RawStartIndex',NaN,'RawBitCount',numel(raw),'DestuffedBits',false(0,1));
hits=find(conv(double(raw),flipud(2*flag-1),'valid')==6);
if ~isempty(hits)
    opening=hits(end); data=raw(opening+8:end);
    [bits,valid,removed]=forward(data);
    if valid && mod(numel(bits),8)==0 && numel(bits)>=16
        r.Bytes=pack(bits); r.Complete=true; r.StopReason="opening_flag";
        r.StuffedZeros=removed; r.RawStartIndex=opening+8; r.DestuffedBits=bits;
        expected=bitor(uint16(r.Bytes(end-1)),bitshift(uint16(r.Bytes(end)),8));
        if scionx_crc16_x25(r.Bytes(1:end-2))==expected, r.CRC="PASS"; else, r.CRC="FAIL"; end
        result=r; return
    end
end
% Process zero and its left-hand run together. A zero preceded by exactly
% five known ones is stuffed. Six or more ones / unknown left context stop.
i=numel(raw); output=false(size(raw)); cursor=numel(raw)+1;
while i>=1
    if raw(i)
        j=i; while j>=1 && raw(j), j=j-1; end
        count=i-j;
        if j==0, r.StopReason="left_boundary_unknown"; break; end
        if count>=5, r.StopReason="missing_stuffed_zero_or_abort"; break; end
        output(cursor-count:cursor-1)=true; cursor=cursor-count; i=j;
    else
        j=i-1; while j>=1 && raw(j), j=j-1; end
        count=i-j-1;
        if j==0, r.StopReason="left_boundary_unknown"; break; end
        if count>=6, r.StopReason="stuffing_violation"; break; end
        if count==5
            output(cursor-5:cursor-1)=true; cursor=cursor-5; r.StuffedZeros=r.StuffedZeros+1;
        else
            group=[true(count,1);false];
            output(cursor-numel(group):cursor-1)=group; cursor=cursor-numel(group);
        end
        i=j;
    end
end
bits=output(cursor:end); r.RawStartIndex=i+1;
r.LeadingPartialBits=mod(numel(bits),8); r.DestuffedBits=bits;
r.Bytes=pack(bits(r.LeadingPartialBits+1:end)); result=r;
end

function [output,valid,removed]=forward(input)
output=false(size(input)); n=0; ones=0; removed=0; valid=true;
for j=1:numel(input)
    if ones==5
        if input(j), valid=false; return; end
        removed=removed+1; ones=0; continue
    end
    n=n+1; output(n)=input(j);
    if input(j), ones=ones+1; else, ones=0; end
end
valid=ones~=5; output=output(1:n);
end

function bytes=pack(bits)
bytes=uint8((2.^(0:7))*reshape(double(bits),8,[]));
end
