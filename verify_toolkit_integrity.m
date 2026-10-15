function report=verify_toolkit_integrity()
% Check every shipped file against the SHA-256 manifest (including raw IQ).
root=string(fileparts(mfilename('fullpath')));
assert(usejava('jvm'),'SCIONX:Portable:JVM','MATLAB must run with Java enabled.');
manifest=jsondecode(fileread(fullfile(root,'PACKAGE_MANIFEST.json')));
total=0;
for k=1:numel(manifest.files)
    item=manifest.files(k);
    file=fullfile(root,string(item.path));
    info=dir(file);
    assert(isfile(file) && isscalar(info) && info.bytes==item.bytes, ...
        'SCIONX:Portable:File','Missing or changed file: %s',file);
    engine=java.security.MessageDigest.getInstance('SHA-256');
    fid=fopen(file,'rb'); assert(fid>=0,'Cannot read %s',file);
    closer=onCleanup(@() fclose(fid));
    while true
        block=fread(fid,1048576,'*uint8');
        if isempty(block), break; end
        engine.update(typecast(block,'int8'));
    end
    clear closer
    bytes=typecast(engine.digest(),'uint8');
    actual=lower(reshape(dec2hex(bytes,2).',1,[]));
    assert(strcmp(actual,item.sha256),'SCIONX:Portable:Hash','SHA-256 mismatch: %s',file);
    total=total+item.bytes;
end
report=struct('Status',"PASS",'Files',numel(manifest.files),'Bytes',total);
fprintf('PASS: %d package files verified by SHA-256.\n',report.Files);
end
