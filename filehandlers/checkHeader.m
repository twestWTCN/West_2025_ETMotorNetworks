function flag = checkHeader(ftdata)
flag = 0;
A = numel(ftdata.label);
if numel(ftdata.hdr.label) ~= A; flag = 1; end
if ftdata.hdr.nChans ~= A; flag = 1; end
if numel(ftdata.hdr.chantype) ~= A; flag = 1; end
if numel(ftdata.hdr.chanunit) ~= A; flag = 1; end

if flag
    warning('Your data header is corrupt!')
end