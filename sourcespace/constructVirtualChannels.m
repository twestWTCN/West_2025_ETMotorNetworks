function VCData = constructVirtualChannels(beamformer,ftdata,chansel,labellist)
% if isempty(labellist)
%     labellist =
% end

VCData = [];
VCData.label = labellist;
VCData.time = ftdata.time;
VCData.fsample = ftdata.fsample;
 svdtmptmp = {};
parfor vc = 1:numel(beamformer)
    vctmp = [];
    vctmp.label = {'x', 'y', 'z'};
    vctmp.time = ftdata.time;
    for i=1:length(ftdata.trial)
        vctmp.trial{i} = beamformer{vc} * ftdata.trial{i}(chansel,:);
    end
    vctmp = cat(2, vctmp.trial{:});

    u1 = svd(vctmp, 'econ');
    svdtmp = {};
    for k = 1:length(ftdata.trial)
        svdtmp{k}  = u1(:,1)' * beamformer{vc} * ftdata.trial{k}(chansel,:);
    end

    svdtmptmp{vc} = svdtmp;
        disp(vc./numel(beamformer))

end

% Now assign in
% VCData.trial = cell(1,numel(beamformer));
for vc = 1:numel(beamformer)
    for k = 1:length(ftdata.trial)
        VCData.trial{k}(vc,:) = svdtmptmp{vc}{k}; 
    end
end

a = [];
