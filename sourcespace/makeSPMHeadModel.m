function makeSPMHeadModel(R,spmDataPath,MRIPath,Fid)

if ~isempty(MRIPath)
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {MRIPath};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {MRIPath};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'Z:\TimWest\GITHUB\spm12\tpm\TPM.nii'};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    
    
    matlabbatch{2}.spm.meeg.source.headmodel.D = {spmDataPath};
    matlabbatch{2}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{2}.spm.meeg.source.headmodel.comment = '';
    [pathstr, name, ext] = fileparts(MRIPath);
    matlabbatch{2}.spm.meeg.source.headmodel.meshing.meshes.mri = {[pathstr '\w' name ext]};
    matlabbatch{2}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = Fid{1};
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type =  Fid{2};
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = Fid{3};
    matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{2}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
% %     matlabbatch{2}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    
else
    % Use Template MRI
    matlabbatch{1}.spm.meeg.source.headmodel.D = {spmDataPath};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = Fid{1};
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type =  Fid{2};
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = Fid{3};
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
%     matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    
end

spm_jobman('run',matlabbatch);