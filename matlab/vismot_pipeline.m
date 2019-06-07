load list
subj = 1;

subjectname = list{subj};

%% MRI preprocessing

vismot_anatomy_dicom2mgz(subjectname)
vismot_anatomy_mgz2bti(subjectname) % this step is interactive. will not be reproduced, but is reproducible (tested)

vismot_anatomy_headmodel(subjectname)
vismot_anatomy_sourcemodel3d(subjectname, 4)

%% MEG preprocessing

% setting up subject info.
vismot_subjinfo(subjectname);

% artifact selection? 
% work in progress: vismot_preproc_artifact

% preprocess emptyroom data
vismot_execute_pipeline('vismot_emptyroom_script', subjectname);

% defining trials and rejecting artifacts
vismot_execute_pipeline('vismot_preproc_script', subjectname)

%% behavioral effects
vismot_rt

%% post cue power
frequencies = [10 22 38 42 58 62];
smoothing  = [2, 8, 8, 8, 8, 8];
for f=1:numel(frequencies)
vismot_execute_pipeline('vismot_bf_script', subjectname, {'latoi', 'post'}, ...
  {'frequency', frequencies(f)}, {'subjectname', subjectname}, {'smoothing', smoothing(f)},...
  {'lambda', []}, {'prewhiten', true})
end

vismot_execute_function('vismot_script_Tstat_post', [])



%% pre cue power
frequencies = [10 22 38 42 58 62];
smoothing  = [2, 8, 8, 8, 8, 8];
for f=1:numel(frequencies)
vismot_execute_pipeline('vismot_bf_script', subjectname, {'latoi', 'pre'}, ...
  {'frequency', frequencies(f)}, {'subjectname', subjectname}, {'smoothing', smoothing(f)},...
  {'lambda', []}, {'prewhiten', true})
end

vismot_execute_function('vismot_script_Tstat_pre', [])


%% pre cue coherence
vismot_execute_pipeline('vismot_bf_coh_roi_script', subjectname) % check this




