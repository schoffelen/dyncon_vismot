load list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE SUBJECT ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subj = input('please input subject number (1-19) and press enter')

subjectname = list{subj};

%% MRI preprocessing
vismot_anatomy_dicom2mgz(subjectname)
vismot_anatomy_mgz2bti(subjectname) % this step is interactive. will not be reproduced, but is reproducible (tested)

vismot_anatomy_headmodel(subjectname)
vismot_anatomy_sourcemodel3d(subjectname, 4)

%% MEG preprocessing
% setting up subject info.
vismot_subjinfo(subjectname);

% artifact selection
vismot_execute_pipeline('vismot_preproc_artifact', subjectname)

% preprocess emptyroom data
vismot_execute_pipeline('vismot_emptyroom_script', subjectname);

% defining trials and rejecting artifacts
vismot_execute_pipeline('vismot_preproc_script', subjectname)


%% post cue power
frequencies = [6 10 22 38 42 58 62, 78, 82];
smoothing  = [4, 4, 8, 8, 8, 8, 8, 8, 8];
% stratified for RT
for f=1:numel(frequencies)
  vismot_execute_pipeline('vismot_bf_script', subjectname, {'latoi', 'post'}, ...
    {'frequency', frequencies(f)}, {'subjectname', subjectname}, {'smoothing', smoothing(f)},...
    {'lambda', []}, {'prewhiten', true}, {'stratifyflag', true})
end

% not stratified
for f=1:numel(frequencies)
  vismot_execute_pipeline('vismot_bf_script', subjectname, {'latoi', 'post'}, ...
    {'frequency', frequencies(f)}, {'subjectname', subjectname}, {'smoothing', smoothing(f)},...
    {'lambda', []}, {'prewhiten', true})
end


%% pre cue power
frequencies = [6, 10 22 38 42 58 62, 78, 82];
smoothing  = [2, 2, 8, 8, 8, 8, 8, 8, 8];
for f=1:numel(frequencies)
  vismot_execute_pipeline('vismot_bf_script', subjectname, {'latoi', 'pre'}, ...
    {'frequency', frequencies(f)}, {'subjectname', subjectname}, {'smoothing', smoothing(f)},...
    {'lambda', []}, {'prewhiten', true})
end

%% pre cue coherence
frequencies = [6, 10 22 38 42 58 62, 78 82];
smoothing  = [2, 2, 8, 8, 8, 8, 8, 8, 8];
for f = 1:numel(frequencies)
  vismot_execute_pipeline('vismot_bf_coh_roi_script', subjectname, {'toi', 'pre'},...
    {'conditions', 'previous'}, {'frequency', frequencies(f), {'smoothing', smoothing(f)}}, ...
    {'lambda', '100%'}, {'nrand', 0}, {'roi_to', 'roi'}) 
end

%%%%%%%%%%%%%%%%%%%
%% GROUP ANALYSIS %
%%%%%%%%%%%%%%%%%%%
% behavioral effects
vismot_rt

% post cue power
vismot_execute_pipeline('vismot_script_Tstat_wholebrain', [], {'stratifyflag', true}, {'toi', 'post'}) % stratified for RT
vismot_execute_pipeline('vismot_script_Tstat_wholebrain', [], {'stratifyflag', false}, {'toi', 'post'}) % not stratified
vismot_execute_pipeline('vismot_define_roi', []);

% pre cue power
vismot_execute_pipeline('vismot_script_Tstat_wholebrain', [],  {'stratifyflag', false}, {'toi', 'pre'})
vismot_execute_pipeline('vismot_script_Tstat_roi', [])

% correlation between sensorimotor power and behavior
vismot_singletrialcorr

% pre cue coherence
vismot_execute_pipeline('vismot_script_Tstat_coh', [])


%%%%%%%%%%%%%
%% PLOTTING %
%%%%%%%%%%%%%
vismot_figures
