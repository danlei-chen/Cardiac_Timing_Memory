close all;

%Gets subject/trial/eccentricity information from user
prompt = {'Subject ID:','Shimmer enabled (True/False)','Minimum R Prominence'}; %Sets up a prompt to collect information for the experiment
num_lines = 3;
dialogue_title = 'Experiment Start'; %Dialogue Box Title
def = {'000','true','0.7'}; %Default Responses
answer = inputdlg(prompt,dialogue_title,num_lines,def);
pa.subjid = char(answer(1,1));
pa.shimmer = char(answer(2,1));
pa.Minimum_R_Prominence_mV=str2double(char(answer(3,1)));
% pa.subjid='000';
% pa.shimmer='true';
% pa.Minimum_R_Prominence_mV=0.7;
  
Screen('Preference', 'SkipSyncTests', 1); 
% screenNum = max(Screen('Screens'));
screenNum = 1;
[wPtr, rect] = Screen('OpenWindow', screenNum);
% [wPtr,rect] = Screen('OpenWindow',screenNum,0,[0 0 800 600]);

% Set up dir
if (ismac)
    base_dir = '/Users/chendanlei/Google Drive/Interoceptive Rinstatement/InteroReins_Matlab/';
    output_dir = sprintf('%ssubj-%s_InteroRein_encoding_output/',base_dir,pa.subjid);
elseif (IsWin)
    base_dir = 'C:\Users\LFBadmin\Desktop\InteroReins_Matlab\';
    output_dir = sprintf('%ssubj-%s_InteroRein_encoding_output\\',base_dir,pa.subjid);
end
cd(base_dir)
mkdir(output_dir);

if strcmpi(pa.shimmer, 'true')
    pa.comPort_ECG = sprintf('%d', 8);
    pa.comPort_EDA = sprintf('%d', 9);
    pa.sampling_rate_ECG_Hz = 1024;
    % This is how much time between checking the Shimmer for data (smaller numbers increase temporal precision but risk an error if there are no data available)
    pa.sampling_period_check_Shimmer_sec = 0.050;
    % If time to beep is THIS long then there is clearly a problem
    MAX_TIME_TO_BEEP_SEC = 0.8;
    % For detecting R spikes
    Minimum_RR_Interval_sec = 0.5;
    % Don't wait longer than this time between R spikes. If it's been this long then there is something wrong and the trial will be repeated. DEFAULT = 3 just so there has to be a major error before aborting the trial
    Maximum_RR_Interval_sec = 3;
    Minimum_R_Prominence_mV = pa.Minimum_R_Prominence_mV; % run resting state and see 
    Npeaks_found_that_do_not_count = 0;
end
 
%%
rng('shuffle')
 
%%% Read in pics
%face images
face_stim={};face_label = {};
load([base_dir, '/stim_mat/CFD_faceStim/Face_AF_image.mat'])
face_stim{end+1}=images;
face_label{end+1}=repmat("female",1,size(images,3));
load([base_dir, '/stim_mat/CFD_faceStim/Face_AM_image.mat'])
face_stim{end+1}=images;
face_label{end+1}=repmat("male",1,size(images,3));
load([base_dir, '/stim_mat/CFD_faceStim/Face_BF_image.mat'])
face_stim{end+1}=images;
face_label{end+1}=repmat("female",1,size(images,3));
load([base_dir, '/stim_mat/CFD_faceStim/Face_BM_image.mat'])
face_stim{end+1}=images;
face_label{end+1}=repmat("male",1,size(images,3));
load([base_dir, '/stim_mat/CFD_faceStim/Face_WF_image.mat'])
face_stim{end+1}=images;
face_label{end+1}=repmat("female",1,size(images,3));
load([base_dir, '/stim_mat/CFD_faceStim/Face_WM_image.mat'])
face_stim{end+1}=images;
face_label{end+1}=repmat("male",1,size(images,3));
% load([base_dir, '/stim_mat/CFD_faceStim/Face_LF_image.mat'])
% face_stim{end+1}=images;
% load([base_dir, '/stim_mat/CFD_faceStim/Face_LM_image.mat'])
% face_stim{end+1}=images;
%scene images
scene_stim={};scene_label = {};
load([base_dir, '/stim_mat/LMSun_sceneStim/Scene_indoor_image.mat']);
scene_stim{end+1}=images;
scene_label{end+1}=repmat("indoor",1,size(images,3));
load([base_dir, '/stim_mat/LMSun_sceneStim/Scene_outdoor_image.mat']);
scene_stim{end+1}=images;
scene_label{end+1}=repmat("outdoor",1,size(images,3));

%block info and other parameters
pa.enc.mini_blocks = true;
pa.enc.num_blocks = 24; %needs to be a common multiple of length(pa.enc.num_block_images), length(pa.enc.cross_block_HB_ITI), and pa.enc.num_face_subcategory
pa.enc.num_block_images = [4,5,6];
pa.enc.in_block_HB_ITI = [1,2];
pa.enc.in_block_HB_rep_ITI = 1;
pa.enc.cross_block_HB_ITI = [3,4,5];
pa.enc.image_duration = 0.12;
pa.enc.systole_timing = 0.18; %timing from Yang et al., 2017
pa.enc.diastole_timing = 0.48;
pa.enc.image_repetition = 4;
pa.enc.num_image_category = 2; %faces and scenes
pa.enc.num_face_subcategory = length(face_stim); 
pa.enc.num_scene_subcategory = length(scene_stim); 
pa.enc.num_blocks_each_category = pa.enc.num_blocks/pa.enc.num_image_category;
active_keys = [KbName('k'),KbName('l'),KbName('n'),KbName('m')];
 
%encoding stim randomization
pa.enc.num_total_images = pa.enc.num_blocks * mean(pa.enc.num_block_images);
pa.enc.num_total_trials = pa.enc.num_total_images;
pa.enc.num_image_in_each_category = pa.enc.num_total_images/pa.enc.num_image_category;
pa.enc.num_face_in_each_subcategory = pa.enc.num_image_in_each_category/pa.enc.num_face_subcategory;
pa.enc.num_scene_in_each_subcategory = pa.enc.num_image_in_each_category/pa.enc.num_scene_subcategory;
 
pa.rec.num_new_images = 24;
pa.rec.num_total_images = pa.rec.num_new_images+pa.enc.num_total_images;% show both old and new images in recognition
pa.rec.num_total_trials = pa.rec.num_total_images;
pa.rec.num_new_image_in_each_category = pa.rec.num_new_images/pa.enc.num_image_category;
pa.rec.num_new_face_in_each_subcategory = pa.rec.num_new_image_in_each_category/pa.enc.num_face_subcategory;
pa.rec.num_new_scene_in_each_subcategory = pa.rec.num_new_image_in_each_category/pa.enc.num_scene_subcategory;
pa.rec.in_block_ITI = [3,4,5];
pa.rec.image_repetition = pa.enc.image_repetition;
pa.rec.image_duration = pa.enc.image_duration;

pa.total_image_used = pa.enc.num_total_images+pa.rec.num_new_images;
pa.total_image_used_in_each_category = pa.enc.num_image_in_each_category + pa.rec.num_new_image_in_each_category;
pa.total_face_used_in_each_subcategory = pa.enc.num_face_in_each_subcategory + pa.rec.num_new_face_in_each_subcategory;
pa.total_scene_used_in_each_subcategory = pa.enc.num_scene_in_each_subcategory + pa.rec.num_new_scene_in_each_subcategory;

%randomly selecting from subcategory
stim.enc.face=[];stim.rec.face=[];
stim.enc.face_label=[];stim.rec.face_label=[];
stim.rec.face_novelty_marker=[];stim.rec.face_congruency_marker=[];
for n=1:pa.enc.num_face_subcategory
    num_image_this_subcategory = size(face_stim{n},3);
    random_selection = randperm(num_image_this_subcategory, pa.total_face_used_in_each_subcategory);
    random_img = face_stim{n}(:,:,random_selection);
    random_label = face_label{n}(random_selection);
    %part of the randomized images are assined as old, and the rest and marked as novel
    stim.enc.face = cat(3,stim.enc.face,random_img(:,:,1:pa.enc.num_face_in_each_subcategory));
    stim.rec.face = cat(3,stim.rec.face,random_img);
    stim.enc.face_label = [stim.enc.face_label,random_label(1:pa.enc.num_face_in_each_subcategory)];
    stim.rec.face_label = [stim.rec.face_label,random_label];
    stim.rec.face_novelty_marker = [stim.rec.face_novelty_marker,[repmat("old",1,pa.enc.num_face_in_each_subcategory),repmat("new",1,pa.rec.num_new_face_in_each_subcategory)]];
    stim.rec.face_congruency_marker = [stim.rec.face_congruency_marker,[repmat("congruent",1,pa.enc.num_face_in_each_subcategory/2),repmat("incongruent",1,pa.enc.num_face_in_each_subcategory/2),repmat("congruent",1,pa.rec.num_new_face_in_each_subcategory/2),repmat("incongruent",1,pa.rec.num_new_face_in_each_subcategory/2)]];
end
stim.enc.scene=[];stim.rec.scene=[];
stim.enc.scene_label=[];stim.rec.scene_label=[];
stim.rec.scene_novelty_marker=[];stim.rec.scene_congruency_marker=[];
for n=1:pa.enc.num_scene_subcategory
    num_image_this_subcategory = size(scene_stim{n},3);
    random_selection = randperm(num_image_this_subcategory, pa.total_scene_used_in_each_subcategory);
    random_img = scene_stim{n}(:,:,random_selection);
    random_label = scene_label{n}(random_selection);
    %part of the randomized images are assined as old, and the rest and marked as novel
    stim.enc.scene = cat(3,stim.enc.scene,random_img(:,:,1:pa.enc.num_scene_in_each_subcategory));
    stim.rec.scene = cat(3,stim.rec.scene,random_img);
    stim.enc.scene_label = [stim.enc.scene_label,random_label(1:pa.enc.num_scene_in_each_subcategory)];
    stim.rec.scene_label = [stim.rec.scene_label,random_label];
    stim.rec.scene_novelty_marker = [stim.rec.scene_novelty_marker,[repmat("old",1,pa.enc.num_scene_in_each_subcategory),repmat("new",1,pa.rec.num_new_scene_in_each_subcategory)]];
    stim.rec.scene_congruency_marker = [stim.rec.scene_congruency_marker,[repmat("congruent",1,pa.enc.num_scene_in_each_subcategory/2),repmat("incongruent",1,pa.enc.num_scene_in_each_subcategory/2),repmat("congruent",1,pa.rec.num_new_scene_in_each_subcategory/2),repmat("incongruent",1,pa.rec.num_new_scene_in_each_subcategory/2)]];
end  
%randomize encoding order
random_num = randperm(size(stim.enc.face,3));
stim.enc.face = stim.enc.face(:,:,random_num);
stim.enc.face_label = stim.enc.face_label(random_num);
random_num = randperm(size(stim.enc.scene,3));
stim.enc.scene = stim.enc.scene(:,:,random_num);
stim.enc.scene_label = stim.enc.scene_label(random_num);

%within block ITI sequence (not counterbalanced within block, but within each category blocks)
category_block_ITI_sequence = repmat(pa.enc.in_block_HB_ITI, 1, (pa.enc.num_total_trials/length(pa.enc.in_block_HB_ITI)/pa.enc.num_image_category));
face_block_ITI_sequence = category_block_ITI_sequence(randperm(length(category_block_ITI_sequence)));
scene_block_ITI_sequence = category_block_ITI_sequence(randperm(length(category_block_ITI_sequence)));

%generalize presentation sequence
num_image_per_face_block = repmat(pa.enc.num_block_images,1,(pa.enc.num_blocks/pa.enc.num_image_category/length(pa.enc.num_block_images)));
num_image_per_face_block = num_image_per_face_block(randperm(length(num_image_per_face_block)));
stim.enc.face_block={};
stim.enc.face_label_block={};
seq.enc.face_block_ITI={};
seq.enc.face_block_cumulative_ITI={};
 
num_image_per_scene_block = repmat(pa.enc.num_block_images,1,(pa.enc.num_blocks/pa.enc.num_image_category/length(pa.enc.num_block_images)));
num_image_per_scene_block = num_image_per_scene_block(randperm(length(num_image_per_scene_block)));
stim.enc.scene_block={};
stim.enc.scene_label_block={};
seq.enc.scene_block_ITI={};
seq.enc.scene_block_cumulative_ITI={};
 
for n = 1:pa.enc.num_blocks_each_category
    for m = 1:num_image_per_face_block(n)
        stim.enc.face_block{n} = stim.enc.face(:,:,(sum(num_image_per_face_block(1:n-1))+1:sum(num_image_per_face_block(1:n))));
        stim.enc.face_label_block{n} = stim.enc.face_label(sum(num_image_per_face_block(1:n-1))+1:sum(num_image_per_face_block(1:n)));
        seq.enc.face_block_ITI{n} = face_block_ITI_sequence((sum(num_image_per_face_block(1:n-1))+1:sum(num_image_per_face_block(1:n))));
        seq.enc.face_block_cumulative_ITI{n} = cumsum(face_block_ITI_sequence((sum(num_image_per_face_block(1:n-1))+1:sum(num_image_per_face_block(1:n)))));
    end
    for m = 1:num_image_per_scene_block(n)
        stim.enc.scene_block{n} = stim.enc.scene(:,:,(sum(num_image_per_scene_block(1:n-1))+1:sum(num_image_per_scene_block(1:n))));
        stim.enc.scene_label_block{n} = stim.enc.scene_label(sum(num_image_per_scene_block(1:n-1))+1:sum(num_image_per_scene_block(1:n)));
        seq.enc.scene_block_ITI{n} = scene_block_ITI_sequence((sum(num_image_per_scene_block(1:n-1))+1:sum(num_image_per_scene_block(1:n))));
        seq.enc.scene_block_cumulative_ITI{n} = cumsum(scene_block_ITI_sequence((sum(num_image_per_scene_block(1:n-1))+1:sum(num_image_per_scene_block(1:n)))));
    end
end
 
%cross block ITI sequence
seq.enc.cross_block_ITI_sequence = repmat(pa.enc.cross_block_HB_ITI,1,(pa.enc.num_blocks/length(pa.enc.cross_block_HB_ITI)));
seq.enc.cross_block_ITI_sequence = seq.enc.cross_block_ITI_sequence(randperm(length(seq.enc.cross_block_ITI_sequence)));
 
%block category sequence
seq.enc.block_category_sequence = repmat(["face","scene"], 1, (pa.enc.num_blocks/pa.enc.num_image_category));
%pseudo-randomize such that a category block doesn't repeat more than twice
redo_randomization=true;
while redo_randomization == true
    redo_randomization=false;
    seq.enc.block_category_sequence = seq.enc.block_category_sequence(randperm(length(seq.enc.block_category_sequence)));
    for n = 1:pa.enc.num_blocks
        if n < pa.enc.num_blocks-1 && seq.enc.block_category_sequence(n) == seq.enc.block_category_sequence(n+1) 
            if seq.enc.block_category_sequence(n) == seq.enc.block_category_sequence(n+2)
                redo_randomization = true;
                break;
            end
        end
    end
end
 
%even number participants associate face with systole, scene with diastole
seq.enc.cardiac_cycle_block_sequence = seq.enc.block_category_sequence;
if mod(pa.subjid,2) == 0
    seq.enc.cardiac_cycle_block_sequence(strcmp(seq.enc.block_category_sequence,'face')) = 'systole';
    seq.enc.cardiac_cycle_block_sequence(strcmp(seq.enc.block_category_sequence,'scene')) = 'diastole';
else
    seq.enc.cardiac_cycle_block_sequence(strcmp(seq.enc.block_category_sequence,'face')) = 'diastole';
    seq.enc.cardiac_cycle_block_sequence(strcmp(seq.enc.block_category_sequence,'scene')) = 'systole';
end 
seq.enc.cardiac_timing_block_sequence(strcmp(seq.enc.cardiac_cycle_block_sequence,'systole')) = pa.enc.systole_timing;
seq.enc.cardiac_timing_block_sequence(strcmp(seq.enc.cardiac_cycle_block_sequence,'diastole')) = pa.enc.diastole_timing;
 
%integrate different category blocks into all continuous blocks
stim.enc.block_stim = cell(1,pa.enc.num_blocks);
stim.enc.block_stim(strcmp(seq.enc.block_category_sequence,'face')) = stim.enc.face_block;
stim.enc.block_stim(strcmp(seq.enc.block_category_sequence,'scene')) = stim.enc.scene_block;
seq.enc.block_stim = stim.enc.block_stim;

stim.enc.block_stim_label = cell(1,pa.enc.num_blocks);
stim.enc.block_stim_label(strcmp(seq.enc.block_category_sequence,'face')) = stim.enc.face_label_block;
stim.enc.block_stim_label(strcmp(seq.enc.block_category_sequence,'scene')) = stim.enc.scene_label_block;
seq.enc.block_stim_label = stim.enc.block_stim_label;

%integrate all block ITI
seq.enc.in_block_ITI = cell(1,pa.enc.num_blocks);
seq.enc.in_block_ITI(strcmp(seq.enc.block_category_sequence,'face')) = seq.enc.face_block_ITI;
seq.enc.in_block_ITI(strcmp(seq.enc.block_category_sequence,'scene')) = seq.enc.scene_block_ITI;
seq.enc.in_block_cumulative_ITI = cell(1,pa.enc.num_blocks);
seq.enc.in_block_cumulative_ITI(strcmp(seq.enc.block_category_sequence,'face')) = seq.enc.face_block_cumulative_ITI;
seq.enc.in_block_cumulative_ITI(strcmp(seq.enc.block_category_sequence,'scene')) = seq.enc.scene_block_cumulative_ITI;
 
%integrate num images per block
seq.enc.num_images_per_block = zeros(1,pa.enc.num_blocks);
seq.enc.num_images_per_block(strcmp(seq.enc.block_category_sequence,'face')) = num_image_per_face_block;
seq.enc.num_images_per_block(strcmp(seq.enc.block_category_sequence,'scene')) = num_image_per_scene_block;

%randomize recognition trial order
perm_rec = randperm(size(stim.rec.face,3));
stim.rec.face = stim.rec.face(:,:,perm_rec);
stim.rec.face_novelty_marker = stim.rec.face_novelty_marker(perm_rec); 
stim.rec.face_congruency_marker = stim.rec.face_congruency_marker(perm_rec); 
perm_rec = randperm(size(stim.rec.scene,3));
stim.rec.scene = stim.rec.scene(:,:,perm_rec);
stim.rec.scene_novelty_marker = stim.rec.scene_novelty_marker(perm_rec); 
stim.rec.scene_congruency_marker = stim.rec.scene_congruency_marker(perm_rec); 

seq.rec.trial = cat(3,stim.rec.face,stim.rec.scene);
seq.rec.trial_category = [repmat("face",1,size(stim.rec.face,3)),repmat("scene",1,size(stim.rec.scene,3))];
seq.rec.trial_label = [stim.rec.face_label,stim.rec.scene_label];
seq.rec.trial_novelty_marker = [stim.rec.face_novelty_marker,stim.rec.scene_novelty_marker];
seq.rec.trial_congruency_marker = [stim.rec.face_congruency_marker,stim.rec.scene_congruency_marker];

%save stim.mat
save(sprintf('%sstimuli.mat',output_dir), 'stim')

perm_rec = randperm(size(seq.rec.trial,3));
seq.rec.trial = seq.rec.trial(:,:,perm_rec);
seq.rec.trial_novelty_marker = seq.rec.trial_novelty_marker(perm_rec);
seq.rec.trial_congruency_marker = seq.rec.trial_congruency_marker(perm_rec);
seq.rec.trial_category = seq.rec.trial_category(perm_rec);
seq.rec.trial_label = seq.rec.trial_label(perm_rec);

seq.rec.cardiac_cycle_sequence = seq.rec.trial_category;
if mod(pa.subjid,2) == 0
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'face') & strcmp(seq.rec.trial_congruency_marker,'congruent')) = 'systole';
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'face') & strcmp(seq.rec.trial_congruency_marker,'incongruent')) = 'diastole';
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'scene') & strcmp(seq.rec.trial_congruency_marker,'congruent')) = 'diastole';
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'scene') & strcmp(seq.rec.trial_congruency_marker,'incongruent')) = 'systole';
else
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'face') & strcmp(seq.rec.trial_congruency_marker,'congruent')) = 'diastole';
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'face') & strcmp(seq.rec.trial_congruency_marker,'incongruent')) = 'systole';
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'scene') & strcmp(seq.rec.trial_congruency_marker,'congruent')) = 'systole';
    seq.rec.cardiac_cycle_sequence(strcmp(seq.rec.trial_category,'scene') & strcmp(seq.rec.trial_congruency_marker,'incongruent')) = 'diastole';

end 
seq.rec.cardiac_timing_sequence(strcmp(seq.rec.cardiac_cycle_sequence,'systole')) = pa.enc.systole_timing;
seq.rec.cardiac_timing_sequence(strcmp(seq.rec.cardiac_cycle_sequence,'diastole')) = pa.enc.diastole_timing;
 
seq.rec.in_block_ITI = repmat(pa.rec.in_block_ITI, 1, pa.total_image_used/length(pa.rec.in_block_ITI));
seq.rec.in_block_ITI = seq.rec.in_block_ITI(randperm(length(seq.rec.in_block_ITI)));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Presentation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%black background
Screen('FillRect', wPtr, [0 0 0], rect);
% HideCursor; 
Screen('TextFont', wPtr, 'Helvetica'); 
xCenter = rect(3)/2;
yCenter = rect(4)/2;

%draw fixation
Screen('Flip',wPtr);
Screen('TextSize', wPtr, 60);
waitMessage = 'Please wait... Connecting...';
DrawFormattedText(wPtr, waitMessage, 'center', 'center',[255,255,255]);
Screen('Flip',wPtr);
 
%%
%%%%%%%%%%%%%%%%%%%%% Connect to Shimmer device %%%%%%%%%%%%%%%%%%%%%
if strcmpi(pa.shimmer, 'true')
    pa.comPort_ECG = sprintf('%d', 8);
    pa.comPort_ECG = sprintf('%d', 8);
    
    pa.sampling_rate_ECG_Hz = 1024;
    % This is how much time between checking the Shimmer for data (smaller numbers increase temporal precision but risk an error if there are no data available)
    pa.sampling_period_check_Shimmer_sec = 0.050;
    % If time to beep is THIS long then there is clearly a problem
    MAX_TIME_TO_BEEP_SEC = 0.8;
    % For detecting R spikes
    Minimum_RR_Interval_sec = 0.5;
    % Don't wait longer than this time between R spikes. If it's been this long then there is something wrong and the trial will be repeated. DEFAULT = 3 just so there has to be a major error before aborting the trial
    Maximum_RR_Interval_sec = 3;
    Minimum_R_Prominence_mV = pa.Minimum_R_Prominence_mV; % run resting state and see 

    % Connect to Shimmer at specified COM port
    fprintf('\n\nAttempting to connect to Shimmer, this may take 10 sec...\n');
    shimmer_ECG = ShimmerHandleClass(pa.comPort_ECG);
    % assign user friendly macros for setenabledsensors
    SensorMacros = SetEnabledSensorsMacrosClass;  
    
    % Ensure connection
    if (shimmer_ECG.connect)
        fprintf('\n\nSuccessfully connected to Shimmer device!\n');
        % Select sampling rate
        shimmer_ECG.setsamplingrate(pa.sampling_rate_ECG_Hz); 
        % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
        shimmer_ECG.setinternalboard('ECG');                                       
        % Disable other sensors
        shimmer_ECG.disableallsensors;                                             
        % Enable SENSOR_EXG1 and SENSOR_EXG2
        shimmer_ECG.setenabledsensors(SensorMacros.ECG,1);
    else
        error('Could not connect to Shimmer device');
    end
 
    % ECG Initialize
    % This is how much time between checking the Shimmer for data (smaller
    % numbers increase temporal precision but risk an error if there are no
    % data available)
    sampling_period_check_Shimmer_sec = 0.050;
 
    % ECG channel used to detect R spikes
    string_ECG_signal_HBD = 'ECG LA-RA';
    
    % ECG output files
    FILE_OUT_ECG = fopen(sprintf('%sencoding-ECG.csv',output_dir),'w');
    % make header
    cHeader = {'Time(Sec)' 'ECG_Raw(mV)'}; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    % write header to file
    fprintf(FILE_OUT_ECG,'%s\n',textHeader);
    
    % ECG output files
    FILE_OUT_ECG_Adj = fopen(sprintf('%sencoding-ECG-Baseline_Adjusted.csv',output_dir),'w');
    % make header
    cHeader = {'Time(Sec)' 'ECG_Raw(mV)','Peaks(mV)'}; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';   
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    % write header to file
    fprintf(FILE_OUT_ECG_Adj,'%s\n',textHeader);
    
    % ECG output files
    FILE_OUT_ECG_RSPIKES = fopen(sprintf('%sencoding-ECG-Rspikes.csv',output_dir),'w');
    % make header
    cHeader = {'Time(Sec)' 'ECG_Raw(mV)' 'Prominence(mV)' 'IBI(sec)'}; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    % write header to file
    fprintf(FILE_OUT_ECG_RSPIKES,'%s\n',textHeader);
end

%% encoding

%block level 
seq.enc.block_start_time = zeros(1,pa.enc.num_blocks); 
seq.enc.trial_start_time = cell(1,pa.enc.num_blocks);
seq.enc.stimOnset = cell(1,pa.enc.num_blocks);
seq.enc.real_trial_onset_relative_to_block_start_time = cell(1,pa.enc.num_blocks);
 
RestrictKeysForKbCheck(active_keys);
rsp.enc.response = cell(1,pa.enc.num_blocks); 
rsp.enc.RT = cell(1,pa.enc.num_blocks); 
rsp.enc.keyName = cell(1,pa.enc.num_blocks);
rsp.enc.accuracy = cell(1,pa.enc.num_blocks);
rsp.enc.time_beep_late_sec_array = cell(1,pa.enc.num_blocks);
for b = 1: 1:pa.enc.num_blocks
    rsp.enc.response{b} = [rsp.enc.response{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    rsp.enc.RT{b} = [rsp.enc.RT{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    rsp.enc.keyName{b} = [rsp.enc.keyName{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    rsp.enc.accuracy{b} = [rsp.enc.accuracy{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    rsp.enc.time_beep_late_sec_array{b} = [rsp.enc.time_beep_late_sec_array{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    
    seq.enc.trial_start_time{b} = [seq.enc.trial_start_time{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    seq.enc.stimOnset{b} = [seq.enc.stimOnset{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
    seq.enc.real_trial_onset_relative_to_block_start_time{b} = [seq.enc.real_trial_onset_relative_to_block_start_time{b}; repelem(nan,1,seq.enc.num_images_per_block(b))];
end

% loop_index=0;
% while loop_index ~= 1
indexPressed = 0;
while indexPressed == 0
    [keyIsDown, pressedSecs, keyCodes] = KbCheck(-1);
    if mod(str2double(pa.subjid),4) == 0 %(K)indoor/male (L) outdoor/female
        instructions = 'For this task, you will view some images.\n Please pay close attention to these images, \nas you will be asked questions about them later. \nYou will see each image repeats a couple of times. \n\n Some of these images are scenes, others are faces. \nPlease press the (K) button if the image is an indoor scene or a male face. \nPlease press the (L) button if the image in an outdoor scene or a female face. \nNote that you only need to press the button ONCE for the same image. \n\n Please let the experimenter know if you have any questions. \nNow press either (K) or (L) button to proceed. ';
        key_press_label = ["indoor","male";"outdoor","female"];
    elseif mod(str2double(pa.subjid),4) == 1%(K)indoor/female (L) outdoor/male
        instructions = 'For this task, you will view some images.\n Please pay close attention to these images, \nas you will be asked questions about them later. \nYou will see each image repeats a couple of times. \n\n Some of these images are scenes, others are faces. \nPlease press the (K) button if the image is an indoor scene or a female face. \nPlease press the (L) button if the image in an outdoor scene or a male face. \nNote that you only need to press the button ONCE for the same image. \n\n Please let the experimenter know if you have any questions. \nNow press either (K) or (L) button to proceed. ';
        key_press_label = ["indoor","female";"outdoor","male"];
    elseif mod(str2double(pa.subjid),4) == 2%(K)outdoor/male (L) indoor/female
        instructions = 'For this task, you will view some images.\n Please pay close attention to these images, \nas you will be asked questions about them later. \nYou will see each image repeats a couple of times. \n\n Some of these images are scenes, others are faces. \nPlease press the (K) button if the image is an outdoor scene or a male face. \nPlease press the (L) button if the image in an indoor scene or a female face. \nNote that you only need to press the button ONCE for the same image. \n\n Please let the experimenter know if you have any questions. \nNow press either (K) or (L) button to proceed. ';
        key_press_label = ["outdoor","male";"indoor","female"];
    elseif mod(str2double(pa.subjid),4) == 3%(K)outdoor/female (L) indoor/male
        instructions = 'For this task, you will view some images.\n Please pay close attention to these images, \nas you will be asked questions about them later. \nYou will see each image repeats a couple of times. \n\n Some of these images are scenes, others are faces. \nPlease press the (K) button if the image is an outdoor scene or a female face. \nPlease press the (L) button if the image in an indoor scene or a male face. \nNote that you only need to press the button ONCE for the same image. \n\n Please let the experimenter know if you have any questions. \nNow press either (K) or (L) button to proceed. .';
        key_press_label = ["outdoor","female";"indoor","male"];
    end
    key_press_button = ['k';'l'];
    
    Screen('TextSize', wPtr, 60);
    DrawFormattedText(wPtr, instructions, 'center', 'center',[255,255,255]);
    Screen('Flip',wPtr);
    if keyIsDown == 1 && KbName(find(keyCodes))~='n' && KbName(find(keyCodes))~='m'
        indexPressed = 1;
%             if KbName(find(keyCodes)) ~='n' && KbName(find(keyCodes))~='m'
%                 loop_index = 1;
%             elseif KbName(find(keyCodes)) ~='n' || KbName(find(keyCodes))~='m'
%                 loop_index = 0;
%             end
    end
end
Screen('Flip',wPtr);
% end
    
KbQueueCreate();%%make cue

% Start getting data from Shimmer
shimmer_ECG.start;
for block = 1:pa.enc.num_blocks
    disp('*********************************************');
    disp(['Start of Block ' num2str(block) ])
    disp('*********************************************');

    block_category = char(seq.enc.block_category_sequence(block));
    images_this_block = seq.enc.block_stim(block); 
    label_this_block = seq.enc.block_stim_label(block);
    in_block_ITI_this_block = seq.enc.in_block_ITI(block);
    cross_block_ITI = seq.enc.cross_block_ITI_sequence(block);
    cardiac_timing_this_block = seq.enc.cardiac_timing_block_sequence(block);
    
    % instruction
    % start instuctions

    %show block category prompt
    indexPressed = 0;
    while indexPressed == 0
        [keyIsDown, pressedSecs, keyCodes] = KbCheck(-1);

        if block_category=="face"
            if mod(str2double(pa.subjid),4) == 0  || mod(str2double(pa.subjid),4) == 2 %(K)male (L) female 
                instructions = ['You will see ',block_category,'s in this block. \n\nPlease press the (K) button if the image is a male face. \nPlease press the (L) button if the image a female face.\n\n To start this block, press (K) or (L) button.'];
            else %(K)female (L) female 
                instructions = ['You will see ',block_category,'s in this block. \n\nPlease press the (K) button if the image is a female face. \nPlease press the (L) button if the image a male face.\n\n To start this block, press (K) or (L) button.'];
            end
        elseif block_category=="scene"
            if mod(str2double(pa.subjid),4) == 0  || mod(str2double(pa.subjid),4) == 1 %(K)indoor (L) outdoor 
                instructions = ['You will see ',block_category,'s in this block. \n\nPlease press the (K) button if the image is an indoor scene. \nPlease press the (L) button if the image an outdoor scene.\n\n To start this block, press (K) or (L) button.'];
            else %(K)outdoor (L) indoor 
                instructions = ['You will see ',block_category,'s in this block. \n\nPlease press the (K) button if the image is an outdoor scene. \nPlease press the (L) button if the image an indoor scene.\n\n To start this block, press (K) or (L) button.'];
            end
        end

        Screen('TextSize', wPtr, 60);
        DrawFormattedText(wPtr, instructions, 'center', 'center',[255,255,255]);
        Screen('TextSize', wPtr, 60);
        Screen('Flip',wPtr);
        if keyIsDown == 1 
            indexPressed = 1;
            Screen('Flip',wPtr);
        end           
    end
    Screen('Flip',wPtr);
    
    seq.enc.block_start_time(block) = GetSecs;    
    
%     %cross block ITI
%     fixation_absolute_onset_time=seq.enc.block_start_time(block);
%     fixation_absolute_offset_time=seq.enc.block_start_time(block) + cross_block_ITI;
        
    %draw fixation
    fixation = '+';
    Screen('TextSize', wPtr, 80);
    DrawFormattedText(wPtr, fixation, 'center', 'center');
    Screen('Flip',wPtr);
    WaitSecs(1);
%     [~, stimOnset] = Screen('Flip',wPtr,image_absolute_onset_time);
%     [~, stimOffset] = Screen('Flip',wPtr, image_absolute_offset_time);
%     WaitSecs('UntilTime',fixation_absolute_offset_time);
%     Screen('Flip',wPtr);

    KbQueueStart();%%start listening
    
    %trial level
    for trial = 1:size(images_this_block{1},3)
        disp(['trial ' num2str(trial) ]);        
        seq.enc.trial_start_time{block}(trial) = GetSecs;
 
        %image
        image_to_show=images_this_block{1}(:,:,trial);
        [imageHeight, imageWidth, colorChannels] = size(image_to_show);
        pointer_to_offscreen_window_with_image_on_it = Screen('MakeTexture',wPtr,image_to_show);
        
        %ITI this trial
        in_block_ITI_this_trial = in_block_ITI_this_block{1}(trial);
        
        trial_label = label_this_block{1}(trial);
        
        % Initialize timers and other variables

        Npeaks_counted = 0;
        Nbeeps_delivered = 0;
        Npeaks_found_that_do_not_count = 0;
        FIRST_READ_HAS_OCCURRED = false;
        FOUND_PEAKS_ON_THIS_TRIAL = false;

        X_time_total_sec=[];
        Y_ECG_total=[];
        all_peaks_X=[];
        all_peaks_Y=[];
        
        X_time_beep_sec_array = nan(1,pa.enc.image_repetition);
        X_time_beep_late_sec_array = nan(1,pa.enc.image_repetition);
        X_time_detected_peak_sec_array = nan(1,pa.enc.image_repetition);
        Y_ECG_detected_peak_mV_array = nan(1,pa.enc.image_repetition);
        Y_peak_prominance_array = nan(1,pa.enc.image_repetition);

        tic_at_last_sample = tic;
        % Reset sampling bufer for first heartbeat
        sampling_buffer_first_HB_sec = 0.1;

        hfig = figure;

        KbQueueFlush();%%removes all previous keyboard presses
        
        %set priority for screen to optimize display timing once peaks are found
%         Priority(MaxPriority(wPtr));
        Priority(topPriorityLevel);

        image_delivered = false;
                                    
        while Npeaks_counted < pa.enc.image_repetition + in_block_ITI_this_trial
%                 disp(['Npeaks_counted: ' num2str(Npeaks_counted) ]);
                    
            % find R peak to determine picture onset
            if toc(tic_at_last_sample) >= (sampling_period_check_Shimmer_sec + sampling_buffer_first_HB_sec)         
                
                % Get time when data arrive
                time_current_packet_arrival_pre_GetSecs = GetSecs();

                % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
                [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer_ECG.getdata('c');

                time_current_packet_arrival_post_GetSecs = GetSecs();

                NnewData = size(newData,1);
                     
                % Process and display data
                if (NnewData >= 2)
                    
                    if( ~FIRST_READ_HAS_OCCURRED  )
%                         disp('first reading!');
                        k_ECG_Signal_HBD = find(strcmp(string_ECG_signal_HBD, signalNameArray));

                        k_time = find(strcmp('Time Stamp', signalNameArray));
                        %k_ECG_LL_RA = find(strcmp('ECG LL-RA', signalNameArray));
                        %k_ECG_Signal_HBD = find(strcmp('ECG LA-RA', signalNameArray));

                        FIRST_READ_HAS_OCCURRED = true;

                        % Reset this to zero now that it has been used
                        sampling_buffer_first_HB_sec = 0;
                    end

                    X_time_current_packet_sec               = newData(:, k_time) / 1e3;
                    X_time_total_sec(end+1 : end+NnewData)  = newData(:, k_time) / 1e3;                                
                    Y_ECG_total(end+1 : end+NnewData)         = newData(:, k_ECG_Signal_HBD);

                    Y_ECG_total_adj = ECG_adjust_baseline_spline( X_time_total_sec, Y_ECG_total, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );

%                     % Real-time plot of trial 
%                     hold on;
%                     plot(X_time_total_sec(end-NnewData+1:end), Y_ECG_total_adj(end-NnewData+1:end), '-k');
%                     drawnow();

                    %------------------------------------------------------
                    % Detection of all R spikes for this trial
                    % so far (between 0 and Ntones)

                    % Ensure there are enough data for peak detection                            
                    time_of_total_recording_sec = X_time_total_sec(end) - X_time_total_sec(1);

                    % If there is at least one peak already
                    if( length(X_time_detected_peak_sec_array) >= 1 )
                        time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_detected_peak_sec_array(end);
                    else
                        % No actual peaks yet so get total duration of recording
                        time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_total_sec(1);
                    end
    %                         fprintf('\nTime since last R spike = %0.4f sec \n', time_since_last_Rspike_sec);

                    if( time_of_total_recording_sec > Minimum_RR_Interval_sec )
                        % Peak detection
                        [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                            findpeaks(Y_ECG_total_adj, X_time_total_sec, ...
                            'MinPeakDistance', Minimum_RR_Interval_sec, ...
                            'MinPeakProminence', Minimum_R_Prominence_mV);
                        
%                         % Update plot
%                         hold on;
%                         plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);
%                         drawnow()

                        Npeaks_found = length(peak_X_array);
                        
                        all_peaks_X(end+1:end+Npeaks_found) = peak_X_array;
                        all_peaks_Y(end+1:end+Npeaks_found) = peak_Y_array;
                        
                        % First peak disvorey has to be ONE peak (not more than that). If there
                        % are 2 peaks on first search, then the first one doesn't count
                        if( ~FOUND_PEAKS_ON_THIS_TRIAL && Npeaks_found > 0 )                                        
                            Npeaks_found_that_do_not_count = Npeaks_found - 1;
                            FOUND_PEAKS_ON_THIS_TRIAL = true;
                        end

                        %------------------------------------------------------
                        % If (1) it has been enough time since last R spike and (2) there is a new peak that hasn't previously been counted
                        %if( (toc(tic_at_last_Rspike) >= Minimum_RR_Interval_sec) && Npeaks_found > Npeaks_counted )
                        if( Npeaks_found - Npeaks_found_that_do_not_count > Npeaks_counted )
                            
                            Npeaks_counted = Npeaks_counted + 1;

                            if Nbeeps_delivered < pa.enc.image_repetition
            
                                % Time from packet start to R spike
                                time_from_packet_start_to_peak_sec = peak_X_array(end) - X_time_current_packet_sec(1);
                                time_beep_intended_sec = time_current_packet_arrival_pre_GetSecs + time_from_packet_start_to_peak_sec;

                                % If it LOOKS like we are waiting TOO LONG to beep, then just set it to the
                                % max waiting duration (the delay itself)
                                time_to_beep_sec = time_beep_intended_sec - GetSecs();
                                if( time_to_beep_sec > cardiac_timing_this_block)
                                    time_beep_intended_sec = GetSecs();
                                end
                                
                                % image_onset = R spike + cardiac timing offset
                                image_absolute_onset_time = time_beep_intended_sec + cardiac_timing_this_block;
                            
                                % image offset = onset time + image duration
                                image_absolute_offset_time = image_absolute_onset_time + pa.enc.image_duration;

%                                      %minor fading in effect
%                                     for thisContrast = [20,50,80]
%                                         Screen('DrawTexture',wPtr,pointer_to_offscreen_window_with_image_on_it, [], [], 0, [], thisContrast);
%                                         DrawFormattedText(wPtr, fixation, 'center', 'center');
%                                         Screen('Flip',wPtr);
%                                     end

                                %draw image and fixation
                                Screen('DrawTexture',wPtr,pointer_to_offscreen_window_with_image_on_it);
                                DrawFormattedText(wPtr, fixation, 'center', 'center');

                                %onset when timing right and haven't
                                %presented yet
                                if (GetSecs() < image_absolute_onset_time) && (~image_delivered) 
% %                                     (GetSecs() >= image_absolute_onset_time) && (GetSecs() < image_absolute_offset_time) && (~image_delivered) 
                                    [~, stimOnset] = Screen('Flip',wPtr,image_absolute_onset_time);
                                    seq.enc.stimOnset{block}(trial) = stimOnset;
                                    seq.enc.real_trial_onset_relative_to_block_start_time{block}(trial) = stimOnset-seq.enc.block_start_time(block);

                                    Nbeeps_delivered = Nbeeps_delivered+1;
                                    image_delivered = true;

                                    duration_beep_late = stimOnset - image_absolute_onset_time;
                                    beep_delay = stimOnset - time_beep_intended_sec;

                                    % Save time of beep
                                    X_time_beep_sec_array(Nbeeps_delivered)            = stimOnset;
                                    X_time_beep_late_sec_array(Nbeeps_delivered)       = duration_beep_late;  
                                    X_time_beep_delay_sec_array(Nbeeps_delivered)      = beep_delay;  
                                    X_time_detected_peak_sec_array(Nbeeps_delivered)   = peak_X_array(end);                                        
                                    Y_ECG_detected_peak_mV_array(Nbeeps_delivered)     = peak_Y_array(end);                                        
                                    Y_peak_prominance_array(Nbeeps_delivered)          = peak_prom_array(end);
                                    
%                                     % Update plot
%                                     hold on;
%                                     plot(X_time_detected_peak_sec_array(end) + rsp.enc.time_beep_late_sec_array(end), peak_Y_array(end), 'sb', 'LineWidth', 3);
%                                     % Update the plot
%                                     drawnow()
                                    
                                end
                                
                                %offset when timing right and images
                                %already presented
                                if (GetSecs() < image_absolute_offset_time) && (image_delivered) 
                                    [~, stimOffset] = Screen('Flip',wPtr, image_absolute_offset_time);
                                    Screen('Flip',wPtr);
                                    image_delivered = false;
                                    
                                    DrawFormattedText(wPtr, fixation, 'center', 'center');
                                    Screen('Flip',wPtr);
                                end

                            else
                                %draw + during in block ITI time after a image is repeated 4 times
                                DrawFormattedText(wPtr, fixation, 'center', 'center');
                                Screen('Flip',wPtr);

%                                 % Update plot
%                                 hold on;
%                                 plot(peak_X_array(end)+0.05, peak_Y_array(end), 'sc', 'LineWidth', 3);
%                                 % Update the plot
%                                 drawnow()

                            end

%                             %always show fixation after image is off screen
%                             DrawFormattedText(wPtr, fixation, 'center', 'center');
%                             Screen('Flip',wPtr);
% 
%                             %release priority on screen once display of image is off
%                             Priority(0);
                        end  
                        
                        % Write data to file
                        time_ECG_sec = newData(:, k_time) / 1e3;
                        Y_ECG_mV = newData(:, k_ECG_Signal_HBD);
                        out = [time_ECG_sec, Y_ECG_mV];
                        % write data to end of file
                        dlmwrite(sprintf('%sencoding-ECG.csv',output_dir),out,'delimiter',',','precision',6,'-append');

%                         X_peaks_adj = nan(1,length(X_time_total_sec));
%                         if ~isempty(all_peaks_X)
%                             X_peaks_adj(X_time_total_sec==all_peaks_X') = all_peaks_X'; 
%                         end
%                         out = [X_time_total_sec; Y_ECG_total_adj;X_peaks_adj];
                        out = [X_time_total_sec; Y_ECG_total_adj];
                        % write data to end of file
                        dlmwrite(sprintf('%sencoding-ECG-Baseline_Adjusted.csv',output_dir),out,'delimiter',',','-append');

                        IBI_sec_array = diff( [nan peak_X_array] );
                        out = [peak_X_array; peak_Y_array; peak_prom_array; IBI_sec_array]';
                        % write data to end of file
                        dlmwrite(sprintf('%sencoding-ECG-Rspikes.csv',output_dir),out,'delimiter',',','-append');

                    end
                    
                end
                tic_at_last_sample = tic;
                
            end
 
        end

        %release priority
        Priority(0);
        
        close(hfig);
        
        hfig2 = figure;
        plot(X_time_total_sec, Y_ECG_total_adj, '-k');
        hold on;
        plot(all_peaks_X, all_peaks_Y, 'or', 'LineWidth', 3);
        hold on;
        plot(X_time_detected_peak_sec_array + X_time_beep_delay_sec_array, repmat(peak_Y_array(end),1,length(X_time_detected_peak_sec_array)), 'sb', 'LineWidth', 3);
        title(sprintf('Subj %s, Block %d, Trial %d, Delay=%0.4f s', pa.subjid, block, trial, cardiac_timing_this_block));
        xlabel('Time (sec)');
        ylabel('Adjusted ECG (mv)');
        grid('on');
        print(hfig2, sprintf('%sECG-Block_%02d_Trial_%02d.png', output_dir, block, trial), '-dpng');   
        hgsave(hfig2, sprintf('%sECG-Block_%02d_Trial_%02d.fig', output_dir, block, trial));
        drawnow()
        close(hfig2);
        
        [pressed, firstpress] = KbQueueCheck(); %check response
        if isnan(rsp.enc.response{block}(trial))
            %get earliest time stamp of active key pressed
            active_keys_pressed = firstpress(active_keys);
            active_keys_pressed(active_keys_pressed == 0 ) = nan;
            [first_press_time_relative_to_trial_start,first_press_index] = min(active_keys_pressed);
            if pressed && (first_press_time_relative_to_trial_start > 0)
                rsp.enc.response{block}(trial) = 1;
                rsp.enc.RT{block}(trial)       = first_press_time_relative_to_trial_start;
                key_pressed                    = active_keys(first_press_index);
                rsp.enc.keyName{block}(trial)  = KbName(key_pressed);
                if convertCharsToStrings(block_category) == "face"
                    if find(KbName(key_pressed)==key_press_button) == find(trial_label==key_press_label(:,2))
                        rsp.enc.accuracy{block}(trial) = 1;
                    else
                        rsp.enc.accuracy{block}(trial) = 0;
                    end
                elseif convertCharsToStrings(block_category) == "scene"
                    if find(KbName(key_pressed)==key_press_button) == find(trial_label==key_press_label(:,1))
                        rsp.enc.accuracy{block}(trial) = 1;
                    else
                        rsp.enc.accuracy{block}(trial) = 0;
                    end
                end 
                disp(rsp.enc.accuracy{block}(trial))
            end
        end
        
        KbQueueFlush();
%         KbQueueStop();
        
        disp(['pressed: ', rsp.enc.keyName{block}(trial), ' acc: ', num2str(rsp.enc.accuracy{block}(trial))]);
            
        rsp.enc.time_beep_late_sec_array{block}(end+1:end+pa.enc.image_repetition)=X_time_beep_late_sec_array;
    end
    
    Screen('TextSize', wPtr, 80);
    DrawFormattedText(wPtr, fixation, 'center', 'center');
    Screen('Flip',wPtr);
        
%     %save file
%     save(sprintf('%ssequence.mat',output_dir), 'seq')
% %     save(sprintf('%sstimuli.mat',output_dir), 'stim')
%     save(sprintf('%sresponse.mat',output_dir), 'rsp')
    
    disp('*********************************************');
    disp(['End of Block ' num2str(block) ])
    disp('*********************************************');
    
end

%save csv
trial_num = [1:pa.enc.num_total_trials]';
block_num = []; cardiac_cycle=[]; cardiac_timing=[]; trial_label=[]; 
response=[]; RT=[]; trial_accuracy=[];key=[];
stim_onset=[]; trial_start_time=[]; real_trial_onset_relative_to_block_start_time=[];
for b = 1: 1:pa.enc.num_blocks
    block_num = [block_num;repelem(b,seq.enc.num_images_per_block(b),1)];
    cardiac_cycle = [cardiac_cycle;repelem(seq.enc.cardiac_cycle_block_sequence(b),seq.enc.num_images_per_block(b),1)];
    cardiac_timing = [cardiac_timing;repelem(seq.enc.cardiac_timing_block_sequence(b),seq.enc.num_images_per_block(b),1)];
    trial_label = [trial_label;seq.enc.block_stim_label{b}'];
    response = [response;rsp.enc.response{b}'];
    RT = [RT;rsp.enc.RT{b}'];
    trial_accuracy = [trial_accuracy;rsp.enc.accuracy{b}'];
    key = [key;rsp.enc.keyName{b}'];
    stim_onset = [stim_onset;seq.enc.stimOnset{b}'];
    trial_start_time = [trial_start_time;seq.enc.trial_start_time{b}'];
    real_trial_onset_relative_to_block_start_time = [real_trial_onset_relative_to_block_start_time;seq.enc.real_trial_onset_relative_to_block_start_time{b}'];
end
T = table(trial_num,block_num,cardiac_cycle,cardiac_timing,trial_label,RT,trial_accuracy,key,response,stim_onset,trial_start_time,real_trial_onset_relative_to_block_start_time);
writetable(T,sprintf('%ssubj-%s_Encoding_Response_by_Trial.csv', output_dir, pa.subjid));

% end instuctions
keyIsDown = 0;
indexPressed = 0;
while keyIsDown == 0 && indexPressed == 0
    [keyIsDown, pressedSecs, keyCodes] = KbCheck(-1);

    instructions = 'You have completed the first part of the experiment.\n Please take a short break.\n\n Press (K) or (L) key to proceed to the next part.';
      
    Screen('TextSize', wPtr, 60);
    DrawFormattedText(wPtr, instructions, 'center', 'center',[255,255,255]);
    
    Screen('Flip',wPtr);
    if keyIsDown == 1
        indexPressed = 1;
        keyIsDown = 0;
    end           
end
 
if strcmpi(pa.shimmer, 'true')
    shimmer_ECG.stop;
    shimmer_ECG.disconnect;
    fclose(FILE_OUT_ECG);
    fclose(FILE_OUT_ECG_Adj);
    fclose(FILE_OUT_ECG_RSPIKES);
end 
 
%save file
save(sprintf('%ssequence.mat',output_dir), 'seq')
save(sprintf('%sstimuli.mat',output_dir), 'stim')
save(sprintf('%sresponse.mat',output_dir), 'rsp')
save(sprintf('%sparameter.mat',output_dir), 'pa')

Screen('Flip',wPtr);
sca;

InteroceptiveReinstatement_Recognition;
rsp.enc.time_beep_late_sec_array
