function []=conn_batch_execute(varargin)
%%vzvz
%Arguments;
% 1: subject folder name
% 2: string defining run/batch number
%
if(length(varargin)>0)
    Sub_IDs=varargin{1};
end

if(length(varargin)>1)
    iter_no=varargin{2};
end
iter_no = [num2str(1),'/','Node_',num2str(0)]
%% Change data, run and results directories
% run_dir: a new dir where Conn saves intermediate files that can be used
% by SPM later
% data_dir: directory where T1 and bold files have to be located
% 

% Go through subejcts in a  folders,check if the results folder exist. If
% exists, then go to next subject. If not exist, then start preprocess.

ids_array = cell2mat(Sub_IDs(1))
id = num2str(ids_array)

run_dir = strcat('[Insert filepath to run_dir]',id,'/')

data_dir=['[Insert filepath to data_dir]']
roifile=[cd,'/gordon333MNI.nii'];

new_b=create_batch_data(iter_no,run_dir,data_dir,Sub_IDs);
DataDir=new_b;

% ******************************************************************************************************
%% STEP 0: INITIALIZING/COPYING NECESSARY DATA
% ******************************************************************************************************

if(isnumeric(iter_no))
    iter=num2str(iter_no);
else
    iter=iter_no;
end

% ******************************************************************************************************
%% STEP 1: FUNCTIONAL PROCESSING PIPELINE
% ******************************************************************************************************


DOPARALLEL=false;   % set to true/false to run in parallel or locally


%% FIND functional/structural files

NSUBJECTS=0;
cwd=pwd;
FUNCTIONAL_FILE={};
STRUCTURAL_FILE={};

data=dir(DataDir);% Refers to NKI rockland data set

% saving summary data
jid= str2num(getenv('SLURM_JOB_ID'));
SumLog = fopen('SummaryLog.txt','w');

NSUBS = 1

fprintf(SumLog,'Detailed log can be found in file: \n');
fprintf(SumLog,['output.',num2str(jid),'.out \n']);

fprintf(SumLog,'STARTED RUNNING CONN BATCH \n');
fprintf(SumLog,'Attempting to read %d subjects \n',NSUBS);

for n=1:NSUBS,
    filename=[DataDir,num2str(cell2mat(Sub_IDs(n))),'/'];
    bold=conn_dir(fullfile(filename,['60','*','bold','*','.nii']));
    t1=conn_dir(fullfile(filename,['60','*','T','*','.nii']));
    if isfile(bold)
        fprintf('The bold file has been found');
    end
    if isfile(t1)
        fprintf('The t1 file has been found');
    end
    if(~isempty(bold)&&~isempty(t1))
        tFUNCTIONAL_FILE=cellstr(conn_dir(fullfile(filename,['60','*','bold','*','.nii'])));
        tSTRUCTURAL_FILE=cellstr(conn_dir(fullfile(filename,['60','*','T','*','.nii'])));
        FUNCTIONAL_FILE=[FUNCTIONAL_FILE;tFUNCTIONAL_FILE(:)];
        STRUCTURAL_FILE=[STRUCTURAL_FILE;tSTRUCTURAL_FILE(:)];
        fprintf('\t Subject: %s \n',num2str(cell2mat((Sub_IDs(n)))));
        fprintf('\t \t Functional file: %s \n',char(tFUNCTIONAL_FILE));
        fprintf('\t \t Structural file: %s \n',char(tSTRUCTURAL_FILE));
        
        fprintf(SumLog,'\t Subject: %s \n',num2str(cell2mat((Sub_IDs(n)))));
        fprintf(SumLog,'\t \t Functional file: %s \n',char(tFUNCTIONAL_FILE));
        fprintf(SumLog,'\t \t Structural file: %s \n',char(tSTRUCTURAL_FILE));
    else
        nf= [isempty(bold),isempty(t1)];Fnd={'Found','Not Found'};
        fprintf(SumLog,'\t Subject: %s \n',num2str(cell2mat((Sub_IDs(n)))));
        fprintf(SumLog,'\t \t NOT FOUND Functional file: %s Structural file: %s \n',char(Fnd(nf(1)+1)),char(Fnd(nf(2)+1)));
    end
end
fprintf('The NSUBJECTS is ',NSUBJECTS);
NSUBJECTS = 1;
if ~NSUBJECTS, NSUBJECTS=length(STRUCTURAL_FILE); end
if rem(length(FUNCTIONAL_FILE),NSUBJECTS),error('mismatch number of functional files %n', length(FUNCTIONAL_FILE));end
if rem(length(STRUCTURAL_FILE),NSUBJECTS),error('mismatch number of anatomical files %n', length(FUNCTIONAL_FILE));end
nsessions=length(FUNCTIONAL_FILE)/NSUBJECTS;
FUNCTIONAL_FILE=reshape(FUNCTIONAL_FILE,[NSUBJECTS,nsessions]);
STRUCTURAL_FILE={STRUCTURAL_FILE{1:NSUBJECTS}};
disp([num2str(size(FUNCTIONAL_FILE,1)),' subjects']);
disp([num2str(size(FUNCTIONAL_FILE,2)),' sessions']);

disp([num2str(size(STRUCTURAL_FILE,1)),' subjects']);
disp([num2str(size(STRUCTURAL_FILE,2)),' sessions']);

TR=3; % Repetition time = 3s PLEASE NOTE THIS!!!!!!!!!!!!!!!!!!!!!!!
aTime=tic;

%% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
%% Prepares batch structure
clear batch;
batch.filename=fullfile(cwd,'conn_NKI.mat');            % New conn_*.mat experiment name

%% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;

batch.Setup.preprocessing.steps={'functional_label_as_original','functional_realign&unwarp','functional_center','functional_art',...
    'functional_segment&normalize_direct','functional_label_as_mnispace','structural_center','structural_segment&normalize',...
    'functional_smooth','functional_label_as_smoothed'};% skip slice timing

batch.Setup.RT=TR;                                        % TR (seconds)
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);       % Point to functional volumes for each subject/session
for nsub=1:NSUBJECTS,
    for nses=1:nsessions,
        batch.Setup.functionals{nsub}{nses}{1}=FUNCTIONAL_FILE{nsub,nses};
    end
end % note: the data from each subject is defined by three sessions and one single (4d) file per session

batch.Setup.structurals=STRUCTURAL_FILE;                % Point to anatomical volumes for each subject
nconditions=nsessions;                                  % treats each session as a different condition (comment the following three lines and lines 84-86 below if you do not wish to analyze between-session differences)
if nconditions==1
    batch.Setup.conditions.names={'rest'};
    for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
else
    batch.Setup.conditions.names=[{'rest'}, arrayfun(@(n)sprintf('Session%d',n),1:nconditions,'uni',0)];
    for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
    for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=1:nsessions,  batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=[];batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=[]; end;end;end
    for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=ncond,        batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=inf;end;end;end % session-specific conditions
end

connDir='[Insert filepath to the directory containing the CONN resources]';
 
for ii=1:333
    Gordon_rois2.files{ii}=[connDir, 'Gordon_rois/Gordon_ROI',num2str(ii),'.nii'];
end

batch.Setup.rois=Gordon_rois2;
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';
batch.Setup.preprocessing.art_thresholds(1)= 3; % threshold value for global-signal (z-value; default 5) % PLEASE NOTE THIS!!!!!!!!!!!!!!!!!!!!!!!
batch.Setup.preprocessing.art_thresholds(2)= 0.5; % threshold value for subject-motion (mm; default .9) % PLEASE NOTE THIS!!!!!!!!!!!!!!!!!!!!!!!

PreprocTime=0;

%% DENOISING step
% CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options
batch.Denoising.filter=[0.01, 0.1];                 % frequency filter (band-pass values, in Hz) % PLEASE NOTE THIS!!!!!!!!!!!!!!!!!!!!!!!
batch.Denoising.done=1;
batch.Denoising.overwrite='Yes';

%% FIRST-LEVEL ANALYSIS step
% CONN Analysis                                     % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';

%% Run all analyses
tic;
conn_batch(batch);
AnalysisTime=toc;

TotalTime=toc(aTime);

save('Timings.mat','AnalysisTime','PreprocTime','TotalTime');

fprintf(SumLog,'Started %d subjects \n',length(Sub_IDs));

fprintf('FINISHED RUNNING CONN BATCH \n');
fprintf(SumLog,'FINISHED RUNNING CONN BATCH \n');

matrix_file = strcat(run_dir, '[Insert filepath to the directory containing the results]/resultsROI_Condition001.mat');

if isfile(matrix_file)
    fprintf('THE PREPROCESS IS SUCCESSFUL AND CONNECTIVITY MATRIX EXISTS');
    status = 'Successful';
else
    fprintf('THE PREPREOCESS IS FAILED');
    status = 'Failed';
end


fclose(SumLog);

end