function run_boost_Israel()

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================ by Rotem Botvinik November-December 2014 ===============
% =================== Modified by Tom Salomon, April 2015 =================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% Runs the cue-approach task

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'binary_ranking'
% % %   'binary_ranking_demo'
% % %   'sort_binary_ranking'
% % %   'trainingDemo_Israel'
% % %   'training_Israel' (if NO EYETRACKING)
% % %   'organizeProbe_Israel'
% % %   'probeDemo_Israel'
% % %   'probe_Israel'
% % %   'probeResolve_Israel'
% % %   'recognitionNewOld_Israel'
% % %   'recognitionGoNoGo_Israel'

% % %   --- Other codes: ---
% % %  'CenterText'

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FOLDERS REQUIRED TO RUN PROPERLY: ------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'Misc': a folder with the audio file.
% % %   'Onset_files': a folder with the onset files for the training and
% % %    for the probe.
% % %   'Output': the folder for the output files- results.
% % %   'stim': with the image files of all the stimuli for the cue-approach
% % %    task (old stimuli).
% % %   'stim/recognitionNew': with the image files of the new stimuli
% % %   (stimuli that are not included in the cue-approach tasks, only in the
% % %   recognitionNewOld task, as new stimuli).

tic

% =========================================================================
%% Get input args and check if input is ok
% =========================================================================

% %---dummy info for debugging purposes --------
% subjectID =  'BM_001';

% timestamp
c = clock;
hr = sprintf('%02d', c(4));
min = sprintf('%02d', c(5));
timestamp=[date,'_',hr,'h',min,'m'];

% Increase the size of the UI fonts
old_UI_font_size=get(0,'defaultUicontrolFontSize');
set(0,'defaultUicontrolFontSize', 14);

% essential for randomization
rng('shuffle');

% input checkers
subjectID = input('Subject code: ','s');
[subjectID_num,okID]=str2num(subjectID(end-2:end));
while okID==0
    disp('ERROR: Subject code must contain 3 characters numeric ending, e.g "BMI_bf_001". Please try again.');
    subjectID = input('Subject code:','s');
    [subjectID_num,okID]=str2num(subjectID(end-2:end));
end

% Assign order
% --------------------------
% give order value of '1' or '2' for subjects with odd or even ID, respectively
if mod(subjectID_num,2) == 1 % subject code is odd
    order = 1;
else % subject code is even
    order = 2;
end

sessionNum=1;

% set the computer and path
% --------------------------
test_comp=0; % 1 MRI, 0 if testooom
mainPath = pwd; % Change if you don't run from the experimnet folder - not recommended.
outputPath = [mainPath '/Output'];

% open a txt file for crashing logs
fid_crash = fopen([outputPath '/' subjectID '_crashingLogs' num2str(sessionNum) '_' timestamp '.txt'], 'a');

% =========================================================================

%% Ask if you want Eye Tracker
% =========================================================================
%use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Eye Tracker', 'Yes', 'No', 'Yes');
%if strcmp(ask_if_want_eyetracker, 'Yes')
%    use_eyetracker = 1;
%end

% =========================================================================
%% Personal Details
% =========================================================================
personal_details(subjectID, order, outputPath, sessionNum)

% =========================================================================
%% Part 1 - Binary Ranking (including demo)
% =========================================================================

Crashesd_binary_ranking_demo = 0;
keepTrying = 1;
while keepTrying < 10
    try
        binary_ranking_demo(subjectID,test_comp,mainPath);
        
        % Ask if subject wanted another demo
        % ----------------------------------
        demo_again = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        if strcmp(demo_again,'Yes')
            binary_ranking_demo(subjectID,test_comp,mainPath);
        end
        keepTrying = 10;
    catch
        sca;
        Crashesd_binary_ranking_demo = Crashesd_binary_ranking_demo + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - BINARY RANKING DEMO!');
    end
end
fprintf(fid_crash,'Binary ranking demo crashed:\t %d\n', Crashesd_binary_ranking_demo);
use_eyetracker = 0;
ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Eye Tracker', 'Yes', 'No', 'Yes');
if strcmp(ask_if_want_eyetracker, 'Yes')
    use_eyetracker = 1;
end

Crashesd_binary_ranking = 0;
keepTrying = 1;
while keepTrying < 10
    try
        binary_ranking_withMix(subjectID,test_comp,mainPath,use_eyetracker);
        keepTrying = 10;
    catch
        sca;
        Crashesd_binary_ranking = Crashesd_binary_ranking + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - BINARY RANKING!');
    end
end
fprintf(fid_crash,'Binary ranking crashed:\t %d\n', Crashesd_binary_ranking);

% Sort stimuli according to the binary ranking
% -------------------------------------------------
sort_binary_ranking(subjectID,order,outputPath);

% Break    
questdlg('Thank you. Press when you are ready to continue'...
    ,'Break','Continue','Continue');
WaitSecs(0.5);
    
% =========================================================================
%% Part 2 - Training (including demo)
% =========================================================================
% Set number of runs
% -------------------------------------------------

crashedDemoTraining = 0;
keepTrying = 1;
while keepTrying < 10
    try
        trainingDemo_Israel(subjectID,order,mainPath,test_comp);
        
        % Ask if subject wanted another demo
        % ----------------------------------
        demo_again = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        if strcmp(demo_again,'Yes')
            trainingDemo_Israel(subjectID,order,mainPath,test_comp);
        end
        keepTrying = 10;
    catch
        sca;
        crashedDemoTraining = crashedDemoTraining + 1;
        keepTrying = keepTrying + 1;
        if keepTrying==9
            skip_synctest=1;
        end
    end
end
Screen('Preference', 'SkipSyncTests', 0);
fprintf(fid_crash,'trainingDemo crashed:\t %d\n', crashedDemoTraining);

Ladder1IN = 750;
Ladder2IN = 750;

total_num_runs_training = 20;
%total_num_runs_training = 4; % for debugging
num_parts_of_training = 2; % we use eyetracking and want to calibrate and validate it once at the middle of training
parts = linspace(1,total_num_runs_training,num_parts_of_training+1);
runInd = ceil(parts(1:end-1));
last_run = [ceil(parts(2:end-1)-1) parts(end)];

for part_num = 1:num_parts_of_training
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Eye Tracker', 'Yes', 'No', 'Yes');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
    
    crashedTraining = 0;
    keepTrying = 1;
    
    while keepTrying < 10
        try
            training_Israel(subjectID,order,mainPath,test_comp,runInd(part_num),total_num_runs_training,Ladder1IN,Ladder2IN,use_eyetracker, last_run(part_num));
            keepTrying = 10;
        catch
            sca;
            crashedTraining = crashedTraining + 1;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - TRAINING!');
            if keepTrying==9
                skip_synctest=1;
            end
        end
    end
    Screen('Preference', 'SkipSyncTests', 0);
    
    fprintf(fid_crash,'training crashed:\t %d\n', crashedTraining);
    
    if part_num ~= num_parts_of_training
        questdlg('Time for a short break. Please call the experimenter.','','OK','OK');
    end
end
fprintf(fid_crash,'training crashed:\t %d\n', crashedTraining);

%==========================================================
%% Part 3 - Break
%==========================================================

questdlg('For the next step, please take the time to fill a few questionnaires. Please call the experimenter'... '
    ,'part 3','Continue','Continue');
WaitSecs(1);
% Mood scale
system(['/usr/local/bin/python2.7 ',pwd,'/mood_scale.py ' subjectID]);

% Fractal rankings
keepTrying = 0;
while keepTrying < 10
    try
        ScaleRanking(subjectID, 0,1);
        keepTrying = 10;
    catch
        keepTrying = keepTrying + 1;
        sca
        WaitSecs(0.1);
    end
end
WaitSecs(1);

% Break    
questdlg('Thank you. Press when you are ready to continue'...
    ,'Break','Continue','Continue');
WaitSecs(0.5);

%==========================================================
%% Part 4 - probe_demo & probe
%==========================================================

numBlocks = 2; % Define how many blocks for the probe session. Each block includes all comparisons, one time each.
numRunsPerBlock = 1; % Define the required number of runs per block

crashedDemoProbe = 0;
keepTrying = 1;
while keepTrying < 10
    try
        probeDemo_Israel(subjectID, order, mainPath, test_comp);
        
        % Ask if subject wanted another demo
        % ----------------------------------
        demo_again = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        if strcmp(demo_again,'Yes')
            probeDemo_Israel(subjectID, order, mainPath, test_comp);
        end
        keepTrying = 10;
    catch
        sca;
        crashedDemoProbe = crashedDemoProbe + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - PROBE DEMO!');
    end
end
fprintf(fid_crash,'probe demo crashed:\t %d\n', crashedDemoProbe);


% Run blocks. Before each block, stimuli need to be organized in a list and
% divided to the required number of runs
for ind = 1:numBlocks
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Eye Tracker', 'Yes', 'No', 'Yes');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
    
    block = (sessionNum-1)*numBlocks + ind;
    % Organize the stimuli for the probe
    % ===================================
    [trialsPerRun] = organizeProbe_Israel(subjectID, order, mainPath, block, numRunsPerBlock);
    crashedProbe=0;
    for numRun = 1:numRunsPerBlock
        keepTrying = 1;
        while keepTrying < 10
            try
                probe_Israel(subjectID, order, mainPath, test_comp, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun, use_eyetracker);
                keepTrying=10;
            catch
                sca;
                keepTrying = keepTrying+1;
                crashedProbe = crashedProbe + 1;
                disp('CODE HAD CRASHED - PROBE!');
            end
        end
    end % end for numRun = 1:numRunsPerBlock
end % end for block = 1:numBlocks
fprintf(fid_crash,'probe crashed:\t %d\n', crashedProbe);

% Break    
questdlg('Thank you. You can continue to the next part when you are ready. Please call the experimenter if you have questions'...
    ,'Break','Continue','Continue');
WaitSecs(0.5);

%==========================================================
%% Part 5 - Recognition Go/NoGo
%==========================================================
crashedRecognitionGoNogo = 0;
keepTrying = 1;
while keepTrying < 10
    try
        recognitionGoNoGo_new_Israel(subjectID, test_comp, mainPath, order, sessionNum);
        keepTrying = 10;
    catch
        sca;
        crashedRecognitionGoNogo = crashedRecognitionGoNogo + 1;
        keepTrying = keepTrying +1;
        disp('CODE HAD CRASHED - RECOGNITION GO NOGO!');
    end
end
fprintf(fid_crash,'RecognitionNewOld crashed:\t %d\n', crashedRecognitionGoNogo);

%==========================================================
%%   'post-task'
%==========================================================
% Break    
questdlg('Thank you. The experiment is over!'...
    ,'End of experiment','End','End');
WaitSecs(0.5);

fclose(fid_crash);
sca;

Copy2Dropbox(subjectID);

WaitSecs(5);
% quit %quits matlab

% Return to the defualt UI font size
set(0,'defaultUicontrolFontSize', old_UI_font_size)
end % end function
