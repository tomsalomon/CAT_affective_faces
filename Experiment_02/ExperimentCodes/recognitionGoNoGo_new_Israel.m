function recognitionGoNoGo_new_Israel(subjectID,test_comp,mainPath,order, sessionNum)

% function recognitionNewOld_Israel(subjectID,test_comp,mainPath,order, sessionNum)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik December 2014 ====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the first part of the memory (recognition) session of the
% boost (cue-approach) task.
% Subjects are presented with the stimuli from the previous sessions (but
% only those that were part of the probe comparisons- as defined hard-coded in this function)
% as well as some items that were not included in the training sesuiiiiision, in a random order, and
% should answer whether they recognize each stimuli from the previous
% sessions or not.
% In the next memory (recognition) session, run by the
% 'recognitionGoNoGo_Israel' function, they are presented with the items
% they said were present during the training session, and are asked whether
% each item was paired with a beep (GO item) or not (NOGO item).

% Old stimuli should be located in the folder [mainPath'/stim/']
% New stimuli should be located in the folder [mainPath'/stim/recognitionNew/']


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order%d.txt', order --> Created by the 'sortBdm_Israel' function


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''recognitionNewOld_results' num2str(sessionNum) '_' timestamp '.txt''
%   --> Which includes the variables: trialNum (in the newOld recognition task), index
%   of the item (ABC order, old and than new), name of item, whether the item
%   is old (1) or not (0), the subject's answer (0- not included; 1-
%   included)

%   ''recognitionNewOld' num2str(sessionNum) '_' timestamp '.txt''


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'BM_9001';
% order = 1;
% test_comp = 3;
% sessionNum = 1;
% mainPath = 'D:/Rotem/Matlab/Boost_Israel_New_Rotem';

tic

%==========================================================
%% 'INITIALIZE Screen variables to be used in each task'
%==========================================================

Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize = 32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 1000 700],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

% % Here Be Colors
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
green = [0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% % set up screen positions for stimuli
% [wWidth, wHeight] = Screen('WindowSize', w);
% xcenter = wWidth/2;
% ycenter = wHeight/2;

% text settings
theFont = 'Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

%---------------------------------------------------------------
%% 'GLOBAL VARIABLES'
%---------------------------------------------------------------
c = clock;
hr = sprintf('%02d', c(4));
min = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',min,'m'];

outputPath = [mainPath '/Output'];

% essential for randomization
rng('shuffle');

%---------------------------------------------------------------
%% 'Assign response keys'
%---------------------------------------------------------------
KbName('UnifyKeyNames');

if test_comp == 1
    leftresp = 'b';
    rightresp = 'g';
    %     badresp = 'x';
else
    leftresp = 'u';
    rightresp = 'i';
    %     badresp = 'x';
end

%---------------------------------------------------------------
%% 'DEFINE stimuli of interest'
%---------------------------------------------------------------

switch order
    case 1
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HV_beep =   [7 10 12 13 15 18]; % HV_beep
        HV_nobeep = [8 9 11 14 16 17 19:22]; % HV_nobeep
        
        LV_beep =   [47 50 52 53 55 58]; % LV_beep
        LV_nobeep = [48 49 51 54 56 57 59:62]; % LV_nobeep
        
        
        %  Happy sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Happy_sanityHV_nobeep = [3 4]; % HV_nobeep
        Happy_sanityLV_nobeep = [37 38]; % LV_nobeep
        
        %  Neutral sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Neutral_sanityHV_nobeep = [43 44]; % HV_nobeep
        Neutral_sanityLV_nobeep = [77 78]; % LV_nobeep
        
    case 2
        
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HV_beep =   [8 9 11 14 16 17]; % HV_beep
        HV_nobeep = [7 10 12 13 15 18 19:22]; % HV_nobeep
        
        
        LV_beep =   [48 49 51 54 56 57]; % LV_beep
        LV_nobeep = [47 50 52 53 55 58 59:62]; % LV_nobeep
        
        
        %  Happy sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Happy_sanityHV_nobeep = [3 4]; % HV_nobeep
        Happy_sanityLV_nobeep = [37 38]; % LV_nobeep
        
        %  Neutral sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Neutral_sanityHV_nobeep = [43 44]; % HV_nobeep
        Neutral_sanityLV_nobeep = [77 78]; % LV_nobeep
        
end % end switch order


%---------------------------------------------------------------
%% 'READ data about the stimuli - Go \ NoGo
%---------------------------------------------------------------

%   'Reading in the sorted BDM list - defines which items will be GO\NOGO'
% - - - - - - - - - - - - - - -
file = dir([outputPath '/' subjectID sprintf('_stopGoList_allstim_order%d.txt', order)]);
fid = fopen([outputPath '/' sprintf(file(length(file)).name)]);
vars = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
fclose(fid);

allOldStimName = vars{1};
allGoNoGo = vars{2};
allIsBeep = zeros(length(allGoNoGo),1);
allIsBeep(allGoNoGo == 11 | allGoNoGo == 22) = 1;
allBidInd = vars{3};
% bidValue = vars{4};
oldStimName = allOldStimName([HV_beep HV_nobeep LV_beep LV_nobeep Happy_sanityHV_nobeep Happy_sanityLV_nobeep Neutral_sanityHV_nobeep Neutral_sanityLV_nobeep]);
% goNoGo = allGoNoGo([HV_beep HV_nobeep LV_beep LV_nobeep]);
isBeep = allIsBeep([HV_beep HV_nobeep LV_beep LV_nobeep   Happy_sanityHV_nobeep Happy_sanityLV_nobeep Neutral_sanityHV_nobeep Neutral_sanityLV_nobeep]);
bidInd = allBidInd([HV_beep HV_nobeep LV_beep LV_nobeep   Happy_sanityHV_nobeep Happy_sanityLV_nobeep Neutral_sanityHV_nobeep Neutral_sanityLV_nobeep]);

[oldStimName, indSortedOldStimName] = sort(oldStimName); % sort old stimuli ABC

%---------------------------------------------------------------
%% 'LOAD image arrays'
%---------------------------------------------------------------

newStimName = dir([mainPath '/stim/recognitionNew/*.jpg']); % Read new stimuli

% Read old images to a cell array
imgArraysOld = cell(1,length(oldStimName));
for i = 1:length(oldStimName)
    imgArraysOld{i} = imread([mainPath '/stim/' oldStimName{i}]);
end
Old = ones(length(oldStimName),1);

% Read new images to a cell array
imgArraysNew = cell(1,length(newStimName));
for i = 1:length(newStimName)
    imgArraysNew{i} = imread([mainPath '/stim/recognitionNew/' newStimName(i).name]);
end
New = zeros(length(newStimName),1);

isOld = [Old; New]; % Create an array indicating whether an item is old (1) or not (0)

% Load Hebrew instructions image files
if sessionNum==1
    Instructions=dir([mainPath '/Instructions/*recognitionGoNoGo.JPG' ]);
    % Instructions_name=struct2cell(rmfield(Instructions,{'date','bytes','isdir','datenum'}));
    Instructions_name={Instructions.name};
    Instructions_image=imread([mainPath '/Instructions/' sprintf(Instructions_name{1})]);
else
    Instructions=dir([mainPath '/Instructions/recall/*recognitionGoNoGo.JPG' ]);
    % Instructions_name=struct2cell(rmfield(Instructions,{'date','bytes','isdir','datenum'}));
    Instructions_name={Instructions.name};
    Instructions_image=imread([mainPath '/Instructions/recall/' sprintf(Instructions_name{1})]);
end

%---------------------------------------------------------------
%% 'ORGANIZE data about the stimuli - Go \ NoGo
%---------------------------------------------------------------

% Add zeros for the goNoGo, isBeep, bidInd and bidValue of the new items
% - - - - - - - - - - - - - - -

sortedIsBeep = isBeep(indSortedOldStimName);
sortedIsBeep(length(oldStimName)+1:length(oldStimName)+length(newStimName)) = 0;

sortedBidInd = bidInd(indSortedOldStimName);
sortedBidInd(length(oldStimName)+1:length(oldStimName)+length(newStimName)) = 0;

%---------------------------------------------------------------
%% 'SHUFFLE data about the stimuli - Go \ NoGo
%---------------------------------------------------------------

% Merge the old and new lists, and shuffle them

% Add the names of the new stimuli to stimName
stimName = cell(1, length(oldStimName) + length(newStimName));
stimName(1:length(oldStimName)) = oldStimName;

for newStimInd = 1:length(newStimName)
    stimName{length(oldStimName)+newStimInd} = newStimName(newStimInd).name;
end

[shuffledlist, shuffledlistInd] = Shuffle(stimName);
imgArrays = [imgArraysOld imgArraysNew];
imgArrays = imgArrays(shuffledlistInd);
shuffledIsOld = isOld(shuffledlistInd);
shuffledSortedIsBeep = sortedIsBeep(shuffledlistInd);
shuffledSortedBidInd = sortedBidInd(shuffledlistInd);


% r = Shuffle(1:4);
% onsetlist = load(['Onset_files\sweet_salty_onset_' num2str(r(1)) '.mat']);
% onsetlist = onsetlist.onsetlist;

% ListenChar(2); % suppresses terminal ouput

%---------------------------------------------------------------
%% 'Write output file header'
%---------------------------------------------------------------

fid1 = fopen([outputPath '/' subjectID '_recognitionGoNoGo' num2str(sessionNum) '_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\t order\t itemIndABC\t runtrial\t isGO?\t subjectAnswer\t onsettime\t Name\t resp_choice\t RT\t bidInd\n'); %write the header line

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------


Screen('PutImage',w,Instructions_image);
Screen(w,'Flip');

HideCursor; % Make a comment in debugging mode

if test_comp == 1
    CenterText(w,'GET READY!', white, 0, 0);
    Screen('Flip',w);
    escapeKey = KbName('space');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
else
    noresp = 1;
    while noresp
        [keyIsDown] = KbCheck; % deviceNumber=keyboard
        if keyIsDown && noresp
            noresp = 0;
        end
    end
end % end if test_comp == 1

WaitSecs(0.001);
% anchor = GetSecs;

Screen('TextSize',w, 60);
CenterText(w,'+', white,0,0);
Screen(w,'Flip');

WaitSecs(1);

KbQueueCreate;

%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------
Included = zeros(1, length(stimName)); % An array for the results. 1 represent items that the subject said were present

runStart = GetSecs;

%for trial = 1:6 % for debugging
for trial = 1:length(stimName)
    
    colorright = white;
    colorleft = white;
    
    %-----------------------------------------------------------------
    % display image
    % self-paced; next image will only show after response
    
    Screen('PutImage',w, imgArrays{trial});
    Screen('TextSize',w, 40);
       
    if order==1
        CenterText(w,'u - Beep',colorleft,-250,300);
        CenterText(w,'i - No beep',colorright,250,300);
    else
        CenterText(w,'u - No beep',colorleft,-250,300);
        CenterText(w,'i - Beep',colorright,250,300);
    end
    
    Screen(w,'Flip');
    StimOnset = GetSecs;
    
    KbQueueStart;
    %-----------------------------------------------------------------
    % get response
    
    
    noresp = 1;
    while noresp
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck(-1);
        
        if keyIsDown && noresp
            findfirstPress = find(firstPress);
            respTime = firstPress(findfirstPress(1))-StimOnset;
            tmp = KbName(findfirstPress);
            if ischar(tmp)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                tmp = char(tmp);
            end
            response = tmp(1);
            if response=='u'||response=='i'
                noresp = 0;
            end
        end % end if keyIsDown && noresp
        
    end % end while noresp
    
    
    
    %-----------------------------------------------------------------
    
    % determine what bid to highlight and write the subject's response to
    % the Included array
    if order ==1
        switch response
            case leftresp
                colorleft = green;
                Included(trial) = 1;
            case rightresp
                colorright = green;
        end
    else % order == 2
        switch response
            case leftresp
                colorleft = green;
            case rightresp
                colorright = green;
                Included(trial) = 1;
        end
    end % end if order == 1
    
    %-----------------------------------------------------------------
    % redraw text output with the appropriate colorchanges
    Screen('PutImage',w,imgArrays{trial});
    Screen('TextSize',w, 40);

    
    if order==1
        CenterText(w,'u - Beep',colorleft,-250,300);
        CenterText(w,'i - No beep',colorright,250,300);
    else
        CenterText(w,'u - No beep',colorleft,-250,300);
        CenterText(w,'i - Beep', colorright,250,300);
    end
    
    Screen(w,'Flip');
    WaitSecs(0.5);
    
    %-----------------------------------------------------------------
    % show fixation ITI
    Screen('TextSize',w, 60);
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip');
    WaitSecs(0.5);
    
    %-----------------------------------------------------------------
    % Write to output files
    
    fprintf(fid1,'%s\t %d\t %d\t %d\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t \n', subjectID, order, shuffledlistInd(trial), trial, shuffledSortedIsBeep(trial), Included(trial), StimOnset-runStart, shuffledlist{trial}, response, respTime, shuffledSortedBidInd(trial));    
    KbQueueFlush;
    
end % end loop for trial = 1:length(stimName);

% Close open files
fclose(fid1);

% Save variables to mat file
outfile = strcat(outputPath,'/', sprintf('%s_recognitionGoNoGo_%s.mat',subjectID, timestamp));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;

run_info.script_name = mfilename;
clear imgArrays imgArraysNew imgArraysOld Instructions_image;
save(outfile);

% End of session screen
%{
if sessionNum == 1
    Screen('TextSize',w, 40);
    CenterText(w,'Thank you! Please read the instructions and proceed with PART 6.', white,0,0);
    Screen(w,'Flip');
else
    Screen('TextSize',w, 40);
    CenterText(w,'Thank you! Please read the instructions and proceed with PART 3.', white,0,0);
    Screen(w,'Flip');
end
%}
% Closing
WaitSecs(4);
toc
ShowCursor;
Screen('CloseAll');

end % end function