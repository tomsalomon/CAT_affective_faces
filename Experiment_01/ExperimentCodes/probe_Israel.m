function probe_Israel(subjectID, order, mainPath, test_comp, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun)

% function probe_Israel(subjectID, order, mainPath, test_comp, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik November 2014 ====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% This function runs the probe session of the boost (cue-approach) task. In
% this session, the subject is presented with two items in each trial, and
% should choose which one he prefers within 1.5 seconds. At the end of the
% experiment the function 'probeResolve_Israel' chooses a random trial from
% this probe session and the subject is given the item he chose on that
% chosen trial.
% This function runs one run each time. Each block is devided to
% 'numRunsPerBlock' runs. The stimuli for each run are shuffled, chosen and
% organized in the 'organizeProbe_Israel' function.


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''stopGoList_allstim_order*.txt'' --> created by sortBDM_Israel
%   ''stimuliForProbe_order%d_block_%d_run%d.txt'' --> Created by
%   organizeProbe_Israel

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''probe_block_' block '_' timestamp '.txt''

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'tal_900';
% order = 2;
% test_comp = 0;
% sessionNum = 1;
% mainPath = pwd;
% block = 1;
% numRun = 1;

 %trialsPerRun = 8; % for debugging

tic

%==============================================
%% 'GLOBAL VARIABLES'
%==============================================

% 'block' indicates how many times the subject has completed all
% trials of the probe experiment. Therefore, the first time the subject completes
% a probe block on his/her 2nd session, it is actually probe block 3 for
% that person:

if test_comp == 1 % If it's an fMRI experiment, let the experimenter insert the block number, for him to control when to start
    okblock = [1 2];
    which_block = input('Enter Block 1 or 2 ');
    while isempty(which_block) || sum(okblock == which_block) ~=1
        disp('ERROR: input must be 1 or 2. Please try again.');
        which_block=input('Enter Block 1 or 2 ');
    end
    block = (sessionNum-1)*2 + which_block;
end % end if test_comp == 1

outputPath = [mainPath '/Output'];

%   'timestamp'
% - - - - - - - - - - - - - - - - -
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

%   'set phase times'
% - - - - - - - - - - - - - - - - -
maxtime = 1.5;      % 1.5 second limit on each selection
baseline_fixation_dur = 2; % Need to modify based on if first few volumes are saved or not
afterrunfixation = 6;

% essential for randomization
rng('shuffle');

tic

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_stopGoList_allstim_order*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
stimName = data{1};
bidIndex = data{3};
bidValue = data{4};
fclose(fid);

%   'read in organized list of stimuli for this run'
% - - - - - - - - - - - - - - - - - - - - - - - - - -
fid = fopen([outputPath '/' sprintf('%s_stimuliForProbe_order%d_block_%d_run%d.txt',subjectID,order,block,numRun)]);
stimuli = textscan(fid, '%d %d %d %d %s %s') ;% these contain everything from the organizedProbe
stimnum1 = stimuli{1};
stimnum2 = stimuli{2};
leftGo = stimuli{3};
pairType = stimuli{4};
leftname = stimuli{5};
rightname = stimuli{6};
fclose(fid);

%   'load image arrays'
% - - - - - - - - - - - - - - -
Images = cell(1,length(stimName));
for i = 1:length(stimName)
    Images{i}=imread([mainPath,'/stim/',stimName{i}]);
end

%  Load Hebrew instructions image files
if sessionNum==1
    Instructions=dir([mainPath '/Instructions/*probe.JPG' ]);
    Instructions_name=struct2cell(rmfield(Instructions,{'date','bytes','isdir','datenum'}));
    Instructions_image=imread([mainPath '/Instructions/' sprintf(Instructions_name{1})]);
else
    Instructions=dir([mainPath '/Instructions/recall/*probe.JPG' ]);
    Instructions_name=struct2cell(rmfield(Instructions,{'date','bytes','isdir','datenum'}));
    Instructions_image=imread([mainPath '/Instructions/recall/' sprintf(Instructions_name{1})]);
end
%   'load onsets'
% - - - - - - - - - - - - - - -
r = Shuffle(1:4);
onsetlist = load([mainPath '/Onset_files/probe_onset_length_68_' num2str(r(1)) '.mat']);
onsetlist = onsetlist.onsetlist;

%==============================================
%% 'INITIALIZE Screen variables'
%==============================================
Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize=32;
%[w] = Screen('OpenWindow',screennum,[],[0 0 1000 700],pixelSize/4);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

HideCursor;

% Define Colors
% - - - - - - - - - - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
green = [0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% setting the text
% - - - - - - - - - - - - - - -
theFont = 'Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

% stack locations
% - - - - - - - - - - - - - - -
[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;

% Define image scale - Change according to your stimuli
% - - - - - - - - - - - - - - -
stackH= size(Images{1},1);
stackW= size(Images{1},2);

leftRect = [xcenter-stackW-300 ycenter-stackH/2 xcenter-300 ycenter+stackH/2];
rightRect = [xcenter+300 ycenter-stackH/2 xcenter+stackW+300 ycenter+stackH/2];

penWidth = 10;

HideCursor;

%==============================================
%% 'ASSIGN response keys'
%==============================================
KbName('UnifyKeyNames');
%MRI=0;
switch test_comp
    case {0,2,3}
        leftstack = 'u';
        rightstack = 'i';
        badresp = 'x';
    case 1
        leftstack = 'b';
        rightstack = 'y';
        badresp = 'x';
end

%-----------------------------------------------------------------
%% 'Write output file header'
%-----------------------------------------------------------------

fid1 = fopen([outputPath '/' subjectID sprintf('_probe_block_%02d_run_%02d_', block, numRun) timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\tscanner\torder\tblock\trun\ttrial\tonsettime\tImageLeft\tImageRight\tbidIndexLeft\tbidIndexRight\tIsleftGo\tResponse\tPairType\tOutcome\tRT\tbidLeft\tbidRight\n'); %write the header line

%==============================================
%% 'Display Main Instructions'
%==============================================
KbQueueCreate;
Screen('TextSize',w, 40);

%%% While they are waiting for the trigger
if test_comp == 1
    Screen('TextSize',w, 60);
    CenterText(w,'PART 4:', white,0,-325);
    
    Screen('TextSize',w, 40);
    CenterText(w,'In this part two pictures of faces will be presented on the screen.', white,0,-270);
    CenterText(w,'For each trial, we want you to choose one of the faces using the keyboard.', white,0,-215);
    CenterText(w,'You will have 1.5 seconds to make your choice on each trial, so please', white,0,-160);
    CenterText(w,'try to make your choice quickly.', white,0,-105);
    
    CenterText(w,'Press the keys on the keypad for the left and right items respectively.', white,0,60);
    CenterText(w,'Waiting for trigger...GET READY....', white, 0, 180);
    
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
    
    DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected
else
    if numRun == 1 && (block ==1||block ==3) % If this is the first run of the first or third block, show instructions
        Screen('PutImage',w,Instructions_image);
        Screen('Flip',w);
        
        noresp = 1;
        while noresp,
            [keyIsDown] = KbCheck(-1);%deviceNumber=keyboard
            if keyIsDown && noresp,
                noresp = 0;
            end;
        end;
    else % if this is not the first run of the first block
        Screen('TextSize',w, 40);
        CenterText(w,'Another run begins now.', white,0,-200);
        CenterText(w,'Press any key to continue.', green,0,0);
        Screen('Flip',w);
        WaitSecs(2);
        
        noresp = 1;
        while noresp,
            [keyIsDown] = KbCheck(-1);%deviceNumber=keyboard
            if keyIsDown && noresp,
                noresp = 0;
            end;
        end;
        
    end % end if numRun == 1 && block ==1
    
end % end if test_comp==1

%   baseline fixation cross
% - - - - - - - - - - - - -
prebaseline = GetSecs;
% baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
while GetSecs < prebaseline+baseline_fixation_dur
    %    Screen(w,'Flip', anchor);
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
    
end
postbaseline = GetSecs;
baseline_fixation = postbaseline - prebaseline;


%==============================================
%% 'Run Trials'
%==============================================
runtrial = 1;
runStart = GetSecs;

%for trial = 1:8 % for debugging
for trial = 1:trialsPerRun
    
    
    % initial box outline colors
    % - - - - - - -
    out = 999;
    
    %-----------------------------------------------------------------
    % display images
    %-----------------------------------------------------------------
    if leftGo(trial) == 1
        Screen('PutImage',w,Images{stimnum1(trial)}, leftRect);
        Screen('PutImage',w,Images{stimnum2(trial)}, rightRect);
    else
        Screen('PutImage',w,Images{stimnum2(trial)}, leftRect);
        Screen('PutImage',w,Images{stimnum1(trial)}, rightRect);
    end
    CenterText(w,'+', white,0,0);
    StimOnset = Screen(w,'Flip', runStart+onsetlist(runtrial)+baseline_fixation);
    
    KbQueueFlush;
    KbQueueStart;
    
    %-----------------------------------------------------------------
    % get response
    %-----------------------------------------------------------------
    
    noresp = 1;
    goodresp = 0;
    while noresp
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck;
        
        if keyIsDown && noresp
            keyPressed = KbName(firstPress);
            if ischar(keyPressed) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed = char(keyPressed);
                keyPressed = keyPressed(1);
            end
            switch keyPressed
                case leftstack
                    respTime = firstPress(KbName(leftstack))-StimOnset;
                    noresp = 0;
                    goodresp = 1;
                case rightstack
                    respTime = firstPress(KbName(rightstack))-StimOnset;
                    noresp = 0;
                    goodresp = 1;
            end
        end % end if keyIsDown && noresp
        
        % check for reaching time limit
        if noresp && GetSecs-runStart >= onsetlist(runtrial)+baseline_fixation+maxtime
            noresp = 0;
            keyPressed = badresp;
            respTime = maxtime;
        end
    end % end while noresp
    
    %-----------------------------------------------------------------
    % determine what bid to highlight
    %-----------------------------------------------------------------
    
    switch keyPressed
        case leftstack
            if leftGo(trial) == 0
                out = 0;
            else
                out = 1;
            end
        case rightstack
            if leftGo(trial) == 1
                out = 0;
            else
                out = 1;
            end
    end
    
    if goodresp==1
        if leftGo(trial)==1
            Screen('PutImage',w,Images{stimnum1(trial)}, leftRect);
            Screen('PutImage',w,Images{stimnum2(trial)}, rightRect);
        else
            Screen('PutImage',w,Images{stimnum2(trial)}, leftRect);
            Screen('PutImage',w,Images{stimnum1(trial)}, rightRect);
        end
        
        if keyPressed=='u'
            Screen('FrameRect', w, green, leftRect, penWidth);
        elseif keyPressed=='i'
            Screen('FrameRect', w, green, rightRect, penWidth);
        end
        
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip',runStart+onsetlist(trial)+respTime+baseline_fixation);
        
    else % if did not respond: show text 'You must respond faster!'
        CenterText(w,sprintf('You must respond faster!') ,white,0,0);
        %         Screen('DrawText', w, 'You must respond faster!', xcenter-450, ycenter, white);
        Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+baseline_fixation);
    end % end if goodresp==1
    
    %-----------------------------------------------------------------
    % show fixation ITI
    %-----------------------------------------------------------------
    
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+.5+baseline_fixation);
    
    if goodresp ~= 1
        respTime = 999;
    end
    
    %-----------------------------------------------------------------
    % 'Save data'
    %-----------------------------------------------------------------
    if leftGo(trial)==1
        fprintf(fid1,'%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n', subjectID, test_comp, order, sprintf('%02d', block), numRun, runtrial, StimOnset-runStart, char(leftname(trial)), char(rightname(trial)), stimnum1(trial), stimnum2(trial), leftGo(trial), keyPressed, pairType(trial), out, respTime*1000, bidValue(stimnum1(trial)), bidValue(stimnum2(trial)));
    else
        fprintf(fid1,'%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n', subjectID, test_comp, order, sprintf('%02d', block), numRun, runtrial, StimOnset-runStart, char(leftname(trial)), char(rightname(trial)), stimnum2(trial), stimnum1(trial), leftGo(trial), keyPressed, pairType(trial), out, respTime*1000, bidValue(stimnum2(trial)), bidValue(stimnum1(trial)));
    end
    
    runtrial = runtrial+1;
    %     KbQueueFlush;
    
end % loop through trials
fclose(fid1);

Postexperiment = GetSecs;
while GetSecs < Postexperiment + afterrunfixation;
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
    
end

%-----------------------------------------------------------------
%	display outgoing message
%-----------------------------------------------------------------
WaitSecs(2);
Screen('FillRect', w, black);
Screen('TextSize',w, 40);

if test_comp == 1 % An fMRI experiment. Edit so that it would be good for the fmri experiment
    if numRun ~= numRunsPerBlock
        % This is not the last run of the block
        CenterText(w,sprintf('In a moment we will complete another run of the same task') ,white,0,-170);
        Screen('Flip',w);
        WaitSecs(3);
        Screen('CloseAll');
        ShowCursor;
    elseif mod(block,2) == 1 && numRun == numRunsPerBlock
        CenterText(w,sprintf('In a moment we will complete another block of the same task') ,white,0,-170);
        Screen('Flip',w);
        WaitSecs(3);
        Screen('CloseAll');
        ShowCursor;
    elseif mod(block,2) == 0 && sessionNum == 1 && numRun == numRunsPerBlock
        %if block is an even number & it's the first session & the last run
        CenterText(w,sprintf('Part 2 is done') ,white,0,-270);
        CenterText(w,sprintf('Questions?') ,white,0,-170);
        Screen('Flip',w);
        WaitSecs(3);
        Screen('CloseAll');
        ShowCursor;
    elseif mod(block,2) == 1 && sessionNum > 1 && numRun == numRunsPerBlock
        %if block is an even number & it's a follow-up session and the last
        %run
        CenterText(w,sprintf('Part 4 is done') ,white,0,-270);
        CenterText(w,sprintf('Questions?') ,white,0,-170);
        Screen('Flip',w);
        WaitSecs(3);
        Screen('CloseAll');
        ShowCursor;
    end % end if numRun ~= numRunsPerBlock
    
else % Not an fMRI experiment
    if numRun ~= numRunsPerBlock
        % This is not the last run of the block
        CenterText(w,sprintf('Once you are ready we will complete another run of the same task') ,white,0,-270);
        CenterText(w,sprintf('Press any key to continue.') ,white,0,-170);
        Screen('Flip',w);
        
        noresp=1;
        while noresp,
            [keyIsDown,~,~] = KbCheck;
            if keyIsDown && noresp,
                noresp=0;
            end;
        end;
        Screen('CloseAll');
        ShowCursor;
        WaitSecs(3);
        
    elseif mod(block,2) == 1 && numRun == numRunsPerBlock % The last run of the first block of the session
        CenterText(w,sprintf('Once you are ready we will complete another block of the same task') ,white,0,-270);
        CenterText(w,sprintf('Press any key to continue.') ,white,0,-170);
        Screen('Flip',w);
        
        noresp=1;
        while noresp,
            [keyIsDown,~,~] = KbCheck;
            if keyIsDown && noresp,
                noresp=0;
            end;
        end;
        Screen('CloseAll');
        ShowCursor;
        WaitSecs(3);
        
    elseif mod(block,2) == 0 && sessionNum == 1 && numRun == numRunsPerBlock
        %if block is an even number & it's the first session & the last run
        %of the block
        CenterText(w,sprintf('Please read the Part 5 instructions and continue on your own.') ,white,0,-270);
        CenterText(w,sprintf('Press any key to continue.') ,white,0,-170);
        Screen('Flip',w);
        
        noresp=1;
        while noresp,
            [keyIsDown,~,~] = KbCheck;
            if keyIsDown && noresp,
                noresp=0;
            end;
        end;
        Screen('CloseAll');
        ShowCursor;
        
    elseif mod(block,2) == 0 && sessionNum > 1 && numRun == numRunsPerBlock
        %if block is an even number & it's a follow-up session & the last
        %run of the block
        CenterText(w,sprintf('Please read the Part 2 instructions and continue on your own.') ,white,0,-270);
        CenterText(w,sprintf('Press any key to continue.') ,white,0,-170);
        Screen('Flip',w);
        
        noresp=1;
        while noresp,
            [keyIsDown,~,~] = KbCheck;
            if keyIsDown && noresp,
                noresp=0;
            end;
        end;
        Screen('CloseAll');
        ShowCursor;
    end % end if mod(block,2) == 1
end % end if test)comp == 1

WaitSecs(3);

%---------------------------------------------------------------
% create a data structure with info about the run
%---------------------------------------------------------------
outfile = strcat(outputPath,'/', sprintf('%s_probe_block_%d_run_%d_%s.mat',subjectID,block,numRun,timestamp));

% create a data structure with info about the run
run_info.subject=subjectID;
run_info.date=date;
run_info.outfile=outfile;

run_info.script_name=mfilename;
clear ('Images','Instructions_image');
save(outfile);

KbQueueFlush;

end % end function

