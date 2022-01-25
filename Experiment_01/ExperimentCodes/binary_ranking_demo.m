
function binary_ranking_demo(subjectID,test_comp,~)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ======================== by Tom Salomon 2014-2015 =======================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% Ranking using a set of binary choices - DEMO
% This function was written based on the probe code. In each trial two images
% will appear on screen, and the subject is asked to choose one, using the
% 'u' and 'i' keys.
% There are 8 Demo items and 4 trials. so each stimulus will be presented
% exactly one time.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'random_stimlist_generator'

%=========================================================================
%% dummy info for testing purposes
%=========================================================================

% subjectID='BM2_000';
% test_comp=0;
% mainPath=pwd

%==============================================
%% GLOBAL VARIABLES
%==============================================
% default input
% - - - - - - - - - - - - - - - - -
% if main path and test computer were not inserted, the function assume you
% run in the test room and the main path is the current dir (pwd).
if nargin<=2
    if nargin==1
        test_comp=0;
    end
end

% timestamp
% - - - - - - - - - - - - - - - - -
c = clock;
hr = sprintf('%02d', c(4));
min = sprintf('%02d', c(5));
timestamp=[date,'_',hr,'h',min,'m'];
% rand('state',sum(100*clock));

% set phase times
% - - - - - - - - - - - - - - - - -
maxtime=2;      % 1.5 second limit on each selection
baseline_fixation_dur=0.5; % Need to modify based on if first few volumes are saved or not
afterrunfixation=0.5;

tic

%==============================================
%% Read in data
%==============================================

% load image arrays
% - - - - - - - - - - - - - - -
% Read all image files in the 'stim' folder. make sure the ending is suitable
% to your stimuli, and the folder contains only the experiment's 60 images.
stimuli=dir('./stim/demo/*.jpg');
stimname=struct2cell(rmfield(stimuli,{'date','bytes','isdir','datenum'}));
Images=cell(length(stimname),1);
for i=1:length(stimname)
    Images{i}=imread(sprintf('%s','./stim/demo/',stimname{i}));
end

% Load Hebrew instructions image files
Instructions=dir('./Instructions/*binary_ranking_demo.JPG');
Instructions_name=struct2cell(rmfield(Instructions,{'date','bytes','isdir','datenum'}));
Instructions_image=imread(['./Instructions/' sprintf(Instructions_name{1})]);%sprintf(Introduction_name{1})

Introduction=dir('./Instructions/*Introduction.JPG');
Introduction_name=struct2cell(rmfield(Introduction,{'date','bytes','isdir','datenum'}));
Introduction_image=imread(['./Instructions/' sprintf(Introduction_name{1})]); %sprintf(Introduction_name{1})

% load onsets
% - - - - - - - - - - - - - - -
onsetlist=0:2.5:1200;

%==============================================
%% 'INITIALIZE Screen variables'
%==============================================
Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));
%Screen('Preference', 'SkipSyncTests', 1);
pixelSize=32;
%[w] = Screen('OpenWindow',screennum,[],[0 0 1000 700],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

% Define Colors
% - - - - - - - - - - - - - - -
black=BlackIndex(w); % Should equal 0.
white=WhiteIndex(w); % Should equal 255.
green=[0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);


% text stuffs
% - - - - - - - - - - - - - - -
theFont='Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

% Define image scale - Change according to your stimuli
% - - - - - - - - - - - - - - -
stackH= size(Images{1},1);
stackW= size(Images{1},2);

% Frame and stack properties
% - - - - - - - - - - - - - - -
[wWidth, wHeight]=Screen('WindowSize', w);
xcenter=wWidth/2;
ycenter=wHeight/2;
xDist = 300; % distance from center in the x axis. can be changed to adjust to smaller screens
leftRect=[xcenter-stackW-xDist ycenter-stackH/2 xcenter-xDist ycenter+stackH/2]; % left stack location
rightRect=[xcenter+xDist ycenter-stackH/2 xcenter+stackW+xDist ycenter+stackH/2]; % right stack location
penWidth=10; % frame width

HideCursor;

%==============================================
%% 'ASSIGN response keys'
%==============================================
KbName('UnifyKeyNames');
switch test_comp
    case {0,2,3}
        leftstack='u';
        rightstack= 'i';
        badresp='x';
    case 1
        leftstack='b';
        rightstack= 'y';
        badresp='x';
end

%==============================================
%%   'PRE-TRIAL DATA ORGANIZATION'
%==============================================
% Stimuli lists
% - - - - - - - - - - - - - - -
n=length(stimname); % number of stimuli
%number_of_trials=2; %for debug
number_of_trials=8; % desired number of trials (number of unique comparisons)
% IMPORTANT NOTE: every stimulus will be presented exactly 2*number_of_trials/n times. Therefore, number_of_trials should be a multiple of n/2. e.g if n=60, number_of_trials can be 30,60,90,120....
% ==============                                                                                  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[shuffle_stimlist1,shuffle_stimlist2]=random_stimlist_generator(n,number_of_trials);

%-----------------------------------------------------------------
%% 'Write output file header'
%-----------------------------------------------------------------
fid1=fopen([ './Output/' subjectID '_demo_bianry_ranking_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\t runtrial\t onsettime\t ImageLeft\t ImageRight\t StimNumLeft\t StimNumRight\t Response\t Outcome\t RT\n'); %write the header line

%==============================================
%% 'Display Main Instructions'
%==============================================

% Show Introduction
Screen('PutImage',w,Introduction_image);
Screen(w,'Flip');

KbQueueCreate;
KbQueueStart;

noresp=1;
while noresp,
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp,
        noresp=0;
    end;
end;

KbQueueFlush;
KbQueueStart;

% Show Instructions
Screen('PutImage',w,Instructions_image);
Screen(w,'Flip');
KbQueueFlush;
KbQueueStart;
WaitSecs(0.5);
noresp=1;
while noresp,
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp,
        noresp=0;
    end;
end;

% baseline fixation cross
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
runtrial=1;
runStart=GetSecs;

for trial=1:number_of_trials
    
    chose_rand=rand; % randomly chose left/right location od stimuli
    if chose_rand<=0.5
        leftname=stimname(shuffle_stimlist1(trial));
        rightname=stimname(shuffle_stimlist2(trial));
    else
        leftname=stimname(shuffle_stimlist2(trial));
        rightname=stimname(shuffle_stimlist1(trial));
    end
    
    out=999;
    
    %-----------------------------------------------------------------
    % display images
    %-----------------------------------------------------------------
    if chose_rand<=0.5
        Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, leftRect);
        Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, rightRect);
    else
        Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, leftRect);
        Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, rightRect);
    end
    CenterText(w,'+', white,0,0);
    StimOnset=Screen(w,'Flip', runStart+onsetlist(runtrial)+baseline_fixation);
    
    %-----------------------------------------------------------------
    % get response
    %-----------------------------------------------------------------
    KbQueueFlush;
    KbQueueStart;
    noresp=1;
    goodresp=0;
    while noresp
        
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck;
        
        if keyIsDown && noresp
            keyPressed=KbName(firstPress);
            
            if ischar(keyPressed)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed=char(keyPressed);
                keyPressed=keyPressed(1);
            end
            
            switch keyPressed
                case leftstack
                    respTime=firstPress(KbName(leftstack))-StimOnset;
                    noresp=0;
                    goodresp=1;
                case rightstack
                    respTime=firstPress(KbName(rightstack))-StimOnset;
                    noresp=0;
                    goodresp=1;
            end
        end
        
        % check for reaching time limit
        if noresp && GetSecs-runStart >= onsetlist(runtrial)+baseline_fixation+maxtime
            noresp=0;
            keyPressed=badresp;
            respTime=maxtime;
        end
    end
    
    %-----------------------------------------------------------------
    % determine what bid to highlight
    %-----------------------------------------------------------------
    
    switch keyPressed
        case leftstack
            if shuffle_stimlist2(trial)==0
                out=1;
            else
                out=0;
            end
        case rightstack
            if shuffle_stimlist2(trial)==1
                out=1;
            else
                out=0;
            end
    end
    
    if goodresp==1 % if responded: add green rectangle around selected image
        if chose_rand<=0.5
            Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, leftRect);
            Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, rightRect);
        else
            Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, leftRect);
            Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, rightRect);
        end
        
        if keyPressed==leftstack
            Screen('FrameRect', w, green, leftRect, penWidth);
        elseif keyPressed==rightstack
            Screen('FrameRect', w, green, rightRect, penWidth);
        end
        
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip',runStart+onsetlist(trial)+respTime+baseline_fixation);
        
    else % if did not respond: show text 'You must respond faster!'
        CenterText(w,sprintf('You must respond faster!') ,white,0,0);
%         Screen('DrawText', w, 'You must respond faster!', xcenter-450, ycenter, white);
        Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+baseline_fixation);
    end
    
    %-----------------------------------------------------------------
    % show fixation ITI
    %-----------------------------------------------------------------
    
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+.5+baseline_fixation);
    
    if goodresp ~= 1
        respTime=999;
    end
    
    %-----------------------------------------------------------------
    % write to output file
    %-----------------------------------------------------------------
    if chose_rand<=0.5
        fprintf(fid1,'%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n', subjectID, runtrial, StimOnset-runStart, char(leftname), char(rightname), shuffle_stimlist1(trial), shuffle_stimlist2(trial), keyPressed, out, respTime*1000);
    else
        fprintf(fid1,'%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n', subjectID, runtrial, StimOnset-runStart, char(leftname), char(rightname), shuffle_stimlist2(trial), shuffle_stimlist1(trial), keyPressed,  out, respTime*1000);
    end
    
    % end of current trial
    runtrial=runtrial+1;
    KbQueueFlush;
end % loop through trials

fclose(fid1);
Postexperiment=GetSecs;
while GetSecs < Postexperiment+afterrunfixation;
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
end

%-----------------------------------------------------------------
%	display end of part message
%-----------------------------------------------------------------
WaitSecs(2);
Screen('FillRect', w, black);
Screen('TextSize',w, 40);
CenterText(w,sprintf('Thank You') ,white,0,-60);
CenterText(w,sprintf('Any questions?') ,white,0,0);
Screen('Flip',w);
WaitSecs(3);
Screen('CloseAll');

%---------------------------------------------------------------
% create a data structure with info about the run
%---------------------------------------------------------------
outfile=strcat('./Output/', sprintf('%s_demo_binary_ranking_%s',subjectID,timestamp),'.mat');

% create a data structure with info about the run
run_info.subject=subjectID;
run_info.date=date;
run_info.outfile=outfile;

run_info.script_name=mfilename;
clear ('Images','Instructions_image','Introduction_image');
save(outfile);

KbQueueFlush;

end
