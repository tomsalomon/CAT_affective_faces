function [trialsPerRun] = organizeProbe_Israel(subjectID, order, mainPath, block, numRunsPerBlock)

% function [trialsPerRun] = organizeProbe_Israel(subjectID, order, mainPath, block, numRunsPerBlock)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik December 2014 ====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function organizes the matrices for each block of the probe session of the boost
% (cue-approach) task, divided to number of runs as requested (1 or 2 or 4 would
% work. Pay attention that number of runs should be a divisor of number of
% comparisons.

% This function is for the version where only 40 items are being trained,
% and the sanity checks are only on the NOGO items


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order*.txt'' --> created by sortBDM_Israel


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stimuliForProbe_order%d_block_%d_run%d.txt'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'test999';
% order = 1;
% test_comp = 4;
% mainPath = '/Users/schonberglabimac1/Documents/Boost_Israel_New_Rotem_mac';
% numRunsPerBlock = 1;
% block = 1;


tic

%==============================================
%% 'GLOBAL VARIABLES'
%==============================================

% essential for randomization
rng('shuffle');

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_stopGoList_allstim_order*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
stimName = data{1};
% bidIndex = data{3};
% bidValue = data{4};
fclose(fid);

%==============================================
%%   'DATA ORGANIZATION'
%==============================================

% determine stimuli to use based on order number
%-----------------------------------------------------------------
switch order
    case 1
        %   comparisons of interest - Happy/Neutral GO vs. NoGo
        % - - - - - - - - - - - - - - -
        HV_beep =   [7 10 12 13 15 18]; % HV_beep
        HV_nobeep = [8 9 11 14 16 17]; % HV_nobeep
        
        LV_beep =   [47 50 52 53 55 58]; % LV_beep
        LV_nobeep = [48 49 51 54 56 57]; % LV_nobeep
        
        
        %  Happy sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Happy_sanityHV_nobeep = [3 4]; % HV_nobeep
        Happy_sanityLV_nobeep = [37 38]; % LV_nobeep
        
        %  Neutral sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Neutral_sanityHV_nobeep = [43 44]; % HV_nobeep
        Neutral_sanityLV_nobeep = [77 78]; % LV_nobeep
        
        
        %  Happy Fake-Go comparisons
        % - - - - - - - - - - - - - -
        Happy_Fake_Go = [19 22];
        Happy_Fake_NoGo = [20 21];
        
        %  Neutral Fake-Go comparisons
        % - - - - - - - - - - - - - -
        Neutral_Fake_Go = [59 62];
        Neutral_Fake_NoGo = [60 61];
        
    case 2
        
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HV_beep =   [8 9 11 14 16 17]; % HV_beep
        HV_nobeep = [7 10 12 13 15 18]; % HV_nobeep
        
        
        LV_beep =   [48 49 51 54 56 57]; % LV_beep
        LV_nobeep = [47 50 52 53 55 58]; % LV_nobeep
        
        
        %  Happy sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Happy_sanityHV_nobeep = [3 4]; % HV_nobeep
        Happy_sanityLV_nobeep = [37 38]; % LV_nobeep
        
        %  Neutral sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        Neutral_sanityHV_nobeep = [43 44]; % HV_nobeep
        Neutral_sanityLV_nobeep = [77 78]; % LV_nobeep
        
        
        %  Happy Fake-Go comparisons
        % - - - - - - - - - - - - - -
        Happy_Fake_Go = [20 21];
        Happy_Fake_NoGo = [19 22];
        
        %  Neutral Fake-Go comparisons
        % - - - - - - - - - - - - - -
        Neutral_Fake_Go = [60 61];
        Neutral_Fake_NoGo = [59 62];
        
end % end switch order


%   add multiple iterations of each item presentation
%-----------------------------------------------------


%   TRIAL TYPE 1: HighValue Go vs. HighValue NoGo(Stop) - Happy Go vs. NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHVbeepItems = length(HV_beep); % How many happy Go stimuli
numHVnobeepItems = length(HV_nobeep); % How many happy No Go stimuli

HV_beep_new = repmat(HV_beep,numHVbeepItems,1); % Happy Go as a numHVbeepItems X length(HV_beep) matrix
HV_beep_new = HV_beep_new(:)';
HV_nobeep_new = repmat(HV_nobeep,1,numHVnobeepItems); % Happy No Go as a numHVnobeepItems X length(HV_nobeep) matrix
HV_nobeep_new = HV_nobeep_new(:)';

[shuffle_HV_beep_new,shuff_HV_beep_new_ind] = Shuffle(HV_beep_new); % Shuffles each column seperatley and index each column seperatley.
shuffle_HV_nobeep_new = HV_nobeep_new(shuff_HV_beep_new_ind); % Creates maching matrix of No Go stimuli.



%   TRIAL TYPE 2: LowValue Go vs. LowValue NoGo(Stop) - Neutral Go vs. NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numLVbeepItems = length(LV_beep); % same procedure as TRIAL TYPE 1
numLVnobeepItems = length(LV_nobeep);

LV_beep_new = repmat(LV_beep,numLVbeepItems,1);
LV_beep_new = LV_beep_new(:)';
LV_nobeep_new = repmat(LV_nobeep,1,numLVnobeepItems);
LV_nobeep_new = LV_nobeep_new(:)';

[shuffle_LV_beep_new,shuff_LV_beep_new_ind] = Shuffle(LV_beep_new);
shuffle_LV_nobeep_new = LV_nobeep_new(shuff_LV_beep_new_ind);


%   TRIAL TYPE 3: Happy HighValue NoGo(Stop) vs. LowValue NoGo(Stop) -
%   Happy sanity
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHappySanityHVnobeepItems = length(Happy_sanityHV_nobeep); % same procedure as TRIAL TYPE 1
numHappySanityLVnobeepItems = length(Happy_sanityLV_nobeep);

HappysanityHV_nobeep_new = repmat(Happy_sanityHV_nobeep,numHappySanityHVnobeepItems,1);
HappysanityHV_nobeep_new = HappysanityHV_nobeep_new(:)';
HappysanityLV_nobeep_new = repmat(Happy_sanityLV_nobeep,1,numHappySanityLVnobeepItems);
HappysanityLV_nobeep_new = HappysanityLV_nobeep_new(:)';

[shuffle_HappysanityHV_nobeep_new,shuff_HappysanityHV_nobeep_new_ind] = Shuffle(HappysanityHV_nobeep_new);
shuffle_HappysanityLV_nobeep_new = HappysanityLV_nobeep_new(shuff_HappysanityHV_nobeep_new_ind);


%   TRIAL TYPE 4: Neutral HighValue NoGo(Stop) vs. LowValue NoGo(Stop) -
%   Neutral sanity
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numNeutralSanityHVnobeepItems = length(Neutral_sanityHV_nobeep); % same procedure as TRIAL TYPE 1
numNeutralSanityLVnobeepItems = length(Neutral_sanityLV_nobeep);

NeutralsanityHV_nobeep_new = repmat(Neutral_sanityHV_nobeep,numNeutralSanityHVnobeepItems,1);
NeutralsanityHV_nobeep_new = NeutralsanityHV_nobeep_new(:)';
NeutralsanityLV_nobeep_new = repmat(Neutral_sanityLV_nobeep,1,numNeutralSanityLVnobeepItems);
NeutralsanityLV_nobeep_new = NeutralsanityLV_nobeep_new(:)';

[shuffle_NeutralsanityHV_nobeep_new,shuff_NeutralsanityHV_nobeep_new_ind] = Shuffle(NeutralsanityHV_nobeep_new);
shuffle_NeutralsanityLV_nobeep_new = NeutralsanityLV_nobeep_new(shuff_NeutralsanityHV_nobeep_new_ind);


%   TRIAL TYPE 5: Happy Fake Go vs. Happy no GO
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHappyFakeGoItems = length(Happy_Fake_Go); % same procedure as TRIAL TYPE 1
numHappyFakeNoGoItems = length(Happy_Fake_NoGo);

Happy_FakeGo_new = repmat(Happy_Fake_Go,numHappyFakeGoItems,1);
Happy_FakeGo_new = Happy_FakeGo_new(:)';
Happy_FakeNoGo_new = repmat(Happy_Fake_NoGo,1,numHappyFakeNoGoItems);
Happy_FakeNoGo_new = Happy_FakeNoGo_new(:)';

[shuffle_Happy_FakeGo_new,shuffle_Happy_FakeGo_new_ind] = Shuffle(Happy_FakeGo_new);
shuffle_Happy_FakeNoGo_new = Happy_FakeNoGo_new(shuffle_Happy_FakeGo_new_ind);



%   TRIAL TYPE 6: Neutral Fake Go vs. Neutral no GO
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numNeutralFakeGoItems = length(Neutral_Fake_Go); % same procedure as TRIAL TYPE 1
numNeutralFakeNoGoItems = length(Neutral_Fake_NoGo);

Neutral_FakeGo_new = repmat(Neutral_Fake_Go,numNeutralFakeGoItems,1);
Neutral_FakeGo_new = Neutral_FakeGo_new(:)';
Neutral_FakeNoGo_new = repmat(Neutral_Fake_NoGo,1,numNeutralFakeNoGoItems);
Neutral_FakeNoGo_new = Neutral_FakeNoGo_new(:)';

[shuffle_Neutral_FakeGo_new,shuffle_Neutral_FakeGo_new_ind] = Shuffle(Neutral_FakeGo_new);
shuffle_Neutral_FakeNoGo_new = Neutral_FakeNoGo_new(shuffle_Neutral_FakeGo_new_ind);




%   randomize all possible comparisons for all trial types
%-----------------------------------------------------------------
numComparisonsHV = numHVbeepItems^2; % How many Happy Go-No Go comparisons
numComparisonsLV = numLVbeepItems^2; % How many Neutral Go-No Go comparisons
numComparisons_HappyFake = numHappyFakeGoItems^2; % How many Happy Fake Go comparisons
numComparisons_NeutralFake = numNeutralFakeGoItems^2; % How many Neutral Fake Go comparisons
numFake = numComparisons_HappyFake + numComparisons_NeutralFake; %Total number of Fake Go
numComparisons = numComparisonsHV + numComparisonsLV; % How many Happy / neutral Go-No Go
numSanity = numHappySanityHVnobeepItems^2 + numNeutralSanityHVnobeepItems^2; % How many Sanity comparisons
total_num_trials = numComparisons + numSanity + numFake; % total comparisons
trialsPerRun = total_num_trials/numRunsPerBlock; % Trials per run

stimnum1 = zeros(numRunsPerBlock,trialsPerRun); % for us - 2X44 zeros (pairs)
stimnum2 = zeros(numRunsPerBlock,trialsPerRun);
leftname = cell(numRunsPerBlock,trialsPerRun);
rightname = cell(numRunsPerBlock,trialsPerRun);
pairType = zeros(numRunsPerBlock,trialsPerRun);


numComparisonsPerRun = numComparisons/numRunsPerBlock; % number of non-Sanity non-FakeGo comparisons per run
numSanityPerRun = numSanity/numRunsPerBlock; % numbe of sanity comparisons per run
numFakePerRun = numFake/numRunsPerBlock; % numbe of Fake Go comparisons per run
pairType(1:numRunsPerBlock,1:numComparisonsPerRun/2) = 1; % pairs 1:18
pairType(1:numRunsPerBlock,numComparisonsPerRun/2+1:numComparisonsPerRun) = 2; % pairs 19:36
pairType(1:numRunsPerBlock,numComparisonsPerRun+1:numComparisonsPerRun+numSanityPerRun/2) = 3; % pairs 37:38
pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun/2+1:numComparisonsPerRun+numSanityPerRun) = 4; % pairs 39:40
pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun+1:numComparisonsPerRun+numSanityPerRun+numFakePerRun/2) = 5; % pairs 41:42
pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun+numFakePerRun/2+1:numComparisonsPerRun+numSanityPerRun+numFakePerRun) = 6; % pairs 33:44


leftGo = ones(numRunsPerBlock,total_num_trials/numRunsPerBlock);
leftGo(:,1:numComparisonsPerRun/2 + numSanityPerRun/2 + numFakePerRun/2) = 0;
%leftGo(:,[1:numComparisonsPerRun/4 numComparisonsPerRun/2+1:numComparisonsPerRun*3/4 1+numComparisonsPerRun:numComparisonsPerRun+numSanityPerRun/2 numComparisonsPerRun+numSanityPerRun+1:numComparisonsPerRun+numSanityPerRun+numFakePerRun/2]) = 0;


for numRun = 1:numRunsPerBlock
    pairType(numRun,:) = Shuffle(pairType(numRun,:));
    leftGo(numRun,:) = Shuffle(leftGo(numRun,:));
end % end for numRun = 1:numRunsPerBlock

HV_beep = shuffle_HV_beep_new;
HV_nobeep = shuffle_HV_nobeep_new;

LV_beep = shuffle_LV_beep_new;
LV_nobeep = shuffle_LV_nobeep_new;

Happy_sanityHV_nobeep = shuffle_HappysanityHV_nobeep_new;
Happy_sanityLV_nobeep = shuffle_HappysanityLV_nobeep_new;

Neutral_sanityHV_nobeep = shuffle_NeutralsanityHV_nobeep_new;
Neutral_sanityLV_nobeep = shuffle_NeutralsanityLV_nobeep_new;

Happy_Fake_Go = shuffle_Happy_FakeGo_new;
Happy_Fake_NoGo = shuffle_Happy_FakeNoGo_new;

Neutral_Fake_Go = shuffle_Neutral_FakeGo_new;
Neutral_Fake_NoGo = shuffle_Neutral_FakeNoGo_new;

% % Divide the matrices of each comparison to the number of trials
% HV_beep_allRuns = reshape(HV_beep,2,length(HV_beep)/2);
% HV_nobeep_allRuns = reshape(HV_nobeep,2,length(HV_nobeep)/2);
% LV_beep_allRuns = reshape(LV_beep,2,length(LV_beep)/2);
% LV_nobeep_allRuns = reshape(LV_nobeep,2,length(LV_nobeep)/2);
% sanityHV_nobeep_allRuns = reshape(sanityHV_nobeep,2,length(sanityHV_nobeep)/2);
% sanityLV_nobeep_allRuns = reshape(sanityLV_nobeep,2,length(sanityLV_nobeep)/2);

HH = 1;
LL = 1;
HL_S = 1;
HL_G = 1;
HH_F = 1;
NN_F = 1;
for numRun = 1:numRunsPerBlock
    
    % Create stimuliForProbe.txt for this run
    fid1 = fopen(['./Output/' sprintf('%s_stimuliForProbe_order%d_block_%d_run%d.txt',subjectID,order,block,numRun)], 'w');
    
    for trial = 1:trialsPerRun % trial num within block
        switch pairType(numRun,trial)
            case 1
                
                % HighValue Go vs. HighValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = HV_beep(HH);
                stimnum2(numRun,trial) = HV_nobeep(HH);
                HH = HH+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 2
                
                % LowValue Go vs. LowValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = LV_beep(LL);
                stimnum2(numRun,trial) = LV_nobeep(LL);
                LL = LL+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 3
                
                % HighValue Go vs. LowValue Go
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = Happy_sanityHV_nobeep(HL_S);
                stimnum2(numRun,trial) = Happy_sanityLV_nobeep(HL_S);
                HL_S = HL_S+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 4
                
                % HighValue NoGo(Stop) vs. LowValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = Neutral_sanityHV_nobeep(HL_G);
                stimnum2(numRun,trial) = Neutral_sanityLV_nobeep(HL_G);
                HL_G = HL_G+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 5
                
                % Happy Fake Go(Stop) vs. Happy Fake NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = Happy_Fake_Go(HH_F);
                stimnum2(numRun,trial) = Happy_Fake_NoGo(HH_F);
                HH_F = HH_F+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 6
                
                % Neutral Fake Go(Stop) vs. Neutral Fake NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = Neutral_Fake_Go(NN_F);
                stimnum2(numRun,trial) = Neutral_Fake_NoGo(NN_F);
                NN_F = NN_F+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
        end % end switch pairtype
        
        fprintf(fid1, '%d\t%d\t%d\t%d\t%s\t%s\t\n', stimnum1(numRun,trial),stimnum2(numRun,trial),leftGo(numRun,trial),pairType(numRun,trial),leftname{numRun,trial},rightname{numRun,trial});
    end % end for trial = 1:total_num_trials
    
    fprintf(fid1, '\n');
    fclose(fid1);
end % end for numRun = 1:numRunsPerBlocks


%---------------------------------------------------------------------
% create a data structure with info about the run and all the matrices
%---------------------------------------------------------------------
outfile = strcat('./Output/', sprintf('%s_stimuliForProbe_order%d_block_%d_%d_trials_%d_runs_%s.mat',subjectID,order,block,total_num_trials,numRunsPerBlock,date));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
run_info.script_name = mfilename;

save(outfile);


end % end function

