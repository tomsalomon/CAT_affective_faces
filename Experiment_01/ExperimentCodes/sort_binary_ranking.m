function [] = sort_binary_ranking(subjectID,order,outputPath)
% function [] = sort_BDM_Israel(subjectID,order,outputPath)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik November 2014 ====================
% ================= Modified by Tom Salomon, January 2015 =================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function sorts the stimuli according to the Binary ranking results

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   [subjectID,'_ItemRankingResults*']


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order%d.txt', order
%   'stopGoList_trainingstim.txt'
%   'order_testcomp.txt'

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'BM_9001';
% order = 1;
% outputPath = 'D:/Rotem/Matlab/Boost_Israel_New_Rotem/Output'; % On Rotem's PC


%=========================================================================
%%  read in info from BDM1.txt
%=========================================================================

subj_happy_rank_file=dir(fullfile(outputPath,[subjectID,'_Happy_ItemRankingResults*']));
fid=fopen(fullfile(outputPath,subj_happy_rank_file.name));
Happy_data=textscan(fid, '%s %s %f %f %f %f %f', 'HeaderLines', 1); %read in data as new matrix
fclose(fid);

subj_hetral_rank_file=dir(fullfile(outputPath,[subjectID,'_Neutral_ItemRankingResults*']));
fid1=fopen(fullfile(outputPath,subj_hetral_rank_file.name));
Neutral_data=textscan(fid1, '%s %s %f %f %f %f %f', 'HeaderLines', 1); %read in data as new matrix
fclose(fid1);

%=========================================================================
%%  Create matrix sorted by descending bid value
%========================================================================

[bids_sort1,trialnum_sort_bybid1]=sort(Happy_data{4},'descend');
[bids_sort2,trialnum_sort_bybid2]=sort(Neutral_data{4},'descend');
bids_sort = [bids_sort1; bids_sort2];
trialnum_sort_bybid = [trialnum_sort_bybid1; trialnum_sort_bybid2];

stimnames_sorted_by_bid = [Happy_data{2}(trialnum_sort_bybid(1:40,1)); Neutral_data{2}(trialnum_sort_bybid(41:80,1))];


bid_sortedM(:,1)=trialnum_sort_bybid; %trialnums organized by descending bid amt
bid_sortedM(:,2)=bids_sort; %bids sorted large to small
bid_sortedM(:,3)=1:1:80; %stimrank


%=========================================================================
%%   The ranking of the stimuli determine the stimtype
%=========================================================================
if order == 1
    
    bid_sortedM([          7 10 12 13 15 18           ], 4) = 11; % Happy_beep
    bid_sortedM([ 3:4  8 9 11 14 16 17   19:22  37:38 ], 4) = 12; % Happy_nobeep
    bid_sortedM([          47 50 52 53 55 58          ], 4) = 22; % Neutral_beep
    bid_sortedM([ 43:44  48 49 51 54 56 57 59:62 77:78], 4) = 24; % Neutral_nobeep
    bid_sortedM([1:2  5:6  23:36 39:42 45:46  63:76  79:80], 4) = 0;  % notTrained
    
else
    
    bid_sortedM([          8 9 11 14 16 17            ], 4) = 11; % Happy_beep
    bid_sortedM([ 3:4  7 10 12 13 15 18  19:22  37:38 ], 4) = 12; % Happy_nobeep
    bid_sortedM([          48 49 51 54 56 57          ], 4) = 22; % Neutral_beep
    bid_sortedM([ 43:44  47 50 52 53 55 58 59:62 77:78], 4) = 24; % Neutral_nobeep
    bid_sortedM([1:2  5:6  23:36 39:42 45:46 63:76  79:80], 4) = 0;  % notTrained
    
end

itemsForTraining = bid_sortedM(bid_sortedM(:,4)~=0,:);
itemsNamesForTraining = stimnames_sorted_by_bid(bid_sortedM(:,4)~=0);

%=========================================================================
%%  create stopGoList_allstim.txt
%   this file is used during probe
%=========================================================================

fid2 = fopen([outputPath '/' subjectID sprintf('_stopGoList_allstim_order%d.txt', order)], 'w');


for i = 1:length(bid_sortedM)
    fprintf(fid2, '%s\t%d\t%d\t%d\t%d\t\n', stimnames_sorted_by_bid{i},bid_sortedM(i,4),bid_sortedM(i,3),bid_sortedM(i,2),bid_sortedM(i,1));
end

fprintf(fid2, '\n');
fclose(fid2);

fid3 = fopen([outputPath '/' subjectID '_stopGoList_trainingstim.txt'], 'w');

for i = 1:length(itemsForTraining)
    fprintf(fid3, '%s\t%d\t%d\t%d\t%d\t\n', itemsNamesForTraining{i},itemsForTraining(i,4),itemsForTraining(i,3),itemsForTraining(i,2),itemsForTraining(i,1));
end
fprintf(fid3, '\n');
fclose(fid3);

end % end function

