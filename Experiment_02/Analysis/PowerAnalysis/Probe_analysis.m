function [probe_results] = Probe_analysis()
% ~~~ Script for analyzing probe results, modified for face experiment ~~~
% ~~~~~~~~~~~~~~~ Tom Salomon, February 2015  ~~~~~~~~~~~~~~
%
% In order to run the script, you must locate and run the script from within
% "Analysis" folder. The script uses external function, called
% "Probe_recode" which joins all probe output file into one matrix. Please
% make sure that function is also present in the analysis folder.
%
% Note: this script and "Probe_recode" function were written specifically
% for face stimuli in a specific numeric name format. For other stimuli, you need
% to modify the Probe_recode function first.
%
% Enjoy!

% clear


% all subjects
valid_subjects=[113,301,303:304,306:307,309:312,314:316,318:326,328:331,333:346,348:359]; % Define here your subjects' codes.
% group1=[301:303,305,307:308,310,314:316,323,326,328,334:336,338,340];
group2=[113,132,301:359];

%exclude (group):
% 302 (1), 305 (1) - Training: misses
% 308 (1), 317 (2) - Training: false alarm
% not really excluded
% 313 (2) - Missing data. 313 was run as 113 - definitly not
% excluded
% 332 (2) - was coded as 132 - Training: minimal ladder (stopped hearing
% cue).
% 327 (2) - did not run subject with that code

% ---Define which group to analyze---
% subjects=valid_subjects(ismember(valid_subjects,group1)); % analyze group1
subjects=valid_subjects(ismember(valid_subjects,group2)); % analyze group2

data_path=['./../../Output/']; % Analysis folder location
all_subjs_data{length(subjects)}={};
probe_results=zeros(length(subjects),10);
probe_results(:,1)=subjects;

for subjInd=1:length(subjects)
    
    data=Probe_recode(subjects(subjInd),data_path);
    % The Probe_recode function will join all present output file for each subject into a single matrix, these are the matrix columns:
    % 1-subjectID       2-scanner	 3-order        4-block         5-run       6-trial	 7-onsettime    8-ImageLeft	 9-ImageRight	10-bidIndexLeft
    % 11-bidIndexRight	12-IsleftGo	 13-Response    14-PairType     15-Outcome  16-RT	 17-bidLeft     18-bidRight
    
    all_subjs_data{subjInd}=data; %All subjects' data
    order=data(1,3);
    
    PairType=data(:,14);
    Outcome=data(:,15);
    Rank_left=data(:,10);
    Rank_right=data(:,11);
    
    % Organize data in a summary table
    probe_results(subjInd,2)=order;
    
    probe_results(subjInd,3)=sum(PairType==1&Outcome~=999); % Happy faces GO vs NoGo - number of valid trials
    probe_results(subjInd,4)=sum(PairType==2&Outcome~=999); % Neutral faces GO vs NoGo - number of valid trials
    probe_results(subjInd,5)=sum(PairType==4&Outcome~=999); % Neutral Sanity check - number of valid trials
    probe_results(subjInd,6)=sum(Outcome==999); % number of invalid trials
    
    probe_results(subjInd,7)=sum(PairType==1&Outcome==1)/sum(PairType==1&Outcome~=999); % Happy faces GO vs NoGo - Percent chosen Go
    probe_results(subjInd,8)=sum(PairType==2&Outcome==1)/sum(PairType==2&Outcome~=999); % Neutral faces GO vs NoGo - Percent chosen Go
    probe_results(subjInd,9)=sum(PairType==3&Outcome==1)/sum(PairType==3&Outcome~=999); % Happy Sanity check - Percent chosen Sanely
    probe_results(subjInd,10)=sum(PairType==4&Outcome==1)/sum(PairType==4&Outcome~=999); % Neutral Sanity check - Percent chosen Sanely
    
end

Probe_results_table = cell(1+size(probe_results,1),size(probe_results,2));
Titles = {'Subject', 'Order', '#Happy', '#Neutral', '#SanityHappy', '#InvalidTrials', '%HappyChoseGo', '%NeutralChoseGo', '%SanityHappyChoseHigh', '%SanityNeutralGoChoseHigh'};
Probe_results_table(1,:) = Titles;
Probe_results_table(2:end,:) = num2cell(probe_results);


% analyze the data for all subject
means=zeros(1,length(probe_results(1,:))-6);
stddevs=zeros(1,length(probe_results(1,:))-6);
p_values=zeros(1,length(probe_results(1,:))-6);

Text{1}='\nHappy faces GO vs NoGo';
Text{2}='Neutral faces GO vs NoGo';

Text{3}='\nHappy Sanity check';
Text{4}='Neutral Sanity check';

fprintf('\nProbe Results\n')
fprintf('=============\n')
for i=1:length(Text)
    means(i)=mean(probe_results(:,i+6));
    stddevs(i)=std(probe_results(:,i+6));
    [~,p_values(i)]=ttest(probe_results(:,i+6),0.5);
    fprintf([Text{i},': mean=%.2f, p=%.3f\n'],means(i),p_values(i));
end
stderr=stddevs.*(1/sqrt(length(subjects)));

end