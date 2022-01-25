
clear

OutputPath='./../Output/';

binary_ranking_files=dir([OutputPath,'*binary_ranking*.txt']);
for file_id=1:length(binary_ranking_files)
    fid=fopen([OutputPath,binary_ranking_files(file_id).name]);
    Data=textscan(fid,'%s%f%f%s%s%f%f%s%f%f','Headerlines',1);
    fid=fclose(fid);
    SubID=Data{1}{1};
    stimnames=[Data{4};Data{5}];
    
    Data{8}(strcmp(Data{8},'u'))={1}; % response: 1 for left
    Data{8}(strcmp(Data{8},'i'))={0}; % response: 0 for right
    Data{8}(strcmp(Data{8},'x'))={999}; % response: 999 for no response
    
    stim_left=Data{6};
    stim_right=Data{7};
    chose_left=cell2mat(Data{8});
    
    % removed no choices
    missed_trials=chose_left==999;
    stim_left(missed_trials)=[];
    stim_right(missed_trials)=[];
    chose_left(missed_trials)=[];
    stimnames([missed_trials;missed_trials])=[];
    
    num_stimuli=length(unique([stim_left,stim_right]));
    
    Colley_mat_competitions=zeros(num_stimuli);
    Colley_mat_competitions(sub2ind(size(Colley_mat_competitions), stim_left, stim_right))=-1;
    Colley_mat_competitions(sub2ind(size(Colley_mat_competitions), stim_right, stim_left))=-1;
    
    Colley_mat_scores=zeros(num_stimuli,3);
    stimname=cell(num_stimuli,1);
    for stim=1:num_stimuli
        Colley_mat_scores(stim,1)=sum( (stim_left==stim &chose_left) | (stim_right==stim &(~chose_left)));
        Colley_mat_scores(stim,2)=sum( (stim_left==stim &(~chose_left)) | (stim_right==stim &chose_left));
        stimname{stim}=cell2mat(unique(stimnames([stim_left==stim;stim_right==stim])));
    end
    Colley_mat_scores(:,3)=sum(Colley_mat_scores(:,1:2),2);
    Colley_ranks=colley(Colley_mat_competitions,Colley_mat_scores);
    
    fid2=fopen([OutputPath,SubID,'_ItemRankingResults_combined.txt'], 'w');
    fprintf(fid2,'Subject\tStimName\tStimNum\tRank\tWins\tLoses\tTotal\n');
    for j=1:num_stimuli
        fprintf(fid2,'%s\t%s\t%d\t%d\t%d\t%d\t%d\n', SubID, stimname{j}, j, Colley_ranks(j), Colley_mat_scores(j,1), Colley_mat_scores(j,2), Colley_mat_scores(j,3));
    end
    fclose(fid2);  
end
