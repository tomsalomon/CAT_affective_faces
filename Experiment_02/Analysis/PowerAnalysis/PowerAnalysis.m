clear;
close all;
tic;

% Define data
DataPath = './../../Output/';
invalid_subjects = [308, 317, 332, 302, 305, 327, 347, 132];
power = 0.8;
alpha = 0.05;
min_sample_size = 5;
max_sample_size = 150;
num_of_tests_per_sample_size = 1000; % change to 1000 later
consequtive_powered_n=5; % how many times in a row you want to get a powered sample before stopping

% Decide on which variable you want to calculate power.
% Model I: PHQ (rank) ~ PreferenceBias (Rank)
% Model II: PHQ (rank) ~ PreferenceBias (Rank) + CAT_effect (Rank)
% 1 - PreferenceBias in model I; 2 - CAT_effect; 3 -
% PreferenceBias beta; 4 - CAT_effect beta; 5 - both beta values
tested_variable=5;

%%% --------------------- Colley Ranking Analysis --------------------- %%%

% finding all text files that contain "binary_ranking" in their names:
binary_ranking_files = dir([DataPath,'*binary_ranking*.txt']);
% delete hidden files (assuming all hidden files start with '.')
binary_ranking_files = binary_ranking_files(arrayfun(@(x) ~strcmp(x.name(1),'.'),binary_ranking_files));

% Using Colley ranking to calculate the subjects' stim preferences:
% Generating a cell that will keep the colley scores:
colley_scores = cell(1,numel(binary_ranking_files));

% processing the colley algorithm on each "binary_ranking" file
for file = 1:numel(binary_ranking_files)
    % reading the specific subject's file:
    binary_ranking_fileID = fopen([binary_ranking_files(file).folder,'/',binary_ranking_files(file).name]);
    % generating the correct field types and ignoring the first row (of titles)
    binary_ranking_subject_file = textscan(binary_ranking_fileID,'%s %d %n %s %s %d %d %c %d %n','headerlines',1);
    fclose(binary_ranking_fileID);
    
    stim_left = binary_ranking_subject_file{6};
    stim_right = binary_ranking_subject_file{7};
    num_of_stimuli = length(unique([binary_ranking_subject_file{6};binary_ranking_subject_file{7}]));
    chose_left = binary_ranking_subject_file{8} == 'u';
    chose_right = binary_ranking_subject_file{8} == 'i';
    number_of_competitions = length(binary_ranking_subject_file{6});
    
    % generating and filling in the matrix of results (n*3):
    % 1. wins; 2. losses (not including one competition with response 'x');
    % 3. total competitions
    colley_mat_results = zeros(num_of_stimuli,3);
    for stim = 1:num_of_stimuli
        colley_mat_results(stim,1) = sum(stim_left == stim & chose_left) + sum(stim_right == stim & chose_right);
        colley_mat_results(stim,2) = sum(stim_left == stim & chose_right) + sum(stim_right == stim & chose_left);
        colley_mat_results(stim,3) = sum(colley_mat_results(stim,1)) + sum(colley_mat_results(stim,2));
    end
    
    % generating and filling in the symetric matrix of competitions (n*n):
    colley_mat_competitions = zeros(num_of_stimuli);
    for comp = 1:number_of_competitions
        if chose_right(comp) || chose_left(comp)
            colley_mat_competitions(stim_left(comp),stim_right(comp)) = -1;
            colley_mat_competitions(stim_right(comp),stim_left(comp)) = -1;
        end
    end
    
    % calling the colley function to rank all stims
    colley_scores{file} = colley(colley_mat_competitions,colley_mat_results);
end

% generating a matrix from colley_scores cell array
colley_scores_mat = cell2mat(colley_scores);

% calculating the happy stim and the neutral stim means per subject
happy_mean = mean(colley_scores_mat((1:40),:))';
neutral_mean = mean(colley_scores_mat((41:80),:))';

%%% --------------------- Happy Bias Analysis --------------------- %%%

% Happy_mean analysis:
% generating and filling in a matrix of happy means and PHQ scores:
% 1. subject numbers; 2. subject happy mean; 3. subjects' PHQ scores:
Summary_matrix = zeros(numel(binary_ranking_files),5);
% filling the happy_mean_mat's first column with subject numbers:
% generating a cell of file names:
binary_ranking_file_names = {binary_ranking_files.name}';

% extracting the numbers from the file names:
subject_numbers = cell(numel(binary_ranking_file_names),2);
for name = 1:numel(binary_ranking_file_names)
    underline_locations = strfind(binary_ranking_file_names{name},'_');
    subject_numbers{name,1} = str2double(binary_ranking_file_names{name}(underline_locations(1)+1:underline_locations(2)-1));
    subject_numbers{name,2} = subject_numbers{name,1};
    % changing file names that start with 100 to start with 300:
    if subject_numbers{name,1} < 300
        subject_numbers{name,2} = subject_numbers{name,1} + 200;
    end
    % detecting subject number from files with a different prefix:
    if isnan(subject_numbers{name,1})
        subject_numbers{name,1} = str2double(binary_ranking_file_names{name}(underline_locations(2)+1:underline_locations(3)-1));
        subject_numbers{name,2} = subject_numbers{name,1};
    end
end
% filling in subject numbers in happy_mean_mat's 1st column:
Summary_matrix(:,1) = cell2mat(subject_numbers(:,2));
% filling in happy mean in happy_mean_mat's 2nd column:
Summary_matrix(:,3) = happy_mean;
% filing in PHQ scores in happy_mean_mat's 3rd column:
% reading the PHQ_scores file:
binary_ranking_fileID = fopen([DataPath,'PHQ_scores.txt']);
% generating the correct field types and ignoring the first row (of titles)
PHQ_scores = textscan(binary_ranking_fileID,'%s %n','headerlines',1);
fclose(binary_ranking_fileID);
% Filling in the phq scores in happy_mean_mat's 3rd column:
for subject = 1:length(Summary_matrix)
    % Finding the phq score by subject number and copying to happy_mean_mat:
    subject_num_location = contains(PHQ_scores{1},num2str(Summary_matrix(subject,1)));
    Summary_matrix(subject,2) = PHQ_scores{2}(subject_num_location);
end

% Analysis of happy_mean and phq results:
% Defining valid and non-valid subjects:
valid_subject = (~(ismember(Summary_matrix(:,1),invalid_subjects)));

% correlation and scattering:
% Spearman correlation and scattering of happy_mean_score and phq_score
phq_happy_mean_spearman_corr_test = corr(Summary_matrix(valid_subject,2),Summary_matrix(valid_subject,3),'type','Spearman');
figure('Name','Spearman correlation - PHQ and Preference Bias (mean happy rank)');
scatter(tiedrank(Summary_matrix(valid_subject,2)),tiedrank(Summary_matrix(valid_subject,3)));
xlabel('Preferecne bias - Happy mean Colley score (Rank)');
ylabel('PHQ (Rank)');
pause(0.5);
% Dependent t-test:
% Generating a matrix of pair diffs:
pair_diffs_mat = colley_scores_mat(1:40,:) - colley_scores_mat(41:num_of_stimuli,:);
% "Dependent t test" = mean/ (std/sqrt(n)):
Summary_matrix(:,4) = mean(pair_diffs_mat)'./ (std(pair_diffs_mat)' /sqrt(num_of_stimuli/2));

% Spearman correlation and scattering of t test (dependent) and phq score:
figure('Name','Spearman correlation - PHQ and Preference Bias ("dependent t")');
[phq_t_dependent_spearman_corr,phq_t_dependent_spearman_p] = corr(Summary_matrix(valid_subject,4),Summary_matrix(valid_subject,2),'type','Spearman');
scatter(tiedrank(Summary_matrix(valid_subject,4)),tiedrank(Summary_matrix(valid_subject,2)));
xlabel('Preferecne bias - "Dependent t" (Rank)');
ylabel('PHQ (Rank)');
pause(0.5);

%%% --------------------- Power Analysis --------------------- %%%

% analyzing probe restuls of all subjects using Probe_analysis script:
probe_results = Probe_analysis();
% fix miss-coded subjects (100's turn to 300's)
probe_results(probe_results(:,1)<=300,1)=probe_results(probe_results(:,1)<=300,1)+200;
% Copying the CAT_effect (Happy faces GO vs NoGo - Percent chosen Go) from probe_results to CAT_effect_mat:
[~,temp_location] = ismember(Summary_matrix(:,1),probe_results(:,1));
Summary_matrix(temp_location~=0,5)=probe_results(temp_location(temp_location~=0),7);
% convert matrix to table with titles
Summary_table = array2table(Summary_matrix(valid_subject,:));
Summary_table.Properties.VariableNames = {'SubjectID','PHQ','PreferenceBiasMean','PreferenceBias','CAT_effect'};

Summary_table.PHQ_rank=tiedrank(Summary_table.PHQ);
Summary_table.PreferenceBias_rank=tiedrank(Summary_table.PreferenceBias);
Summary_table.CAT_effect_rank=tiedrank(Summary_table.CAT_effect);
my_model = fitlm(Summary_table,'PHQ_rank ~ PreferenceBias_rank + CAT_effect_rank');

Summary_table.z_PHQ_rank=zscore(Summary_table.PHQ_rank);
Summary_table.z_PreferenceBias_rank=zscore(Summary_table.PreferenceBias_rank);
Summary_table.z_CAT_effect_rank=zscore(Summary_table.CAT_effect_rank);
my_model_z = fitlm(Summary_table,'z_PHQ_rank ~ z_PreferenceBias_rank + z_CAT_effect_rank');

figure('Name','GLM PHQ ~ PreferenceBias + CAT_effect')
plot(my_model)
my_p_values=my_model.Coefficients.pValue;
text(min(xlim)+0.1*(max(xlim)-min(xlim)),max(ylim)*0.9,sprintf('PHQ: p = %.3f \nPreferenceBias: p = %.3f',my_p_values(2),my_p_values(3)))
pause(0.5);

% power analysis - analyzing the minimal number of subjects required in
% order to receive a significant effect in 5 different cases:
% 1) PHQ ~ PreferenceBias
% 2) PHQ ~ CAT
% 3) PHQ ~ PreferenceBias + CAT (PreferenceBias beta)
% 4) PHQ ~ PreferenceBias + CAT (CAT beta)
% 5) PHQ ~ PreferenceBias + CAT (Both beta)
% power analysis - analyzing the minimal number of subjects required in
% power_tests is a cell array that contains all power tests:
power_tests = nan(1,4);
rand_sample = cell(1);
power_per_sample_size_mat = zeros(1,6);

row = 0;
sample_size_ind=0;
progress_bar = waitbar(0);
figure('Name','Power Analysis');
for sample_size = min_sample_size:max_sample_size
    rng(1); % for reproducibility
    sample_size_ind=sample_size_ind+1;
    for test = 1:num_of_tests_per_sample_size
        row = row + 1;
        waitbar(test/num_of_tests_per_sample_size,progress_bar,sprintf('Analyzing n = %i, %.1f%%',sample_size,100*test/num_of_tests_per_sample_size));
        rand_sample{row} = datasample(Summary_table{:,1},sample_size);
        % generating a temp mat for the current test:
        % 1. test's subjects; 2. dependent t-test result; 3. CAT effect; 4. phq_score;
        temp_power_test_mat = zeros(sample_size,4);
        temp_power_test_mat(:,1) = cell2mat(rand_sample(row));
        % copying the dependent_t_test, CAT_effect and phq_score of current test samples from dependent_t_test_mat:
        [~,temp_test] = (ismember(temp_power_test_mat(:,1),Summary_matrix(:,1)));
        temp_power_test_mat(:,2:4) = tiedrank(Summary_matrix(temp_test,[2,4:5]));
        
        temp_power_test_spearman_table = array2table(temp_power_test_mat);
        temp_power_test_spearman_table.Properties.VariableNames = {'subject_number','PHQ','PreferenceBias','CAT_effect'};
        
        % calculating p value for each lm in each test
        
        % First model: PHQ ~ dependent t-test bias
        lm1 = fitlm(temp_power_test_spearman_table,'PHQ ~ PreferenceBias');
        power_tests(row,1)=lm1.Coefficients.pValue(2); % p-value of PreferenceBias
        lm2 = fitlm(temp_power_test_spearman_table,'PHQ ~ CAT_effect');
        power_tests(row,2)=lm2.Coefficients.pValue(2); % p-value of PreferenceBias
        % Second model: PHQ ~ dependent t-test bias + CAT_effect
        lm3 = fitlm(temp_power_test_spearman_table,'PHQ~PreferenceBias+CAT_effect');
        power_tests(row,3) = lm3.Coefficients.pValue(2); % p-value of PreferenceBias
        power_tests(row,4) = lm3.Coefficients.pValue(3); % p-value of CAT_effect
        
    end
    
    power_tests(:,5)=(power_tests(:,3)>alpha) | (power_tests(:,4)>alpha); % both betas are significant
    current_sample_size_ind=(1+(sample_size_ind-1)*num_of_tests_per_sample_size):(sample_size_ind)*num_of_tests_per_sample_size;
    current_significance_count = sum(power_tests(current_sample_size_ind,:)<alpha); %CHANGE
    
    % power_per_sample_size_mat: 1. sample size; 2. power
    % calculating the power of each sample size:
    power_per_sample_size_mat(sample_size_ind,1) = sample_size;
    power_per_sample_size_mat(sample_size_ind,2:end) = current_significance_count / num_of_tests_per_sample_size; % power
    % Calculating the number of consecutive sample_sizes with >=80% significant results:
    
    
    
    % plotting the power as a function of sample size
    
    plot(power_per_sample_size_mat(:,1),power_per_sample_size_mat(:,2:end));
    xlim([0,max_sample_size]);
    ylim([0,1])
    xticks(0:10:max_sample_size);
    xlabel('Sample Size')
    ylabel('Power')
    % plot reference required power
    hold on
    plot(xlim,power*ones(size(xlim)),'--k')
    hold off
    legend({'Preference Bias','CAT','\beta Preference Bias','\beta CAT','Both \beta',sprintf('%i%% power',power*100)})
    pause(0.1);
    
    % Stop once you hit the required number of consequtive samples sizes with required power
    try % will not work until you completed a minimal number of sample sized
        stop_condition = power_per_sample_size_mat(end - consequtive_powered_n +1:end,tested_variable+1) > power ;
        if all(stop_condition)
            break
        end
    catch
    end
end

% Writing the N (minimal number of subjects required in order to receive a
% significant effect) on the graph:
num_of_variables=length(power_per_sample_size_mat(1,:)) -1;
minimal_n=nan(num_of_variables,1);
ax = gca;
for test_variable = 1:num_of_variables
    try
        minimal_n(test_variable)=power_per_sample_size_mat(find(power_per_sample_size_mat(:,test_variable+1)>power,1,'first'),1);
        text(minimal_n(test_variable)+2,0.1*test_variable,sprintf('N = %.0f',minimal_n(test_variable)),'FontWeight','bold','color',ax.ColorOrder(test_variable,:))
    catch
    end
end
hold on
ax.ColorOrderIndex = 1;
plot([minimal_n,minimal_n],ylim,'--');
hold off
legend({'Preference Bias','CAT','\beta Preference Bias','\beta CAT','Both \beta',sprintf('%i%% power',power*100)},'Location','southeast')
xlim([0,10*(1+ceil(sample_size/10))])

AnalysisDuration=toc;
timestamp=clock;
save(sprintf('PowerAnalysisResults_%04.f%02.f%02.f_%02.f%02.f.mat',timestamp(1),timestamp(2),timestamp(3),timestamp(4),timestamp(5)))

