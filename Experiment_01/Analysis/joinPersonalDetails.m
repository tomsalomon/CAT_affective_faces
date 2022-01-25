function [personalDetailsAllSubjects] = joinPersonalDetails(experimentName,mainPath,Subjects)
% This function joins the personal details of the subjects to one matrix
% experimentName: the prefix of subjectID of this experiment (e.g.
% BMI_bs_40)
% mainPath: the path of the main folder of this experiment, in which there
% is the output folder
% Subjects: the number of subjects to be analyzed (e.g. 101:110)

if nargin < 3
    subjects_group1=[301,303,307,310,314:316,323,326,328,334:336,338,340,342,345,346,351:355] ;
    subjects_group2=[113,304,306,309,311:312,318:322,324:325,329:331,333,337,339,341,343,344,348:350,356:359];
    Subjects = sort([subjects_group1,subjects_group2]); % Define here your subjects' codes.
% exclude (group):
% 308 (1), 317 (2) - Training: false alarm
% 132/332 (2), 302 (1), 305 (1) - Training: minimal ladder; 132 - also Training: misses.
% 356 (PHQ=7), 359 (PHQ=6) - can't assign group
% 
% not really excluded
% 313 (2) - was coded as 113 - not excluded
% 332 (2) - was coded as 132. Training: minimal ladder (stopped hearing cue).
% 327 (2) - did not run subject with that code
% 347 (1) - did not complete the experiment
end

if nargin < 2
    experimentName = '*';
end

if nargin < 1
    mainPath = './../';
end

outputPath = [mainPath '/Output'];

personalDetailsAllSubjects = cell(length(Subjects)+1,9);
titles = {'subjectID','order','date','gender(1=female)','age','dominantHand(1=right)','height','weight','occupation'};
personalDetailsAllSubjects(1,:) = titles;

for subjectInd = 1:length(Subjects)
    subjectNum = Subjects(subjectInd);
    filename = strcat(outputPath,sprintf('/%s%d',experimentName,subjectNum));
    logs = dir(strcat(filename, '_personalDetails','*.txt')) ;
    fid = fopen(strcat(outputPath,'/',logs(end).name));
    Data = textscan(fid, '%s %d %s %d %d %d %d %d %s' , 'HeaderLines', 1); % read in personal details output file into Data ;
    fclose(fid);
    for ind = 1:9 % put the details inside the matrix for all subjects
        personalDetailsAllSubjects{subjectInd+1,ind} = Data{ind};
    end
end

[personalDetailsAllSubjects] = fixPersonalDetailsMatrix(personalDetailsAllSubjects);
is_female=cell2mat(personalDetailsAllSubjects(2:end,4))==1;
age=double(cell2mat(personalDetailsAllSubjects(2:end,5)));
fprintf('Number of females: %i, Prop. of females: %.4f\n',sum(is_female),mean(is_female))
fprintf('ages %i - %i, mean = %.2f, SD = %.2f\n',min(age),max(age),mean(age),std(age))
end % end main function

function [personalDetailsAllSubjects] = fixPersonalDetailsMatrix(personalDetailsAllSubjects)
for row = 2:size(personalDetailsAllSubjects,1)
    for column = 1:size(personalDetailsAllSubjects,2)
        if length(personalDetailsAllSubjects{row,column})>1
            personalDetailsAllSubjects{row,column} = personalDetailsAllSubjects{row,column}(1);
        end
    end
end % end function fixPersonalDetailsMatrix(personalDetailsAllSubjects)
end