function [] = Copy2Dropbox(subjectID,outputPath,targetPath)

% Define defualt variables
if nargin<3
    targetPath = '~/Dropbox/experimentsOutput/TomSalomon/Boost_faces_emotional/';
    if nargin<2
        outputPath = './Output/';
    end
end

copyfile ([outputPath,subjectID,'*'],targetPath);
files = dir([targetPath,subjectID,'*']);
fprintf('\nSuccessfully copied %i files for subject: %s\n',length(files),subjectID)