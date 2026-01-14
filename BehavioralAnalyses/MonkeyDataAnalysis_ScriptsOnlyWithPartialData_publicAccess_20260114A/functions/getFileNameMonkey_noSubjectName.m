function [fileName_fig, shortfileName_fig] = getFileNameMonkey_noSubjectName(currentFigName, monkey_name)
% to generate a unique file name for saving figure

temp_fileCount = 1;

nameCell_fig = cell(1,4);
nameCell_fig{1} = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\';
nameCell_fig{2} = [monkey_name, '_'];
nameCell_fig{3} = currentFigName;
nameCell_fig{4} = num2str(temp_fileCount);
nameCell_fig{5} = '.emf';
% nameCell_fig{5} = '.tif';
fileName_fig = [nameCell_fig{1} nameCell_fig{2} nameCell_fig{3} nameCell_fig{4} nameCell_fig{5}];

while exist(fileName_fig , 'file') == 2
    temp_fileCount = temp_fileCount + 1;
    nameCell_fig{4} = num2str(temp_fileCount);    
    
    fileName_fig = [nameCell_fig{1} nameCell_fig{2} nameCell_fig{3} nameCell_fig{4} nameCell_fig{5}];
end
shortfileName_fig = [nameCell_fig{2} nameCell_fig{3} nameCell_fig{4} nameCell_fig{5}];
end
