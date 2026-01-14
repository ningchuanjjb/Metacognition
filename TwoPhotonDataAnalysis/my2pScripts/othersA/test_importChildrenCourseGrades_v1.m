filename_excel = 'D:\OneDrive\Inner Data\Children\成绩单&名单&游戏得分表\学生成绩单汇总_20231007A_onlyForMatlab.xlsx';
table_Class2020_Course2021S2 = readtable(filename_excel,'sheet',1,'VariableNamingRule','preserve');
table_Class2021_Course2021S1 = readtable(filename_excel,'sheet',2,'VariableNamingRule','preserve');
table_Class2021_Course2021S2 = readtable(filename_excel,'sheet',3,'VariableNamingRule','preserve');

matrix_Class2020_Course2021S2 = str2double(table2array(table_Class2020_Course2021S2(:,4:end)));
matrix_Class2021_Course2021S1 = str2double(table2array(table_Class2021_Course2021S1(:,4:end)));
matrix_Class2021_Course2021S2 = str2double(table2array(table_Class2021_Course2021S2(:,4:end)));

courseName_Class2020_Course2021S2 = table_Class2020_Course2021S2.Properties.VariableNames(4:end);
courseName_Class2021_Course2021S1 = table_Class2021_Course2021S1.Properties.VariableNames(4:end);
courseName_Class2021_Course2021S2 = table_Class2021_Course2021S2.Properties.VariableNames(4:end);

a = 1;

%% childID_Class2020_Course2021S2
temp_childIDA = table2array(table_Class2020_Course2021S2(:,3));
temp_childIDB = num2str(table2array(table_Class2020_Course2021S2(:,2)));
temp_num = size(temp_childIDA,1);
childID_Class2020_Course2021S2 = cell(temp_num,1);
for tempi=1:temp_num
    tempStrA = '2020';    
    tempStrB = temp_childIDA{tempi}(2);
    tempStrC = temp_childIDB(tempi,:);
    if strcmp(tempStrC(1),' ') == 1
        tempStrC(1) = '0';
    end    
    tempID = [tempStrA,tempStrB,tempStrC];
    childID_Class2020_Course2021S2{tempi} = tempID;
end

%% childID_Class2021_Course2021S1
temp_childIDA = table2array(table_Class2021_Course2021S1(:,3));
temp_childIDB = num2str(table2array(table_Class2021_Course2021S1(:,2)));
temp_num = size(temp_childIDA,1);
childID_Class2021_Course2021S1 = cell(temp_num,1);
for tempi=1:temp_num
    tempStrA = '2021';    
    tempStrB = temp_childIDA{tempi}(2);
    tempStrC = temp_childIDB(tempi,:);
    if strcmp(tempStrC(1),' ') == 1
        tempStrC(1) = '0';
    end    
    tempID = [tempStrA,tempStrB,tempStrC];
    childID_Class2021_Course2021S1{tempi} = tempID;
end

%% childID_Class2021_Course2021S2
temp_childIDA = table2array(table_Class2021_Course2021S2(:,3));
temp_childIDB = num2str(table2array(table_Class2021_Course2021S2(:,2)));
temp_num = size(temp_childIDA,1);
childID_Class2021_Course2021S2 = cell(temp_num,1);
for tempi=1:temp_num
    tempStrA = '2021';    
    tempStrB = temp_childIDA{tempi}(2);
    tempStrC = temp_childIDB(tempi,:);
    if strcmp(tempStrC(1),' ') == 1
        tempStrC(1) = '0';
    end    
    tempID = [tempStrA,tempStrB,tempStrC];
    childID_Class2021_Course2021S2{tempi} = tempID;
end

%%
table_Class2020_Course2021S2;
table_Class2021_Course2021S1;
table_Class2021_Course2021S2;

matrix_Class2020_Course2021S2;
matrix_Class2021_Course2021S1;
matrix_Class2021_Course2021S2;

courseName_Class2020_Course2021S2;
courseName_Class2021_Course2021S1;
courseName_Class2021_Course2021S2;

childID_Class2020_Course2021S2;
childID_Class2021_Course2021S1;
childID_Class2021_Course2021S2;

Class2020_Course2021S2 = struct;
Class2020_Course2021S2.table = table_Class2020_Course2021S2;
Class2020_Course2021S2.matrix = matrix_Class2020_Course2021S2;
Class2020_Course2021S2.courseName = courseName_Class2020_Course2021S2;
Class2020_Course2021S2.childID = childID_Class2020_Course2021S2;

Class2021_Course2021S1 = struct;
Class2021_Course2021S1.table = table_Class2021_Course2021S1;
Class2021_Course2021S1.matrix = matrix_Class2021_Course2021S1;
Class2021_Course2021S1.courseName = courseName_Class2021_Course2021S1;
Class2021_Course2021S1.childID = childID_Class2021_Course2021S1;

Class2021_Course2021S2 = struct;
Class2021_Course2021S2.table = table_Class2021_Course2021S2;
Class2021_Course2021S2.matrix = matrix_Class2021_Course2021S2;
Class2021_Course2021S2.courseName = courseName_Class2021_Course2021S2;
Class2021_Course2021S2.childID = childID_Class2021_Course2021S2;

a = 1;


