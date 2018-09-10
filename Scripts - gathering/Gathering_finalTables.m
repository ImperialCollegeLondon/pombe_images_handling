%% COMMON STEPS
input0  = '..\..\..\Cycle sensor data\2. Dati sperimentali\';
%addpath('..\..\CellSegment plotting');
newTitles = [                                                                                   %OLD
    {'cell_area'},    {'cell_length'},  {'cell_width_avg'}, {'cell_width_max'}...   %1:4        1:3
    {'dye_area'},     {'dye_length'},   {'dye_width'},      {'dye_delta'},...       %5:8        4:7
    {'ratio_area'},   {'ratio_length'}, {'ratio_width'},...                         %9:11       8:10
    {'DNA_intensity'},{'DNA_perNuc'},   {'nuclei_delta'},...                        %12:14      11:13
    {'nucleiNr'},     {'mid2_check'},   {'cut3_check'},     {'cycle_phase'},...     %15:18      14:17
    {'growth_rate'},  {'OD_extraction'},{'stain_check'}, ...                        %19:21      18:20
    {'analysis_mode'},{'condition'},    {'replica'},        {'tPoint'}              %22:25      21:24
];



%% MATLAB DYE
run('..\Scripts - auxiliary\outdated_dye_general.m');
clearvars -except input0 newTitles conds16 OD_dye featDye featNamDye grVec16 indSize namDye %qFlagDye wrongDye
cmlSize = [0 cumsum(indSize)];

output1 = [input0 '2. Dati finali\MAT_DYE\'];  mkdir(output1)
for aa=1:16
    matEND  = cell(0,25);
    nameEND = cell2mat(conds16(aa));
    nameEND = ['MAT_Dye_' strrep(nameEND(3:end),' ','_')];
    
    for bb=1:indSize(aa)
        cc = cmlSize(aa)+bb;
        nameIN  = cell2mat(namDye(cc));
        tmatIN  = cell2mat(featDye(cc2));
        tmatIN  = tmatIN(:,[1:13,end]);
        tmatOUT = cell(size(tmatIN,1),25);
        tmatOUT(:) = {NaN};
        tmatOUT(:,[2,3,1, 15,14,16, 6,5,8, 9,10,12, 13,4]) = num2cell(tmatIN);
        tmatOUT(tmatIN(:,6)<=0,16) = {0};    %new entry
        tmatOUT(tmatIN(:,6)>=1,16) = {1};    %new entry
        tmatOUT(tmatIN(:,4)==1,18) = {1};
        tmatOUT(tmatIN(:,4)==2 & tmatIN(:,6)<=0,18) = {3};
        tmatOUT(tmatIN(:,4)==2 & tmatIN(:,6)>=1,18) = {4};

        tmatOUT(:,7)  = {0};                %by design
        tmatOUT(:,11) = {0};                %by design
        tmatOUT(tmatIN(:,7)~=0,7) = num2cell(tmatIN(tmatIN(:,7)~=0,2));
        tmatOUT(tmatIN(:,7)~=0,11) = {1};
        tmatOUT(:,19) = {grVec16(aa)};
        tmatOUT(:,20) = {OD_dye(cc)};
        tmatOUT(:,21) = {1};
        tmatOUT(:,22) = {'MAT'};
        tmatOUT(:,23) = {nameIN(1:3)};
        tmatOUT(:,24) = {nameIN(5)};
        tmatOUT(:,25) = {nameIN(7:end)};

        %converting pixel to micron
        tmatOUT(:,1)    = cellfun(@(x) x/9,tmatOUT(:,1),'un',0); %cell area
        tmatOUT(:,2)    = cellfun(@(x) x/3,tmatOUT(:,2),'un',0); %cell length
        tmatOUT(:,3)    = cellfun(@(x) x/3,tmatOUT(:,3),'un',0); %cell width_avg
        tmatOUT(:,4)    = cellfun(@(x) x/3,tmatOUT(:,4),'un',0); %cell width_max
        tmatOUT(:,5)    = cellfun(@(x) x/9,tmatOUT(:,5),'un',0); %dye area
        tmatOUT(:,6)    = cellfun(@(x) x/3,tmatOUT(:,6),'un',0); %dye length
        tmatOUT(:,7)    = cellfun(@(x) x/3,tmatOUT(:,7),'un',0); %dye_width
        tmatOUT(:,8)    = cellfun(@(x) x/3,tmatOUT(:,8),'un',0); %dye_delta
        tmatOUT(:,14)   = cellfun(@(x) x/3,tmatOUT(:,14),'un',0);%nuclei_delta
        
        matEND = [matEND;tmatOUT];
        disp([aa/16, bb/indSize(aa)])
    end
    tabEND = cell2table(matEND,'VariableNames',newTitles);
    writetable(tabEND,[output1 nameEND '.txt']);
end
clearvars -except input0 newTitles
    


%% MATLAB MICRO
run('..\Scripts - auxiliary\outdated_micro_all.m');
clearvars -except input0 newTitles      featNamMic grVecMic namMic
%featNamMic grVecMic namMic

input1  = [input0 '1. Dati parziali\MAT - dati finali\M_CSVfiles\'];
output1 = [input0 '2. Dati finali\MAT_MICRO\'];  mkdir(output1)
list1 = dir(input1);    list1(1:2) = [];
for aa=1:length(list1)
    nameIN  = cell2mat(namMic(aa));
    nameOUT = ['MAT_Mic' upper(nameIN(1)) '_' nameIN(7:end)];
    
    matIN   = readtable([input1 list1(aa).name],'ReadVariableName',false);
    matIN   = matIN(:,[1:9,end]);
    matIN   = table2array(matIN);
    matOUT  = cell(size(matIN,1),25);
    matOUT(:) = {NaN};    
    matOUT(:,[2,3,1,15,14,16,17,12,13,4]) = num2cell(matIN);
    matOUT(matIN(:,6)<=0,16) = {0};
    matOUT(matIN(:,6)>=1,16) = {1};
    matOUT(matIN(:,7)<=0,17) = {0};
    matOUT(matIN(:,7)>=1,17) = {1};  
    matOUT(matIN(:,4)==1 & matIN(:,7)<=0,18) = {1};
    matOUT(matIN(:,4)==1 & matIN(:,7)>=1,18) = {2};
    matOUT(matIN(:,4)==2 & matIN(:,6)<=0,18) = {3};
    matOUT(matIN(:,4)==2 & matIN(:,6)>=1,18) = {4};
    
    matOUT(:,19) = {grVecMic(aa)};
    matOUT(:,20) = {0.4};
    matOUT(:,21) = {0};
    matOUT(:,22) = {'MAT'};
    matOUT(:,23) = {nameIN(7:9)};
    matOUT(:,24) = {nameIN(11)};
    
    %converting pixel to micron
    matOUT(:,1)    = cellfun(@(x) x*0.0256,matOUT(:,1),'un',0); %cell area
    matOUT(:,2)    = cellfun(@(x) x*0.1600,matOUT(:,2),'un',0); %cell length
    matOUT(:,3)    = cellfun(@(x) x*0.1600,matOUT(:,3),'un',0); %cell width_avg
    matOUT(:,4)    = cellfun(@(x) x*0.1600,matOUT(:,4),'un',0); %cell width_max
    matOUT(:,14)   = cellfun(@(x) x*0.1600,matOUT(:,14),'un',0);%nuclei_delta
    
    tabOUT = cell2table(matOUT,'VariableNames',newTitles);
    writetable(tabOUT,[output1 nameOUT '.txt']);
end
clearvars -except input0 newTitles



%% MATLAB TURBI
run('..\Scripts - auxiliary\outdated_turbi_all.m');
clearvars -except featNamTur input0 newTitles  grVecTur namTur  

input1  = [input0 '1. Dati parziali\MAT - dati finali\T_CSVfiles\'];
output1 = [input0 '2. Dati finali\MAT_TURBI\'];  mkdir(output1)
list1 = dir(input1);    list1(1:2) = [];
namTur(17) = {'T Ser B'};  
for aa=1:3 %1:length(list1)
    nameIN  = cell2mat(namTur(aa));
    nameOUT = ['MAT_Tur_' nameIN(3:end)];
    
    matIN   = readtable([input1 list1(aa).name],'ReadVariableName',false);
    matIN   = matIN(:,[1:9,end]);
    matIN   = table2array(matIN);
    matOUT  = cell(size(matIN,1),25);
    matOUT(:) = {NaN};    
    matOUT(:,[2,3,1, 15,14,16, 17,12,13, 4]) = num2cell(matIN);
    matOUT(matIN(:,6)<=0,16) = {0};
    matOUT(matIN(:,6)>=1,16) = {1}; 
    matOUT(matIN(:,7)<=0,17) = {0};
    matOUT(matIN(:,7)>=1,17) = {1};
    matOUT(matIN(:,4)==1 & matIN(:,7)<=0,18) = {1};
    matOUT(matIN(:,4)==1 & matIN(:,7)>=1,18) = {2};
    matOUT(matIN(:,4)==2 & matIN(:,6)<=0,18) = {3};
    matOUT(matIN(:,4)==2 & matIN(:,6)>=1,18) = {4};

    matOUT(:,19) = {grVecTur(aa)};
    matOUT(:,20) = {0.4};
    matOUT(:,21) = {0};
    matOUT(:,22) = {'MAT'};
    matOUT(:,23) = {nameIN(3:5)};
    matOUT(:,24) = {nameIN(7)};

    %converting pixel to micron
    matOUT(:,1)    = cellfun(@(x) x/9,matOUT(:,1),'un',0); %cell area
    matOUT(:,2)    = cellfun(@(x) x/3,matOUT(:,2),'un',0); %cell length
    matOUT(:,3)    = cellfun(@(x) x/3,matOUT(:,3),'un',0); %cell width_avg
    matOUT(:,4)    = cellfun(@(x) x/3,matOUT(:,4),'un',0); %cell width_max
    matOUT(:,14)   = cellfun(@(x) x/3,matOUT(:,14),'un',0);%nuclei_delta

    tabOUT = cell2table(matOUT,'VariableNames',newTitles);
    writetable(tabOUT,[output1 nameOUT '.txt']);
end
clearvars -except input0 newTitles



%% ISX GENERAL
run('..\Scripts - auxiliary\outdated_dye_general.m');
run('..\Scripts - auxiliary\outdated_turbi_all.m');
clearvars -except input0 newTitles grVecTur conds16 grVec16 OD_DI
%grVecTur(11) =[];   grVecTur = repelem(grVecTur,2,1);
grVecTur(12) =[];   grVecTur = repelem(grVecTur,2,1);   %new, after swapping

input1  = [input0 '1. Dati parziali\ISX - dati finali\'];
list1   = dir(input1);    list1(1:2) = [];
for aa=1:length(list1)
    tName = list1(aa).name;
    tName= strrep(tName,'157_','');
    tName= strrep(tName,'NH4Cl','NH4');
    tName= strrep(tName,'r1_t', 'DyeF_A_t');
    tName= strrep(tName,'r10_t','DyeF_C_t');
    tName= strrep(tName,'_fix_r1',  '_TurF_A____');
    tName= strrep(tName,'_unfix_r1','_TurU_A____');
    tName= strrep(tName,'_fix_r2',  '_TurF_B____');
    tName= strrep(tName,'_unfix_r2','_TurU_B____');
    tName= strrep(tName,'_fix_r3',  '_TurF_C____');
    tName= strrep(tName,'_unfix_r3','_TurU_C____');
    tName= tName(1:14);
    fNames(aa) = {tName};
end
dyeIndex = cellfun(@isempty,strfind(fNames,'Tur')); %#ok<STRCLFH>
listDye = list1(dyeIndex);  listTur = list1(~dyeIndex);
nameDye = fNames(dyeIndex); nameTur = fNames(~dyeIndex);



%% ISX TURBI
output1 = [input0 '2. Dati finali\ISX_TURBI\'];  mkdir(output1)
for aa=1:length(listTur)
    nameIN  = cell2mat(nameTur(aa));
    nameOUT = ['ISX_Tur' nameIN(8) '_' nameIN(1:3) ' ' nameIN(10)];
    
    matIN   = readtable([input1 listTur(aa).name],'ReadVariableName',false, 'Delimiter','\t');
    matIN   = table2array(matIN);
    matOUT  = cell(size(matIN,1),25);
    matOUT(:) = {NaN};
    matOUT(:,[1,2,4,12:17]) = num2cell(matIN(:,[2:6,8,9,11,13]));
    matOUT(matIN(:,9)<0, 15) = {0};    matOUT(matIN(:,9)>2, 15) = {2};
    matOUT(matIN(:,11)<1,16) = {0};    matOUT(matIN(:,11)>1,16) = {1};    matOUT(matIN(:,9)<2, 16) = {0}; 
    matOUT(matIN(:,13)<1,17) = {0};    matOUT(matIN(:,13)>1,17) = {1};  
    
    matOUT(matIN(:,9)==1 & matIN(:,13)<1,18)  = {1};
    matOUT(matIN(:,9)==1 & matIN(:,13)>=1,18) = {2};
    matOUT(matIN(:,9)==2 & matIN(:,11)<1,18)  = {3};
    matOUT(matIN(:,9)==2 & matIN(:,11)>=1,18) = {4};

    matOUT(:,19) = {grVecTur(aa)};
    matOUT(:,20) = {0.4};
    matOUT(:,21) = {0};
    matOUT(:,22) = {'ISX'};
    matOUT(:,23) = {nameIN(1:3)};
    matOUT(:,24) = {nameIN(10)};

    tabOUT = cell2table(matOUT,'VariableNames',newTitles);
    writetable(tabOUT,[output1 nameOUT '.txt']);
end
clearvars matIN matOUT nameIN nameOUT output1 tabOUT



%% ISX DYE
indSize = [5,5, 4,4, 4,4, 5,6, 4,4, 5,5, 5,5, 3,5];
cmlSize = [0 cumsum(indSize)];
grVec16 = reshape(grVec16,2,8);
grVec16([1,2],:) = grVec16([2,1],:);
grVec16 = reshape(grVec16,1,16);

output1 = [input0 '2. Dati finali\ISX_DYE\'];  mkdir(output1)
tConds16 = conds16([2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15]);
namDye(127) = {'Ser A t12'};    namDye(139) = {'Ser B t12'};
for aa=1:16
    matEND  = cell(0,25);
    nameEND = cell2mat(tConds16(aa));
    nameEND = ['ISX_Dye_' strrep(nameEND(3:end),' ','_')];
    
    for bb=1:indSize(aa)
        cc = cmlSize(aa)+bb;
        nameIN  = cell2mat(nameDye(cc));
        tmatIN  = readtable([input1 listDye(cc).name],'ReadVariableName',false, 'Delimiter','\t');
        tmatIN  = table2array(tmatIN);
        tmatOUT = cell(size(tmatIN,1),25);
        tmatOUT(:) = {NaN};        
        tmatOUT(:,[1,2,4:7,12:14,9:11,8,15,16]) = num2cell(tmatIN(:,[2:9,11:16,18]));
        tmatOUT(tmatIN(:,16)<0, 15) = {0};    tmatOUT(tmatIN(:,16)>2, 15) = {2};
        tmatOUT(tmatIN(:,18)<1,16) = {0};     tmatOUT(tmatIN(:,18)>1,16) = {1};    tmatOUT(tmatIN(:,16)<2, 16) = {0}; 
        
        tmatOUT(tmatIN(:,16)==1,18) = {1};
        tmatOUT(tmatIN(:,16)==2 & tmatIN(:,18)<1,18) = {3};
        tmatOUT(tmatIN(:,16)==2 & tmatIN(:,18)>=1,18) = {4};
        
        tmatOUT(:,19) = {grVec16(aa)};
        tmatOUT(:,20) = {OD_DI(cc)};
        tmatOUT(:,21) = {1};
        tmatOUT(tmatIN(:,20)<1,21) = {0};        
        tmatOUT(:,22) = {'ISX'};
        tmatOUT(:,23) = {nameIN(1:3)};
        tmatOUT(:,24) = {nameIN(10)};
        tmatOUT(:,25) = {nameIN(12:end)};
        matEND = [matEND;tmatOUT];
    end
    tabEND = cell2table(matEND,'VariableNames',newTitles);
    writetable(tabEND,[output1 nameEND '.txt']);
end