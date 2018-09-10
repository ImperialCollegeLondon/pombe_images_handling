%% DYE
input0 = '..\Files\treated-data\Dye\';
mkdir([input0 'CSVfiles/']);
list0 = dir(input0); list0(1:3) = []; list0 = list0([list0(:).isdir]);
for aa=1:length(list0)
    input1  = [input0 list0(aa).name '/'];
    list1   = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
    for bb=1:length(list1)
        fileName= list1(bb).name;
        featStr = load([input1 fileName '/' fileName]);
        featStr = featStr.features;
        featMat = struct2array(featStr);
        csvwrite([input0 'CSVfiles/' fileName '.txt'],featMat);
        disp([aa/length(list0) bb/length(list1)]);
    end
end
clear



%% GROWTH
input0 = '..\Files\treated-data\Growth\';
mkdir([input0 'CSVfiles/']);
list0 = dir(input0); list0(1:2)=[]; list0 = list0([list0(:).isdir]);    list0(end)=[];
for aa=1:length(list0)
    fileName= list0(aa).name;
    featStr = load([input0 fileName '/' fileName]);
    featStr = featStr.features;
    featMat = struct2array(featStr);
    csvwrite([input0 'CSVfiles/' fileName '.txt'],featMat);
    disp([aa/length(list0)]);
end
clear



%% MICROSCOPY
input0 = '..\Files\treated-data\Microscopy\';
mkdir([input0 'CSVfiles/']);
list0 = dir(input0); list0(1:3) = []; list0 = list0([list0(:).isdir]);
for aa=1:length(list0)
    input1= list0(aa).name;
    list1 = dir([input0 input1]); list1(1:2) = []; list1 = list1(~[list1(:).isdir]);
    for bb=1:length(list1)
        fileName= list1(bb).name;
        featStr = load([input0 input1 '/' fileName]);
        featStr = featStr.features;
        featMat = struct2array(featStr);
        csvwrite([input0 'CSVfiles/' fileName '.txt'],featMat);
        disp([aa/length(list0) bb/length(list1)]);
    end
end
clear
