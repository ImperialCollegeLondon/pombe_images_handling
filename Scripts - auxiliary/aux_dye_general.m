run('..\Scripts - auxiliary\aux_dye_rateCalc.m')    %extracting data from OD measurements
rootPath0 = '..\..\..\Cycle sensor data\2. Dati sperimentali\2. Dati finali\';
list0 = dir(rootPath0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
%for aa=1:length(list0)
aa=3; %manual override:  leafPath0 will be 'MAT_DYE\'
    leafPath0   = list0(aa).name;
    rootPath1   = [rootPath0 leafPath0 '\'];
    list1       = dir(rootPath1); list1(1:2) = [];
    Dconds16    = cell(length(list1),1); Dconds8=Dconds16;  conds8=Dconds16;
    featDye16   = cell(12,1,length(list1));
    indSize     = zeros(1,length(list1));
    %orSizeDye  = zeros(length(list1),1);
    sizeDye16   = zeros(length(list1),1);
    
    usedLines=0; TPonly=[]; shTPs = [];
    for bb=1:length(list1)
        leafPath1 = list1(bb).name;
        rootPath2 = [rootPath1 leafPath1];
        Dconds16(bb) = {[leafPath1(1) leafPath1(5) ' ' strrep(leafPath1(end-8:end-4),'_',' ')]};
        Dconds8(bb)  = {[leafPath1(1) leafPath1(5) ' ' strrep(leafPath1(end-8:end-6),'_',' ')]};
        conds8(bb)   = {strrep(leafPath1(end-8:end-6),'_',' ')};
        tTab = readtable(rootPath2);    featNamDye = fieldnames(tTab);
        %orSizeDye(bb) = size(tTab,1);
        tTab = tTab(not(table2array(tTab(:,15)) ==0),:);
        sizeDye16(bb) = size(tTab,1);
        tTPs = tTab.tPoint;
        tTab = tTab(:,[2,4,1,15,14,16, 6,5,9,10, 12,13,8,3]);
        
        tX = cellfun('length',tTPs);    tTPs(tX==2) = strrep(tTPs(tX==2),'t','t0');  tTPs = cell2mat(tTPs);
        tTPonly = str2num(tTPs(:,2:end)); %#ok<ST2NM>
        tTPunic = unique(tTPonly);  indSize(bb) = length(tTPunic);
        for cc = 1:length(tTPunic)
            usedLines = usedLines+1;
            featDye(usedLines) = {table2array(tTab(ismember(tTPonly,tTPunic(cc)),:))};
            disp([usedLines bb cc]);
        end              
        for cc=1:14
            featDye16(cc,1,bb) = {table2array(tTab(:,cc))};
        end
        shTPs   = [shTPs;unique(tTPs,'rows')]; %#ok<AGROW>
        TPonly  = [TPonly;tTPunic]; %#ok<AGROW>
    end
    featNamDye = featNamDye([2,4,1,15,14,16, 6,5,9,10, 12,13,8,3]);
    featNamDye = strrep(featNamDye,'_',' ');    featNamDye = strrep(featNamDye,'perNuc', 'per nucleus');
    featNamDye = [featNamDye;{'cell length_end'}];      %new entry 15/7/2018
    Dconds8 = unique(Dconds8);
    conds8  = unique(conds8);
    tempSum = cumsum(indSize);
    indHA2  = [0 tempSum]+6; indHA2(end) = [];
    mSizeDye=[];  condsDye=[];
    for bb=1:length(list1)
        mSizeDye = [mSizeDye 1:indSize(bb)]; %#ok<AGROW>
        condsDye = [condsDye repelem(Dconds16(bb),1,indSize(bb))]; %#ok<AGROW>
    end
    namDye = cell(length(condsDye),1);  shRepl=namDye;
    for bb=1:length(condsDye)
        namDye(bb) = strcat(condsDye(bb),[' ',shTPs(bb,:)]);
        tX = cell2mat(namDye(bb));  shRepl(bb) = {tX(end-4:end)};
    end
%end

%Dconds8    = [{'MD Gly'} {'MD Leu'} {'MD NH4Cl'} {'MD Phe'} {'MD Pro'} {'MD Ser'} {'MD Thr'} {'MD Trp'}]';
clear aa bb cc leafPath0 leafPath1 list0 list1 rootPath0 rootPath1 rootPath2 tempSum tTPonly tTPs tTPunic tTab tX usedLines