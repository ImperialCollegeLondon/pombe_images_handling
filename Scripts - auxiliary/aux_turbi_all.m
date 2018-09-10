%reading inputs
rootPath0 = '..\..\..\Cycle sensor data\2. Dati sperimentali\2. Dati finali\';
list0 = dir(rootPath0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
%for aa=1:length(list0)
aa=5; %manual override:  leafPath0 will be 'MAT_TURBI\'
    leafPath0   = list0(aa).name;
    rootPath1   = [rootPath0 leafPath0 '\'];
    list1       = dir(rootPath1); list1(1:2) = [];
    namTur      = cell(length(list1),1);
    featTur     = cell(10,1,length(list1));
    phaseTur    = cell(length(list1),1);
    %orSizeTur  = zeros(length(list1),1);
    sizeTur     = zeros(length(list1),1);
    for bb=1:length(list1)
        leafPath1   = list1(bb).name;
        rootPath2   = [rootPath1 leafPath1];
        namTur(bb)  = {[leafPath1(1) leafPath1(5) ' ' strrep(leafPath1(end-8:end-4),'_',' ')]};
   
        tTab = readtable(rootPath2);    featNamTur = fieldnames(tTab);
        %orSizeTur(bb)   = size(tTab,1);
        goodInd         = not(table2array(tTab(:,15)) ==0);
        phaseTur(bb)	= {table2array(tTab(goodInd,18))};
        tTab = tTab(goodInd,[2,4,1,15,14,16,17,12,13,3]);
        sizeTur(bb)     = size(tTab,1);
        for cc=1:10
            featTur(cc,1,bb) = {table2array(tTab(:,cc))};
        end
    end
    featNamTur = featNamTur([2,4,1,15,14,16,17,12,13,3]);
    featNamTur = strrep(featNamTur,'_',' ');    featNamTur = strrep(featNamTur,'perNuc', 'per nucleus');
    featNamTur = [featNamTur;{'cell length_end'}];  %new entry 15/7/2018
%end


%extracting data from turbidostat
run('..\Scripts - auxiliary\aux_noD_rateCalc.m')
grVecTur = grVecNoD; grFlagTur = grFlagNoD; clear grVecNoD grFlagNoD;
%grVecTur(17) = []; grFlagTur(17) = [];     %TEMP: removing Ser_r2 from analysis
grVecTur(18) = []; grFlagTur(18) = [];      %NEW, after swapping: removing Ser_r2 from analysis
[~,sort_grIndTur] = sort(grVecTur);


%extracting data from Carlson paper
caVecTur = log(2)*[0.08, 0.22, 0.40, 0.13, 0.27, 0.15, 0.11, 0.07];
caVecTur = repelem(caVecTur,repelem(3,8))'; caVecTur(17)=[];
[~,sort_caIndTur] = sort(caVecTur);


%calculating % of cell subpopulations
qualFeatsTur = zeros(size(featTur,3),8);
for aa = 1:length(phaseTur)
    tVar1 = phaseTur{aa};      
    qualFeatsTur(aa,1) = sum(tVar1==1 | tVar1==2);
    qualFeatsTur(aa,4) = sum(tVar1==2 | tVar1==3);
    qualFeatsTur(aa,6) = sum(tVar1==4);
    %qualFeatsTur(aa,:) = qualFeatsTur(aa,:)/length(tVar1);
    qualFeatsTur(aa,:) = qualFeatsTur(aa,:)/sum(not(isnan(tVar1)));
end
qualFeatsTur(:,[2,3,5]) = 1 - qualFeatsTur(:,[1,4,6]);
qualFeatsTur(:,8) = qualFeatsTur(:,6)./qualFeatsTur(:,2);
qualFeatsTur(:,7) = 1 - qualFeatsTur(:,8);
%qualFeatsTur = [mono;bin;nMit;yMit;nSep;ySep;nSeB;ySeB];
propMatTur  = qualFeatsTur(:,[2,6,4])';
propMatTur2 = qualFeatsTur(:,[2,6,8])';
propNamTur  = [{'Binucleate %'} {'Septated %'} {'Mitotic %'}];
propNamTur2 = [{'Binucleate %'} {'Septated %'} {'Sept on Bin %'}];


%calculating simple statistics
tExpN = size(featTur,3);     
meanMatTur = zeros(length(featNamTur),tExpN);    stdMatTur = meanMatTur;
meanAllTur = zeros(length(featNamTur),3*tExpN); stdAllTur = meanAllTur;
for aa=1:tExpN
    tPhaseTur   = phaseTur{aa};
    tVar0       = featTur(:,:,aa);
    for bb=1:size(featTur,1)
        tVar1 = tVar0{bb};
        meanMatTur(bb,aa) = nanmean(tVar1);
        stdMatTur(bb,aa)  = nanstd(tVar1);
        
        meanAllTur(bb,aa)        = nanmean(tVar1(tPhaseTur <= 2));
        meanAllTur(bb,aa+1*tExpN)= nanmean(tVar1(tPhaseTur == 3));
        meanAllTur(bb,aa+2*tExpN)= nanmean(tVar1(tPhaseTur == 4));
        stdAllTur(bb,aa)         = nanstd(tVar1(tPhaseTur <= 2));
        stdAllTur(bb,aa+1*tExpN) = nanstd(tVar1(tPhaseTur == 3));
        stdAllTur(bb,aa+2*tExpN) = nanstd(tVar1(tPhaseTur == 4));        
    end
    
    tFlag = tVar0{4}==2 & not(tVar0{6}==0);
    %meanMatTur(5,aa)    = nanmean(tVar0{5}(tFlag));     %correction: nucleiD should be calculated only for septated
    meanMatTur(5,aa)    = nanmean(tVar0{5}(tFlag & tVar0{5}>4));     %new: temporary, this change should be done on final files
    %stdMatTur(5,aa)     = nanstd(tVar0{5}(tFlag));      %correction: nucleiD should be calculated only for septated
    stdMatTur(5,aa)     = nanstd(tVar0{5}(tFlag & tVar0{5}>4));     %new: temporary, this change should be done on final files
    meanMatTur(end,aa)	= nanmean(tVar0{1}(tFlag));     %correction: length at end on septated
    stdMatTur(end,aa)	= nanstd(tVar0{1}(tFlag));      %correction: length at end on septated    
end
meanAllTur(end,1:tExpN) = 1;
meanAllTur(end,(tExpN+1):2*tExpN) = 3;
meanAllTur(end,(2*tExpN+1):end) = 4;
clear tFlag

    
%final bits
Tconds8  = [{'MT Gly'} {'MT Leu'} {'MT NH4Cl'} {'MT Phe'} {'MT Pro'} {'MT Ser'} {'MT Thr'} {'MT Trp'}]';
%colors8    = [[1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [0.5 0.5 0.5]; [0 0 0]];
%colorsTur  = repelem(colors8,3,1); colorsTur(17,:)=[]; 
colors8     = [[0.8 0.8 0.8]; [1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [1 0.5 0.5]];
colorsTur   = [colors8;colors8*5/7;colors8*25/49];  
colorsTur   = colorsTur(reshape(reshape(1:24,8,3)',1,24),:);  colorsTur(17,:)=[];
marksTur    = {'s'}; marksTur= repmat(marksTur,length(colorsTur),1);

%flagging wrong files/points
qFlagTur = zeros(1,length(featTur));
for aa=1:length(featTur)
    qFlagTur(aa) = size(cell2mat(featTur(1,1,aa)),1);
end
qFlagTur = (qFlagTur<5000)';
wrongTur = or(grFlagTur,qFlagTur);

%manual flagging (overrides automatic flagging)
%glyC, pheC, serA,serB,(serC), thrC, trpC
%grFlagTur(:) = 0; grFlagTur([3,12,16,17,20,23])=1;
wrongTur = grFlagTur;


clear aa bb cc goodInd input1 leafPath0 leafPath1 list0 list1 names qFlagTur rootPath0 rootPath1 rootPath2 tTab tVar1 tVar2 turbo
