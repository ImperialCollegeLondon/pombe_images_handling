%% INPUTS
useUNFIX = false;
%unfixed files are not considered!
rootPath0 = '..\..\..\Cycle sensor data\2. Dati sperimentali\2. Dati finali\';
list0 = dir(rootPath0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
%for aa=1:length(list0)
aa=4; %manual override:  leafPath0 will be 'MAT_MICRO\'
    leafPath0 = list0(aa).name;
    rootPath1 = [rootPath0 leafPath0 '\'];
    list1   = dir(rootPath1); list1(1:2) = [];
    namMic  = cell(length(list1),1);  conds24 = namMic;
    featMic = cell(9,1,length(list1));
    phaseMic = cell(length(list1),1);
    %orSizeMic= zeros(length(list1),1);
    sizeMic= zeros(length(list1),1);
    for bb=1:length(list1)
        leafPath1   = list1(bb).name;
        rootPath2   = [rootPath1 leafPath1];
        namMic(bb)  = {[leafPath1(1) leafPath1(5) ' ' strrep(leafPath1(end-8:end-4),'_',' ')]};
        conds24(bb) = {strrep(leafPath1(end-8:end-4),'_',' ')};
        
        tTab = readtable(rootPath2);    featNamMic = fieldnames(tTab);
        %orSizeMic(bb)	= size(tTab,1);
        goodInd         = not(table2array(tTab(:,15)) ==0);
        phaseMic(bb)    = {table2array(tTab(goodInd,18))};
        tTab = tTab(goodInd,[2,4,1,15,14,16,17,12,13,3]);
        sizeMic(bb)     = size(tTab,1);
        for cc=1:9
            featMic(cc,1,bb) = {table2array(tTab(:,cc))};
        end
    end
    featNamMic = featNamMic([2,4,1,15,14,16,17,12,13,3]);
    featNamMic = strrep(featNamMic,'_',' ');    featNamMic = strrep(featNamMic,'perNuc', 'per nucleus');
    featNamMic = [featNamMic;{'cell length_end'}];  %new entry 15/7/2018
%end



%% GROWTH RATES
%extracting data from turbidostat
run('..\Scripts - auxiliary\aux_noD_rateCalc.m')
grVecMic = grVecNoD; grFlagMic = grFlagNoD; clear grVecNoD grFlagNoD;
[~,sort_grIndMic] = sort(grVecMic);

%extracting data from Carlson paper
caVecMic = log(2)*[0.08, 0.22, 0.40, 0.13, 0.27, 0.15, 0.11, 0.07];
caVecMic = repelem(caVecMic,repelem(3,8))';
[~,sort_caIndMic] = sort(caVecMic);



%% FINAL STATS
%calculating % of cell subpopulations
qualFeatsMic = zeros(size(featMic,3),8);
for aa = 1:length(phaseMic)
    tVar1 = phaseMic{aa};      
    qualFeatsMic(aa,1) = sum(tVar1==1 | tVar1==2);
    qualFeatsMic(aa,4) = sum(tVar1==2 | tVar1==3);
    qualFeatsMic(aa,6) = sum(tVar1==4);
    qualFeatsMic(aa,:) = qualFeatsMic(aa,:)/length(tVar1);
end
qualFeatsMic(:,[2,3,5]) = 1 - qualFeatsMic(:,[1,4,6]);
qualFeatsMic(:,8) = qualFeatsMic(:,6)./qualFeatsMic(:,2);
qualFeatsMic(:,7) = 1 - qualFeatsMic(:,8);
%qualFeatsMic = [mono;bin;nMit;yMit;nSep;ySep;nSeB;ySeB];
propMatMic  = qualFeatsMic(:,[2,6,4])';
propMatMic2 = qualFeatsMic(:,[2,6,8])';
propNamMic  = [{'Binucleate %'} {'Septated %'} {'Mitotic %'}];
propNamMic2 = [{'Binucleate %'} {'Septated %'} {'Sept on Bin %'}];


%calculating simple statistics
tExpN = size(featMic,3);    meanAllMic = zeros(length(featNamMic),3*tExpN);  
meanMatMic = zeros(length(featNamMic),tExpN);    stdMatMic = meanMatMic;
  
for aa=1:tExpN
    tPhaseMic   = phaseMic{aa};
    tVar0       = featMic(:,:,aa);
    for bb=1:size(featMic,1)
        tVar1 = tVar0{bb};
        meanMatMic(bb,aa) = mean(tVar1);
        stdMatMic(bb,aa)  = std(tVar1);
        
        meanAllMic(bb,aa)        = mean(tVar1(tPhaseMic <= 2));
        meanAllMic(bb,aa+1*tExpN)= mean(tVar1(tPhaseMic == 3));
        meanAllMic(bb,aa+2*tExpN)= mean(tVar1(tPhaseMic == 4));        
    end
    
    tFlag = tVar0{4}==2 & not(tVar0{6}==0);
    meanMatMic(5,aa)    = mean(tVar0{5}(tFlag));	%correction: nucleiD should be calculated only for septated
    stdMatMic(5,aa)     = std(tVar0{5}(tFlag));     %correction: nucleiD should be calculated only for septated
    meanMatMic(end,aa)	= mean(tVar0{1}(tFlag));	%correction: length at end on septated
    stdMatMic(end,aa)	= std(tVar0{1}(tFlag));     %correction: length at end on septated    
end
meanAllMic(end,1:tExpN) = 1;
meanAllMic(end,(tExpN+1):2*tExpN) = 3;
meanAllMic(end,(2*tExpN+1):end) = 4;
clear tFlag


%final bits
Mconds8  = [{'MM Gly'} {'MM Leu'} {'MM NH4Cl'} {'MM Phe'} {'MM Pro'} {'MM Ser'} {'MM Thr'} {'MM Trp'}]';
%colors8    = [[1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [0.5 0.5 0.5]; [0 0 0]];
%colors88   = [colors8;5/7*colors8];
colors8     = [[0.8 0.8 0.8]; [1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [1 0.5 0.5]];
colorsMic   = [colors8;colors8*5/7;colors8*25/49];  
colorsMic   = colorsMic(reshape(reshape(1:24,8,3)',1,24),:);
marksMic    = {'o'}; marksMic= repmat(marksMic,length(colorsMic),1);



%% DOUBLING VECTORS IF UNFIXED ARE CONSIDERED, AND CLOSING
if useUNFIX
    %doubling appropriate vectors (for fixed and unfixed cells)
    caVecMic    = [caVecMic;caVecMic];
    grFlagMic   = [grFlagMic;grFlagMic];
    grVecMic    = [grVecMic;grVecMic];
    sort_grIndMic = [sort_grIndMic;sort_grIndMic];
    sort_caIndMic = [sort_caIndMic;sort_caIndMic];
    colorsMic   = [colorsMic; colorsMic];
    marksMic    = [marksMic; marksMic];
    for aa=1:8
        Mconds16(aa) = {[cell2mat(Mconds8(aa)) ' fix']};
        Mconds16(aa+8)={[cell2mat(Mconds8(aa)) ' unf']};
    end
end

%flagging wrong files/points
qFlagMic = zeros(1,length(featMic));
for aa=1:length(featMic)
    qFlagMic(aa) = size(cell2mat(featMic(1,1,aa)),1);
end
qFlagMic = (qFlagMic<300)';
wrongMic = or(grFlagMic,qFlagMic);

%manual flagging (overrides automatic flagging)
%glyC, pheC, serA,serB,serC, thrC, trpC
grFlagMic(:) = 0; grFlagMic([3,12,16,17,18,21,24])=1;
wrongMic = grFlagMic;


clear aa bb cc goodInd input1 leafPath0 leafPath1 list0 list1 names qFlagMic rootPath0 rootPath1 rootPath2 tTab tVar1 tVar2 useUNFIX turbo