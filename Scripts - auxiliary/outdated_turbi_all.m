%reading inputs
input1  = '..\Files\treated-data\Growth\';
load([input1 'globStructure'],'names','features','counters');
namTur = names;
featTur = features;
featNamTur = fieldnames(featTur);
for aa=1:length(featNamTur);    featNamTur(aa)= strrep(featNamTur(aa),'_',' '); end
featNamTur(1) = strrep(featNamTur(1),'ht','th');       %correcting "length" misspelling
featTur = struct2cell(featTur);
for aa=1:length(namTur)
    namTur(aa)= strrep(namTur(aa),'_',' ');
    namTur(aa)= strrep(namTur(aa),'157 ','T ');
    namTur(aa)= strrep(namTur(aa),' 50','');
    namTur(aa)= strrep(namTur(aa),'NH4Cl','NH4');
    namTur(aa)= strrep(namTur(aa),' d2T','');
    namTur(aa)= strrep(namTur(aa),'r3','C');
    namTur(aa)= strrep(namTur(aa),'r2','B');
    namTur(aa)= strrep(namTur(aa),'r1','A');
    namTur(aa)= strrep(namTur(aa),'fixed','fix');
    namTur(aa)= strrep(namTur(aa),'fix','');  %per semplicita'
    namTur(aa)= strrep(namTur(aa),'  ',' ');
end
namTur = namTur';


%removing wrongly processed images
for aa=1:size(featTur,3)
    tVar1 = cell2mat(featTur(4,:,aa))==0;
    for bb=1:size(featTur,1)
        tVar2 = cell2mat(featTur(bb,:,aa));
        tVar2(tVar1)=[];
        featTur(bb,:,aa) = {tVar2};
    end
end


%extracting data from turbidostat
run('..\Scripts - auxiliary\aux_noD_rateCalc.m')
grVecTur = grVecNoD; grFlagTur = grFlagNoD; clear grVecNoD grFlag NoD;
%grVecTur(17) = []; grFlagTur(17) = [];     %TEMP: removing Ser_r2 from analysis
grVecTur(18) = []; grFlagTur(18) = [];      %NEW, after swapping: removing Ser_r2 from analysis
[~,sort_grIndTur] = sort(grVecTur);


%extracting data from Carlson paper
caVecTur = log(2)*[0.08, 0.22, 0.40, 0.13, 0.27, 0.15, 0.11, 0.07];
caVecTur = repelem(caVecTur,repelem(3,8))'; caVecTur(17)=[];
[~,sort_caIndTur] = sort(caVecTur);


%calculating % of cell subpopulations
for aa = 1:size(featTur,3)
    tVar1 = cell2mat(featTur([4,6,7],1,aa));
    tVar1 = reshape(tVar1,length(tVar1)/3,3);
    tVar1(tVar1(:,1)==0,:) = [];
      
    qualFeatsTur(aa,1) = sum(tVar1(:,1)==1);
    qualFeatsTur(aa,3) = sum(tVar1(:,3)==0);
    qualFeatsTur(aa,5) = sum(tVar1(:,2)==0);
    qualFeatsTur(aa,:) = qualFeatsTur(aa,:)/size(tVar1,1);
end
qualFeatsTur(:,[2,4,6]) = 1 - qualFeatsTur(:,[1,3,5]);
qualFeatsTur(:,8) = qualFeatsTur(:,6)./qualFeatsTur(:,2);
qualFeatsTur(:,7) = 1 - qualFeatsTur(:,8);
%qualFeatsTur = [mono;binu;nMit;yMit;nSep;ySep;nSeB;ySeB];
propMatTur  = qualFeatsTur(:,[2,6,4])';
propMatTur2 = qualFeatsTur(:,[2,6,8])';
propNamTur  = [{'Binucleate %'} {'Septated %'} {'Mitotic %'}];
propNamTur2 = [{'Binucleate %'} {'Septated %'} {'Sept on Bin %'}];


%calculating simple statistics
meanMatTur = zeros(7,23);
stdMatTur  = zeros(7,23);
for aa=1:size(featTur,1)
    for bb=1:size(featTur,3)
        meanMatTur(aa,bb) = mean(cell2mat(featTur(aa,:,bb)));
        stdMatTur(aa,bb)  = std(cell2mat(featTur(aa,:,bb)));
    end
end


%final bits
Tconds8  = [{'T Gly'} {'T Leu'} {'T NH4Cl'} {'T Phe'} {'T Pro'} {'T Ser'} {'T Thr'} {'T Trp'}]';
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
qFlagTur = (qFlagTur<10000)';
wrongTur = or(grFlagTur,qFlagTur);

clear aa bb counters features input1 list1 names tVar1 tVar2 turbo