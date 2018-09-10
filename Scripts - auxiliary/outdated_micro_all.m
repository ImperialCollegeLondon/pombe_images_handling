%% INPUTS
%reading inputs
input1  = '..\Files\treated-data\Microscopy\';
load([input1 'globStructure'],'names','features','counters');
namMic = names;
featMic = features;
featNamMic = fieldnames(featMic);
for aa=1:length(featNamMic);    featNamMic(aa)= strrep(featNamMic(aa),'_',' '); end
featMic = struct2cell(featMic);
for aa=1:length(namMic)   
    namMic(aa)= strrep(namMic(aa),'_',' ');
    namMic(aa)= strrep(namMic(aa),'unfixed','unf M');
    namMic(aa)= strrep(namMic(aa),'fixed','fix M');    
    namMic(aa)= strrep(namMic(aa),'nh4cl','NH4');
    namMic(aa)= strrep(namMic(aa),'r3','C');
    namMic(aa)= strrep(namMic(aa),'r2','B');
    namMic(aa)= strrep(namMic(aa),'r1','A');
    namMic(aa)= strrep(namMic(aa),'.mat','');
    tVar1 = cell2mat(namMic(aa));     tVar1(7)  = upper(tVar1(7));
    namMic(aa)= {tVar1};
end
namMic = namMic';


%removing wrongly processed images
for aa=1:size(featMic,3)
    tVar1 = cell2mat(featMic(4,:,aa))==0;
    for bb=1:size(featMic,1)
        tVar2 = cell2mat(featMic(bb,:,aa));
        tVar2(tVar1)=[];
        featMic(bb,:,aa) = {tVar2};
    end
end



%% GROWTH RATES
%extracting data from turbidostat
run('..\Scripts - auxiliary\aux_noD_rateCalc.m')
grVecMic = grVecNoD; grFlagMic = grFlagNoD; clear grVecNoD grFlag NoD;
[~,sort_grIndMic] = sort(grVecMic);


%extracting data from Carlson paper
caVecMic = log(2)*[0.08, 0.22, 0.40, 0.13, 0.27, 0.15, 0.11, 0.07];
caVecMic = repelem(caVecMic,repelem(3,8))';
caVecMic = [caVecMic;caVecMic];
[~,sort_caIndMic] = sort(caVecMic);


%doubling appropriate vectors (for fixed and unfixed cells)
grFlagMic   = [grFlagMic;grFlagMic];
grVecMic    = [grVecMic;grVecMic];
sort_grIndMic = [sort_grIndMic;sort_grIndMic];



%% FINAL STATS
%calculating % of cell subpopulations
for aa = 1:size(featMic,3)
    tVar1 = cell2mat(featMic([4,6,7],1,aa));
    tVar1 = reshape(tVar1,length(tVar1)/3,3);
    tVar1(tVar1(:,1)==0,:) = [];
      
    qualFeatsMic(aa,1) = sum(tVar1(:,1)==1);
    qualFeatsMic(aa,3) = sum(tVar1(:,3)==0);
    qualFeatsMic(aa,5) = sum(tVar1(:,2)==0);
    qualFeatsMic(aa,:) = qualFeatsMic(aa,:)/size(tVar1,1);
end
qualFeatsMic(:,[2,4,6]) = 1 - qualFeatsMic(:,[1,3,5]);
qualFeatsMic(:,8) = qualFeatsMic(:,6)./qualFeatsMic(:,2);
qualFeatsMic(:,7) = 1 - qualFeatsMic(:,8);
%qualFeatsMic = [mono;bin;nMit;yMit;nSep;ySep;nSeB;ySeB];
propMatMic  = qualFeatsMic(:,[2,6,4])';
propMatMic2 = qualFeatsMic(:,[2,6,8])';
propNamMic  = [{'Binucleate %'} {'Septated %'} {'Mitotic %'}];
propNamMic2 = [{'Binucleate %'} {'Septated %'} {'Sept on Bin %'}];


%calculating simple statistics
meanMatMic = zeros(7,48);
stdMatMic  = zeros(7,48);
for aa=1:size(featMic,1)
    for bb=1:size(featMic,3)
        meanMatMic(aa,bb) = mean(cell2mat(featMic(aa,:,bb)));
        stdMatMic(aa,bb)  = std(cell2mat(featMic(aa,:,bb)));
    end
end


%final bits
Mconds8  = [{'M Gly'} {'M Leu'} {'M NH4Cl'} {'M Phe'} {'M Pro'} {'M Ser'} {'M Thr'} {'M Trp'}]';
for aa=1:8
    Mconds16(aa) = {[cell2mat(Mconds8(aa)) ' fix']};
    Mconds16(aa+8)={[cell2mat(Mconds8(aa)) ' unf']};
end
%colors8    = [[1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [0.5 0.5 0.5]; [0 0 0]];
%colors88   = [colors8;5/7*colors8];
colors8     = [[0.8 0.8 0.8]; [1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [1 0.5 0.5]];
colorsMic   = [colors8;colors8*5/7;colors8*25/49];  
colorsMic   = colorsMic(reshape(reshape(1:24,8,3)',1,24),:);  colorsMic   = [colorsMic; colorsMic];
marksMic    = {'d'}; marksMic= repmat(marksMic,length(colorsMic),1);


%flagging wrong files/points
qFlagMic = zeros(1,length(featMic));
for aa=1:length(featMic)
    qFlagMic(aa) = size(cell2mat(featMic(1,1,aa)),1);
end
qFlagMic = (qFlagMic<300)';
wrongMic = or(grFlagMic,qFlagMic);

clear aa bb counters features input1 list1 names tVar1 tVar2 turbo