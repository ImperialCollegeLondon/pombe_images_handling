%% PLOTTING-ONLY AUXILIARY PART


% PLOTTING FEATURES
%colors8    = rand(8,3);
%colors8x2  = repelem(colors8,2,1); %unused from now on
%colors16   = [colors8;colors8*25/49]; colors16 = colors16(reshape(reshape(1:16,8,2)',1,16),:);
%colors16   = [colors8; colors8*5/7]; colors16 = colors16(reshape(reshape(1:16,8,2)',1,16),:);
colors8     = [[0.8 0.8 0.8]; [1 1 0]; [1 0 1]; [0 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [1 0.5 0.5]];
colors16    = [colors8;colors8*5/7]; colors16 = colors16(reshape(reshape(1:16,8,2)',1,16),:);
colorsDye   = colors16(repelem((1:16)',indSize),:);
marksDye    = {'^'}; marksDye16= repmat(marksDye,16,1);
marksDye    = repmat(marksDye,length(colorsDye),1);


% DATA PROCESSING
tExpN = length(namDye);             meanAllDye	= zeros(length(featNamDye),3*tExpN);
qualFeatsDye    = zeros(tExpN,6);	ratioMatDye = zeros(5,tExpN);
meanMatDye      = zeros(length(featNamDye),tExpN);	stdMatDye   = meanMatDye;
sizeDye         = zeros(tExpN,1);   stdAllDye   = zeros(length(featNamDye),3*tExpN);

for aa=1:tExpN
    tVar1   = cell2mat(featDye(aa));
    sizeDye(aa) = length(tVar1);
    tFlagA  = tVar1(:,4)==1;
    tFlagB  = tVar1(:,4)==2 & tVar1(:,6)==0;
    tFlagC  = tVar1(:,4)==2 & not(tVar1(:,6)==0);

    %calculating % of cell subpopulations
    qualFeatsDye(aa,1) = sum(tFlagA);
    qualFeatsDye(aa,3) = sum(tFlagA)+sum(tFlagB);
    qualFeatsDye(aa,:) = qualFeatsDye(aa,:)/size(tVar1,1);

    %calculating simple statistics
    meanMatDye(1:14,aa) = nanmean(tVar1(:,1:14),1);
    stdMatDye(1:14,aa)  = nanstd(tVar1(:,1:14),1);
    meanMatDye(5,aa)    = nanmean(tVar1(tFlagC,5));     %correction: nucleiD should be calculated only for septated
    stdMatDye(5,aa)     = nanstd(tVar1(tFlagC,5));      %correction: nucleiD should be calculated only for septated
    meanMatDye(end,aa)  = nanmean(tVar1(tFlagC,1));     %correction: length at end on septated
    stdMatDye(end,aa)   = nanstd(tVar1(tFlagC,1));      %correction: length at end on septated
    
    %calculating dye ratio stats
    ratioMatDye(1,aa)   = nanmean(tVar1(:,9));
    ratioMatDye(3,aa)   = sum(tVar1(:,9)==0)/size(tVar1,1);
    ratioMatDye(4,aa)   = nanstd(tVar1(:,13));
    %ratioMatDye(4,aa)  = nanstd([tVar2;-tVar2]);
    ratioMatDye(5,aa)   = nanstd(tVar1(not(tVar1(:,9)==0),13));
    %ratioMatDye(5,aa)  = nanstd([tVar2(not(tVar1(:,9)==0));-tVar2(not(tVar1(:,9)==0))]);
    
    %new matrix for dividing phases
    meanAllDye(1:14,aa)         = nanmean(tVar1(tFlagA,1:14),1);
    meanAllDye(1:14,aa+1*tExpN) = nanmean(tVar1(tFlagB,1:14),1);
    meanAllDye(1:14,aa+2*tExpN) = nanmean(tVar1(tFlagC,1:14),1);
    stdAllDye(1:14,aa)          = nanstd(tVar1(tFlagA,1:14),1);
    stdAllDye(1:14,aa+1*tExpN)  = nanstd(tVar1(tFlagB,1:14),1);
    stdAllDye(1:14,aa+2*tExpN)  = nanstd(tVar1(tFlagC,1:14),1);    
end
ratioMatDye(2,:) = ratioMatDye(1,:)./(1-ratioMatDye(3,:));
ratioNamDye =  {'Ratio tot' 'Ratio stained' '% unstained' 'Delta std tot' 'Delta std stained'};
meanAllDye(end,1:tExpN)          = 1;
meanAllDye(end,(tExpN+1):2*tExpN)= 3;
meanAllDye(end,(2*tExpN+1):end)  = 4;
clear tFlagA tFlagB tFlagC

qualFeatsDye(:,[2,4]) = 1 - qualFeatsDye(:,[1,3]);
qualFeatsDye(:,6) = qualFeatsDye(:,4)./qualFeatsDye(:,2);
qualFeatsDye(:,5) = 1 - qualFeatsDye(:,6);

propMatDye  = qualFeatsDye(:,[2,4,6])';
propNamDye  = [{'Binucleate %'} {'Septated %'} {'Sept on Binucleate %'}];

tR = cumsum(indSize);   tL = [1, 1+tR(1:end-1)];
meanMatDye16old = zeros(15,16);     stdMatDye16old = meanMatDye16old;
propMatDye16old = zeros(3,16);      sizeDye16old = zeros(16,1);
meanAllDye16old = zeros(15,48);     stdAllDye16old = zeros(15,48);
for aa = 1:16
    meanMatDye16old(:,aa) = nanmean(meanMatDye(:,tL(aa):tR(aa)),2);
    propMatDye16old(:,aa) = nanmean(propMatDye(:,tL(aa):tR(aa)),2);
    stdMatDye16old(:,aa)  = nanmean(stdMatDye(:,tL(aa):tR(aa)),2);
    sizeDye16old(aa)      = nanmean(sizeDye(tL(aa):tR(aa)));
    meanAllDye16old(:,aa)    = nanmean(meanAllDye(:,tL(aa):tR(aa)),2);
    meanAllDye16old(:,16+aa) = nanmean(meanAllDye(:,tExpN+(tL(aa):tR(aa))),2);
    meanAllDye16old(:,32+aa) = nanmean(meanAllDye(:,2*tExpN+(tL(aa):tR(aa))),2);
    stdAllDye16old(:,aa)     = nanmean(stdAllDye(:,tL(aa):tR(aa)),2);
    stdAllDye16old(:,16+aa)  = nanmean(stdAllDye(:,183+(tL(aa):tR(aa))),2);
    stdAllDye16old(:,32+aa)  = nanmean(stdAllDye(:,2*183+(tL(aa):tR(aa))),2);
end



% OTHER RATES CALCULATION

%elaborating calculated rates
grVec8 = mean(reshape(grVec16,2,8),1);
grVecDye  = repelem(grVec16,indSize);
grFlagDye = repelem(grFlag16,indSize);
[~,sort_grInd16]= sort(grVec16);
[~,sort_grInd8] = sort(grVec8);
[~,sort_grIndDye]  = sort(grVecDye);

%extracting data from Carlson paper
caVec8  = log(2)*[0.08;0.22;0.40;0.13;0.27;0.15;0.11;0.07];
caVec16 = repelem(caVec8,2,1);
caVecDye= repelem(caVec16,indSize);
[~,sort_caInd8]  = sort(caVec8);
[~,sort_caInd16] = sort(caVec16);
[~,sort_caIndDye]= sort(caVecDye);

%flagging wrong files/points
qFlagDye = zeros(1,length(featDye));
for aa=1:length(featDye);   qFlagDye(aa) = size(cell2mat(featDye(aa)),1);   end
qFlagDye = (qFlagDye<10000)';
wrongDye = or(grFlagDye,qFlagDye);
grFlagDye16  = false(16,1);  y = cumsum(indSize); x = [1 y(1:end-1)+1];    % flagging exps with all insufficient TPs
for ii=1:16;    grFlagDye16(ii) = sum(~wrongDye(x(ii):y(ii)));   end
grFlagDye16  = (grFlagDye16==0);


%manual flagging (overrides automatic flagging)
%all datasets: pheB,serA,serB;
%%grFlag16(:) = 0;  grFlag16([8,11,12]) = 1;
%%grFlagDye   = repelem(grFlag16,indSize);
%%grFlagDye16 = grFlag16; wrongDye = grFlagDye;
%tp specifici:	gly_A_t2, gly_B_t2, leu_B_t2, phe_A_t2, pro_A_t1, pro_A_t8, trp_A_t5
%%wrongDye([2,15,38,71,94,101,168]) = 1;
%con altri tp trovati analizzando gli outliers (funzione isoutlier, niente di campato in aria!
%%wrongDye([2,11,15,36,38,49,60,71,75,94,101,141,150,151,168]) = 1;




%% FINAL BIT
%transposing vectors
condsDye = condsDye';
mSizeDye = mSizeDye';

mainIndDye  = indHA2;


%new mean matrices: must be done after setting wrong timepoints/experiments
tExpN = 16;
qualFeatsDye16  = zeros(tExpN,6);	meanAllDye16	= zeros(15,3*tExpN);
meanMatDye16    = zeros(15,tExpN);	stdMatDye16     = meanMatDye16;

tS = 0;
for aa=1:tExpN
    tE   = tS+indSize(aa);    tInt = (tS+1):tE; 
    if not(all(wrongDye(tInt))); tInt(wrongDye(tInt))=[]; end
    tS = tE;
    tVar1   = vertcat(featDye{tInt});
    tFlagA  = tVar1(:,4)==1;
    tFlagB  = tVar1(:,4)==2 & tVar1(:,6)==0;
    tFlagC  = tVar1(:,4)==2 & not(tVar1(:,6)==0);

    %calculating % of cell subpopulations
    qualFeatsDye16(aa,1)     = sum(tFlagA);
    qualFeatsDye16(aa,3)     = sum(tFlagA)+sum(tFlagB);
    qualFeatsDye16(aa,:)     = qualFeatsDye16(aa,:)/size(tVar1,1);

    %calculating simple statistics
    meanMatDye16(1:14,aa)    = nanmean(tVar1(:,1:14),1);
    stdMatDye16(1:14,aa)     = nanstd(tVar1(:,1:14),1);
    meanMatDye16(5,aa)       = mean(tVar1(tFlagC,5));     %correction: nucleiD should be calculated only for septated
    stdMatDye16(5,aa)        = std(tVar1(tFlagC,5));      %correction: nucleiD should be calculated only for septated
    meanMatDye16(end,aa)     = mean(tVar1(tFlagC,1));     %correction: length at end on septated
    stdMatDye16(end,aa)      = std(tVar1(tFlagC,1));      %correction: length at end on septated
    
    %new matrix for dividing phases
    meanAllDye16(1:14,aa)         = nanmean(tVar1(tFlagA,1:14),1);
    meanAllDye16(1:14,aa+1*tExpN) = nanmean(tVar1(tFlagB,1:14),1);
    meanAllDye16(1:14,aa+2*tExpN) = nanmean(tVar1(tFlagC,1:14),1);
end
meanAllDye16(end,1:tExpN)          = 1;
meanAllDye16(end,(tExpN+1):2*tExpN)= 3;
meanAllDye16(end,(2*tExpN+1):end)  = 4;
clear tFlagA tFlagB tFlagC

qualFeatsDye16(:,[2,4]) = 1 - qualFeatsDye16(:,[1,3]);
qualFeatsDye16(:,6) = qualFeatsDye16(:,4)./qualFeatsDye16(:,2);
qualFeatsDye16(:,5) = 1 - qualFeatsDye16(:,6);
propMatDye16  = qualFeatsDye16(:,[2,4,6])';


clear aa indHA2 qFlagDye tVar1  tL tR
