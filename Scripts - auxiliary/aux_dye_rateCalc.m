%% PART 1: OD reading and rearrangment

%extracting data from OD measurements
GLY = table2array(readtable('..\Files\Dye ODs\GLY_ODs.txt'));
LEU = table2array(readtable('..\Files\Dye ODs\LEU_ODs.txt'));
NH4 = table2array(readtable('..\Files\Dye ODs\NH4_ODs.txt'));
PHE = table2array(readtable('..\Files\Dye ODs\PHE_ODs.txt'));
PRO = table2array(readtable('..\Files\Dye ODs\PRO_ODs.txt'));
SER = table2array(readtable('..\Files\Dye ODs\SER_ODs.txt'));
THR = table2array(readtable('..\Files\Dye ODs\THR_ODs.txt'));
TRP = table2array(readtable('..\Files\Dye ODs\TRP_ODs.txt'));

r01Tab =[{GLY(:,1)} {LEU(:,1)} {NH4(:,1)} {PHE(:,1)} {PRO(:,1)} {SER(:,1)} {THR(:,1)} {TRP(:,1)}];
r10Tab =[{GLY(:,2)} {LEU(:,2)} {NH4(:,2)} {PHE(:,2)} {PRO(:,2)} {SER(:,2)} {THR(:,2)} {TRP(:,2)}];
timTab =[{GLY(:,3)} {LEU(:,3)} {NH4(:,3)} {PHE(:,3)} {PRO(:,3)} {SER(:,3)} {THR(:,3)} {TRP(:,3)}];

OD_dye = [GLY(:,1);GLY(:,2); LEU(:,1);LEU(:,2); NH4(:,1);NH4(:,2); PHE(:,1);PHE(:,2); ...
          PRO(:,1);PRO(:,2); SER(:,1);SER(:,2); THR(:,1);THR(:,2); TRP(:,1);TRP(:,2)];
nonExistant = [26,116+12,116+13+12,179+(8:13)]; OD_dye(nonExistant) = [];
OD_DI= [GLY([1,4,7,10,12],2);    GLY([1,4,7,10,13],1);   LEU([1,4,7,10],2);      LEU([1,4,7,10],1); ...
        NH4([1,4,7,10],2);       NH4([1,4,7,10],1);      PHE([1,4,7,10,12],2);   PHE([1,4,6,7,10,12],1); ...
        PRO([1,4,7,10],2);       PRO([1,4,7,10],1);      SER([1,4,7,10,13],2);   SER([1,4,7,10,13],1); ...
        THR([1,4,7,10,12],2);    THR([1,4,7,10,12],1);   TRP([1,4,7],2);         TRP([1,4,7,10,13],1)];
clear GLY LEU NH4 PHE PRO SER THR TRP nonExistant



%% PART 2: growth rate calculation
tVar1 = [1,2,3,4,5,6,7];
aMat = zeros(16,length(tVar1)); rMat = zeros(16,length(tVar1));

for aa=1:8
    %lettura e rielaborazione input
    tMat01 = [cell2mat(timTab(aa)), cell2mat(r01Tab(aa))]';
    tMat10 = [cell2mat(timTab(aa)), cell2mat(r10Tab(aa))]';
      
    %working on r01
    for bb=1:length(tVar1)
        tVar2  = tVar1(bb):size(tMat01,2);
        aExp   = polyfit(tMat01(1,tVar2), log(tMat01(2,tVar2)),1); aExp=aExp(1);
        aMat(2*aa-1,bb) = aExp;
        yExp = tMat01(2,tVar1(bb))* exp(aExp*tMat01(1,:)- aExp*tMat01(1,tVar1(bb)));
%        rMat(2*aa-1,bb) = sqrt(sum((tMat01(2,:)-yExp).^2));
%        rMat(2*aa-1,bb) = sqrt(sum((tMat01(2,:)-yExp).^2))/mean(tMat01(2,:));
%        rMat(2*aa-1,bb) = sqrt(sum((tMat01(2,:)-yExp).^2))/sum((tMat01(2,:)-mean(tMat01(2,:))).^2);
        rMat(2*aa-1,bb) = sqrt(sum((tMat01(2,:)-yExp).^2))*mean(tMat01(2,:))/sum((tMat01(2,:)-mean(tMat01(2,:))).^2);
    end

    %working on r10
    for bb=1:length(tVar1)
        tVar2  = tVar1(bb):size(tMat10,2);
        aExp   = polyfit(tMat10(1,tVar2), log(tMat10(2,tVar2)),1); aExp=aExp(1);
        aMat(2*aa,bb) = aExp;
        yExp = tMat10(2,tVar1(bb))* exp(aExp*tMat10(1,:)- aExp*tMat10(1,tVar1(bb)));
%        rMat(2*aa,bb) = sqrt(sum((tMat10(2,:)-yExp).^2));
%        rMat(2*aa,bb) = sqrt(sum((tMat10(2,:)-yExp).^2))/mean(tMat10(2,:));
%        rMat(2*aa,bb) = sqrt(sum((tMat10(2,:)-yExp).^2))/sum((tMat10(2,:)-mean(tMat10(2,:))).^2);
        rMat(2*aa,bb) = sqrt(sum((tMat10(2,:)-yExp).^2))*mean(tMat10(2,:))/sum((tMat10(2,:)-mean(tMat10(2,:))).^2);
    end
end

[a,b] = min(rMat,[],2); grVec16= zeros(16,1);
for aa=1:16; grVec16(aa) = aMat(aa,b(aa));   end



%% PART 3: PLOTTING
%figure('units','normalized','outerposition',[0 0 1 1]);
%ha = tight_subplot(2,8,[.06 .025],[.15 .1],[.03 .03]);     %[in vert, in hor] [lower upper] [left right]
%for aa=1:8    
    %working on r01
%    axes(ha(aa*2-1)); %#ok<LAXES>
%    tMat01  = [cell2mat(timTab(aa)), cell2mat(r01Tab(aa))]';
%    bb      = b(2*aa-1);
%    aExp    = grVec16(2*aa-1);
%    yExp    = tMat01(2,bb)* exp(aExp*tMat01(1,:)- aExp*tMat01(1,bb));
%    plot(tMat01(1,:),tMat01(2,:),'Color','blue','Marker','*'); hold on
%    plot(tMat01(1,:),yExp,'Color','red'); hold on

    %working on r10
%    axes(ha(aa*2)); %#ok<LAXES>
%    tMat10  = [cell2mat(timTab(aa)), cell2mat(r10Tab(aa))]'; 
%    bb      = b(2*aa);
%    aExp    = grVec16(2*aa);
%    yExp    = tMat10(2,bb)* exp(aExp*tMat10(1,:)- aExp*tMat10(1,bb));
%    plot(tMat10(1,:),tMat10(2,:),'Color','blue','Marker','*'); hold on
%    plot(tMat10(1,:),yExp,'Color','red'); hold on    
%end



%% FINAL BITS
grVec16  = grVec16/log(2);   %TEMPORARILY CORRECTION!
caVec16  = repelem(log(2)*[0.08;0.22;0.40;0.13;0.27;0.15;0.11;0.07],2,1);
%grFlag16 = zeros(length(grVec16),1);           %all conditions are accepted
%grFlag16 = (a>0.50);                           %removes conditions where grVec calculation was not precise --> don't use: the calculation is good, the OD measurements are shitty! 
grFlag16 = abs(caVec16-grVec16)./caVec16 >0.7;  %removes condition markedly slower than the expected growth rate --> maybe a bit too stringent
clear a aa aExp aMat b bb ha r01Tab r10Tab rMat tMat01 tMat10 tVar1 tVar2 yExp