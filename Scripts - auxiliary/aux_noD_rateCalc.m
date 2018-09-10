%% GROWTH RATES
%extracting data from turbidostat
%inputPath = '..\..\..\Cycle sensor data\0. Turbidostat\';
%addpath([inputPath 'Analysis code']);
%input1 = [inputPath 'Traces\'];
%list1 = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
%grVecNoD = zeros(length(list1),1); %figure;
%for aa=1:length(list1)
%        turbo = load_turbi_data([input1 list1(aa).name],0);
%        if size(turbo.growth.alpha_growth,1)
%            grVecNoD(aa)= mean(turbo.growth.alpha_growth);
%        end
%        %subplot(4,6,aa);plot(turbo.t,turbo.OD);
%        %subplot(4,6,aa);plot(tempTurbo,'LineStyle','--','Marker','s');
%end
%grVecNoD = real(grVecNoD);

%growth rates already extracted and checked, no need to make them again!
grVecNoD =   [0.0403;0.0421;0.0600;  0.1239;0.1474;0.1045;  0.2475;0.2372;0.1961;...
             0.0894;0.0131;0.0713;  0.1513;0.1386;0.1406;  0.0472;0.1073;0.0601;...
             0.0342;0.0470;0.0543;  0.0498;0.0385;0.0587];
grVecNoD([11,12, 17,18, 19,20,21]) = grVecNoD([12,11, 18,17, 20,21,19]);    %swapping per matchare swapped files


grMean=zeros(1,length(grVecNoD)/3);    grStd=grMean;   grFlagNoD=false(length(grVecNoD),1);
for aa=1:length(grMean)    
    %grMean(aa) = mean(grVecNoD(3*aa-2 : 3*aa));
    grMean(aa) = median(grVecNoD(3*aa-2 : 3*aa));
    grStd(aa)  = std(grVecNoD(3*aa-2 : 3*aa));
end
for aa=1:length(grVecNoD)
    tVar1 = (grVecNoD(aa) - grMean(ceil(aa/3)))/grStd(ceil(aa/3)); 
    if abs(tVar1)>1.1; grFlagNoD(aa) = true; end 
end

%grFlagNoD(:) = 0; grFlagNoD([3,12,16,17,18,21,24])=1;
clear aa grMean grStd input1 inputPath list1 turbo tVar1
