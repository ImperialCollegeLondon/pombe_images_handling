%% MULTIPURPOSE AUXILIARY PART


%READING INPUTS
run('..\Scripts - auxiliary\aux_dye_rateCalc.m')    %extracting data from OD measurements
input1  = '..\Files\treated-data\Dye\';
load([input1 'globStructure'],'names','features','counters');
namDye  = names;
featDye = features;
featNamDye = fieldnames(featDye);
for aa=1:length(featNamDye);    featNamDye(aa)= strrep(featNamDye(aa),'_',' '); end


% AUTOMATIC PROCESSING
conds8 = namDye;     shRepl = conds8;    shTPs = shRepl;    TPonly = zeros(1,length(namDye));
for aa=1:length(namDye)
    namDye(aa)  = strrep(namDye(aa),'_',' ');
    namDye(aa)  = strrep(namDye(aa),'157 ','');
    namDye(aa)  = strrep(namDye(aa),'NH4Cl','NH4');
    namDye(aa)  = strrep(namDye(aa),'r10','B');     %inversion, 29/7/2018
    namDye(aa)  = strrep(namDye(aa),'r1','A');
    tVar1       = cell2mat(namDye(aa));
    conds8(aa)  = {tVar1(1:3)};
    shRepl(aa)  = {['D ' tVar1(1:5)]};
    shTPs(aa)   = {tVar1(5:end)};
    TPonly(aa)  = str2double(tVar1(8:end));
end
conds8      = unique(conds8)';
conds16     = unique(shRepl)';
namDye      = namDye';


% CORRECTING FOR SCRAMBLED FILES during files gathering
tVar1=[];    mSize = [];
for aa=1:length(conds16)
    tVar2 = ismember(shRepl, conds16(aa));
    [~,b] = sort(TPonly(tVar2));
    tVar1 = [tVar1, (find(tVar2,1)+b-1)];
    mSize = [mSize,1:sum(tVar2)];
end
counters = counters(tVar1); featDye=featDye(tVar1);
namDye   = namDye(tVar1); shRepl=shRepl(tVar1); shTPs=shTPs(tVar1); TPonly=TPonly(tVar1);

indHA2 = zeros(1,16); indSize = indHA2;
for aa=1:length(conds16)
    tVar2 = ismember(shRepl, conds16(aa));
    indHA2(aa) = find(tVar2,1);   
    indSize(aa)= sum(tVar2);
end
indHA2 = indHA2+5;
condsDye  = repelem(conds16,indSize);


% REMOVING WRONGLY PROCESSED IMAGES
featStr         = struct2cell(featDye);    featDye = cell(length(namDye),1);
for aa=1:size(featStr,3)
    tVar1 = cell2mat(featStr(:,:,aa)');     
    tVar1 = tVar1(tVar1(:,4)~=0,:);
    featDye(aa) = {tVar1};
end