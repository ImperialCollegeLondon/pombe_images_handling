addpath('pombe-scope-tools');

%% MICRO 1st GATHERING
input0 = '..\Files\treated-data\Microscopy\';
list0 = dir(input0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
for aa=1:length(list0)
    fov_folder1 = list0(aa).name;
    input1      = [input0 fov_folder1];
    list1       = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
    for bb=1:length(list1)
        fov_folder2 = list1(bb).name;
        input2      = [input1 '/' fov_folder2 '/segmented/'];
        list2       = dir(input2);  list2(1:2)=[];

        count_BFok  = 0;
        cLen        = zeros(length(list2),1);
        cWid        = zeros(length(list2),1);
        widMax      = zeros(length(list2),1);
        cArea       = zeros(length(list2),1);
        count_GFPok = 0;
        nucleiN     = zeros(length(list2),1);
        nucleiD     = zeros(length(list2),1);
        septum      = zeros(length(list2),1);
        cut3Pos     = zeros(length(list2),1);
        GFP_tot     = zeros(length(list2),1);
        GFP_avg     = zeros(length(list2),1);
        GFP_med     = zeros(length(list2),1);

        for kk=1:length(list2)
            tName = load([input2 list2(kk).name]);
            tName = tName.segmented_cells;
            if abs(tName.status) ~= Inf
                count_BFok  = count_BFok+1;
                cLen(kk)    = tName.geom.l;
                cWid(kk)    = tName.geom.width_avg;
                widMax(kk)  = max(tName.geom.width_vec);
                cArea(kk)   = tName.geom.CS;
                GFP_tot(kk) = tName.gfpInt.totPerinuc;        
                GFP_avg(kk) = tName.gfpInt.meanPerinuc;
                GFP_med(kk) = tName.gfpInt.medPerinuc;
                if strcmp(tName.GFPstatus, 'perfect')
                    count_GFPok  = count_GFPok+1;
                    nucleiN(kk)  = tName.nucleiN;
                    nucleiD(kk)  = tName.nucleiD;
                    septum(kk)   = tName.septated;
                    cut3Pos(kk)  = tName.mitotic;
                end
            end
           if mod(kk,100)==0;    disp(kk);    end
        end
        %features = struct('cell_length',{}, 'cell_width',{}, 'cell_area',{}, 'nucleiN',{}, 'nucleiD',{}, 'septated',{}, 'mitotic',{});
        features.cell_length = cLen;    features.cell_width = cWid; features.cell_area  = cArea;
        features.nucleiN = nucleiN;     features.nucleiD = nucleiD; features.septated = septum;     features.mitotic = cut3Pos;
        features.GFP_tot = GFP_tot;     features.GFP_avg = GFP_avg;	features.GFP_med = GFP_med;     features.cell_MaxWidth = widMax;
        %counters = struct('BF_ok',{}, 'GFP_ok',{}, 'sept',{}, 'mit',{});
        counters.BF_ok = count_BFok;    counters.GFP_ok = count_GFPok; 
        counters.sept = sum(septum > 1.5 & septum <1000);  counters.mit = sum(cut3Pos >= 1 & cut3Pos <1000);
        save([input1 '/' fov_folder1 '_' fov_folder2],'features','counters');
        disp (input2);
    end
end


%% MICRO 2nd GATHERING
input0 = '..\Files\treated-data\Microscopy\';
list0 = dir(input0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
countMat = zeros(4,3*length(list0));
for aa=1:length(list0)
    fov_folder1 = list0(aa).name;
    input1      = [input0 fov_folder1];
    list1       = dir(input1); list1(1:2) = []; list1 = list1(~[list1(:).isdir]);
    %list1       = list1([1,4,7]);
    
    for bb=1:length(list1)
        fov_folder2 = list1(bb).name;
        xx = 3*aa +bb -3;
        names(xx) = {fov_folder2};
        tName = load([input1 '/' fov_folder2]);
        %features = struct('cell_length',{}, 'cell_width',{}, 'cell_area',{}, 'nucleiN',{}, 'nucleiD',{}, 'septated',{}, 'mitotic',{});
        features(xx) = tName.features;
        %counters = struct('BF_ok',{}, 'GFP_ok',{}, 'sept',{}, 'mit',{});
        counters(xx) = tName.counters;
        countMat(:,xx) = cell2mat(struct2cell(tName.counters));
    end
    
end
save([input0 'globStructure'],'names','features','counters');
save([input0 'globCounter'],'names','countMat');



%% GROWTH 1st GATHERING
input1 = '..\Files\treated-data\Growth\';
list1 = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
for jj=1:length(list1)
    fov_folder = list1(jj).name;
    address = [input1 fov_folder '/_segmented/'];
    list2 = dir(address);  list2(1:2)=[];

    count_BFok  = 0;
    cLen        = zeros(length(list2),1);
    cWid        = zeros(length(list2),1);
    widMax      = zeros(length(list2),1);
    cArea       = zeros(length(list2),1);
    count_GFPok = 0;
    nucleiN     = zeros(length(list2),1);
    nucleiD     = zeros(length(list2),1);
    septum      = zeros(length(list2),1);
    cut3Pos     = zeros(length(list2),1);
    GFP_tot     = zeros(length(list2),1);
    GFP_avg     = zeros(length(list2),1);
    GFP_med     = zeros(length(list2),1);

    for kk=1:length(list2)
        tName = load([address list2(kk).name]);
        tName = tName.segmented_cells;
        if abs(tName.status) ~= Inf
            count_BFok      = count_BFok+1;
            if any(strcmp(fieldnames(tName),'geom'))
                cLen(kk)    = tName.geom.l;
                cWid(kk)    = tName.geom.width_avg;
                widMax(kk)  = max(tName.geom.width_vec);
                cArea(kk)   = tName.geom.CS;
            else
                cLen(kk)    = NaN;
                cWid(kk)    = NaN;
                widMax(kk)  = NaN;
                cArea(kk)   = NaN;
            end
            GFP_tot(kk) = tName.gfpInt.totPerinuc;        
            GFP_avg(kk) = tName.gfpInt.meanPerinuc;
            GFP_med(kk) = tName.gfpInt.medPerinuc;
            if strcmp(tName.GFPstatus, 'perfect')
                count_GFPok = count_GFPok+1;
                nucleiN(kk)  = tName.nucleiN;
                if any(strcmp(fieldnames(tName),'nucleiD')) && tName.handpick==0
                        nucleiD(kk)  = tName.nucleiD;
                else;   nucleiD(kk)  = NaN;
                end                
                septum(kk)   = tName.septated;
                if any(strcmp(fieldnames(tName),'mitotic'));    cut3Pos(kk)  = tName.mitotic;
                else;                                           cut3Pos(kk)  = 0;
                end
                
            end
        end
       if mod(kk,100)==0;    disp(kk);    end
    end

    %features = struct('cell_length',{}, 'cell_width',{}, 'cell_area',{}, 'nucleiN',{}, 'nucleiD',{}, 'septated',{}, 'mitotic',{});
    features.cell_length = cLen;    features.cell_width = cWid; features.cell_area  = cArea;
    features.nucleiN = nucleiN;     features.nucleiD = nucleiD; features.septated = septum;     features.mitotic = cut3Pos;
    features.GFP_tot = GFP_tot;     features.GFP_avg = GFP_avg;	features.GFP_med = GFP_med;     features.cell_MaxWidth = widMax;
    %counters = struct('BF_ok',{}, 'GFP_ok',{}, 'sept',{}, 'mit',{});
    counters.BF_ok = count_BFok;    counters.GFP_ok = count_GFPok; 
    counters.sept = sum(septum > 1.5 & septum <1000);  counters.mit = sum(cut3Pos >= 1 & cut3Pos <1000);
    save([input1 fov_folder '/' fov_folder],'features','counters');
    disp (fov_folder);
end


%% GROWTH 2nd GATHERING
input1 = '..\Files\treated-data\Growth\';
list1 = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
countMat = zeros(4,length(list1));
for jj=1:length(list1)
    fov_folder = list1(jj).name;
    names(jj) = {fov_folder};
    tName = load([input1 fov_folder '/' fov_folder]);
    %features = struct('cell_length',{}, 'cell_width',{}, 'cell_area',{}, 'nucleiN',{}, 'nucleiD',{}, 'septated',{}, 'mitotic',{});
    features(jj) = tName.features;
    %counters = struct('BF_ok',{}, 'GFP_ok',{}, 'sept',{}, 'mit',{});
    counters(jj) = tName.counters;
    countMat(:,jj) = cell2mat(struct2cell(tName.counters));
end
save([input1 'globStructure'],'names','features','counters');
save([input1 'globCounter'],'names','countMat');



%% DYE 1st GATHERING: NEWEST - USES BEST DNA INTENSITIES
input0 = '..\Files\treated-data\Dye\';
list0 = dir(input0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
for mm=1:length(list0)
    temp1   = list0(mm).name;
    input1  = [input0 temp1 '/'];
    list1   = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
    
    for jj=1:length(list1)
        fov_folder = list1(jj).name;
        address = [input1 fov_folder '/_segmented/'];
        list2 = dir(address);  list2(1:2)=[];
        count_void = 0;

        count_BFok  = 0;
        cLen        = zeros(length(list2),1);
        cWid        = zeros(length(list2),1);
        widMax      = zeros(length(list2),1);
        cArea       = zeros(length(list2),1);

        count_GFPok = 0;
        nucleiN     = zeros(length(list2),1);
        nucleiD     = zeros(length(list2),1);
        septum      = zeros(length(list2),1);
        GFP_tot     = zeros(length(list2),1);
        GFP_avg     = zeros(length(list2),1);
        GFP_med     = zeros(length(list2),1);

        count_REDok = 0;
        dLen        = zeros(length(list2),1);
        dArea       = zeros(length(list2),1);
        rArea       = zeros(length(list2),1);
        rLen        = zeros(length(list2),1);
        delta       = zeros(length(list2),1);
        
        redINcell   = zeros(length(list2),1);
        redOUTcell  = zeros(length(list2),1);
        redINdye    = zeros(length(list2),1);
        redOUTdye   = zeros(length(list2),1);
        
        for kk=1:length(list2)
            tName = load([address list2(kk).name]);
            tName = tName.segmented_cells;
            if abs(tName.status) ~= Inf
                count_BFok  = count_BFok+1;
                if any(strcmp(fieldnames(tName),'geom'))
                    cLen(kk)    = tName.geom.l;
                    cWid(kk)    = tName.geom.width_avg;
                    widMax(kk)  = max(tName.geom.width_vec);
                    cArea(kk)   = tName.geom.CS;
                else
                    cLen(kk)    = NaN;
                    cWid(kk)    = NaN;
                    widMax(kk)  = NaN;
                    cArea(kk)   = NaN;
                end
                GFP_tot(kk) = tName.gfpInt.totPerinuc;        
                GFP_avg(kk) = tName.gfpInt.meanPerinuc;
                GFP_med(kk) = tName.gfpInt.medPerinuc;
            
                if strcmp(tName.GFPstatus, 'perfect')
                    count_GFPok = count_GFPok+1;
                    nucleiN(kk) = tName.nucleiN;
                    if any(strcmp(fieldnames(tName),'nucleiD'))  && tName.handpick==0
                            nucleiD(kk)  = tName.nucleiD;
                    else;   nucleiD(kk)  = NaN;
                    end                
                    septum(kk)  = tName.septated;
                    if any(strcmp(fieldnames(tName),'REDstatus'))
                        if strcmp(tName.REDstatus, 'perfect')
                            count_REDok = count_REDok+1;
                            dLen(kk)    = tName.REDgeom.l;
                            dArea(kk)   = tName.REDgeom.area;
                            rLen(kk)    = tName.REDgeom.ratio_length;
                            rArea(kk)   = tName.REDgeom.ratio_area;
                            delta(kk)   = tName.REDgeom.deltaYC;

                            redINcell(kk)  = tName.REDint.INcell;
                            redOUTcell(kk)  = tName.REDint.OUTcell;
                            redINdye(kk)   = tName.REDint.INdye;
                            redOUTdye(kk)   = tName.REDint.OUTdye;
                        end
                    end
                end
            end
           if mod(kk,100)==0;    disp(kk);    end
        end

        %features = struct('cell_length',{}, 'cell_width',{}, 'cell_area',{}, 'nucleiN',{}, 'nucleiD',{}, 'septated',{});
        features.cell_length = cLen;    features.cell_width = cWid;     features.cell_area  = cArea;
        features.nucleiN = nucleiN;     features.nucleiD = nucleiD;     features.septated = septum;
        features.dye_length = dLen;     features.dye_area = dArea;      features.dye_delta = delta;
        features.ratio_area = rArea;    features.ratio_length = rLen;
        features.GFP_tot = GFP_tot;     features.GFP_avg = GFP_avg;     features.GFP_med = GFP_med;
        features.RED_out = redOUTcell;  features.RED_cel = redINcell;   features.RED_noDye = redOUTdye;     features.RED_dye = redINdye;
        features.cell_MaxWidth = widMax;
        %counters = struct('BF_ok',{}, 'GFP_ok',{}, 'RED_ok',{},'sept',{});
        counters.BF_ok = count_BFok;    counters.GFP_ok = count_GFPok;  counters.RED_ok = count_REDok;
        counters.void = count_void;     counters.sept = sum(septum > 1.5 & septum <1000);
        save([input1 fov_folder '/' fov_folder],'features','counters');
        disp (fov_folder);
        disp(['tot= ',num2str(length(list2)), ' void= ',num2str(count_void), ' BFok= ',num2str(count_BFok), ...
             ' GFPok= ',num2str(count_GFPok), ' REDok= ',num2str(count_REDok)]);
    end
end


%% DYE 2nd GATHERING
input0 = '..\Files\treated-data\Dye\';
list0 = dir(input0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
for mm=1:length(list0)
    temp1   = list0(mm).name;
    input1  = [input0 temp1 '/'];
    list1   = dir(input1); list1(1:2) = []; list1 = list1([list1(:).isdir]);
    countMat = zeros(5,length(list1));
    
    backbone    = load([input1 list1(1).name '/' list1(1).name]);
    features    = struct(backbone.features);   
    counters    = struct(backbone.counters);        
    names       = cell(1,length(list1));
    for jj=1:length(list1)
        fov_folder = list1(jj).name;
        names(jj) = {fov_folder};
        tName = load([input1 fov_folder '/' fov_folder]);
        %features = struct('cell_length',{}, 'cell_width',{}, 'cell_area',{}, 'nucleiN',{}, 'nucleiD',{}, 'septated',{});
        features(jj) = tName.features;
        %counters = struct('BF_ok',{}, 'GFP_ok',{}, 'RED_ok',{},'sept',{});
        counters(jj) = tName.counters;
        countMat(:,jj) = cell2mat(struct2cell(tName.counters));
    end
    save([input1 temp1 '_globStructure'],'names','features','counters');
    save([input1 temp1 '_globCounter'],'names','countMat');
end


%% DYE 3rd GATHERING
input0 = '..\Files\treated-data\Dye\';
list0 = dir(input0); list0(1:2) = []; list0 = list0([list0(:).isdir]);
names = []; countMat = []; features = []; counters = [];
for jj=1:length(list0)
    fov_folder = list0(jj).name;
    tempStr = load([input0 fov_folder '/' fov_folder '_globStructure']);
    tempMat = load([input0 fov_folder '/' fov_folder '_globCounter']);
    
    names    = horzcat(names,tempMat.names);
    countMat = horzcat(countMat,tempMat.countMat);
    features = horzcat(features,tempStr.features);
    counters = horzcat(counters,tempStr.counters);

end
save([input0 'globStructure'],'names','features','counters');
save([input0 'globCounter'],'names','countMat');


