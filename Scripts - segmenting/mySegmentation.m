%% CALLING FUNCTION
function segmented_cells = mySegmentation(fov,dye,do_example_plots)
    switch dye
        case 'dye';     segmented_cells = dyeSegmentation(fov,do_example_plots,spotRad);
        case 'turbi';   segmented_cells = growthSegmentation(fov, do_example_plots,spotRad);
        case 'micro';   segmented_cells = microSegmentation(fov, do_example_plots);            
    end
end


%% GENERAL FUNCTIONS
%working on the dye  as well, and faster!
function segmented_cells = dyeSegmentation(fov,do_example_plots,spotRad)
    %segmented_cells = struct('mask',{},'area',{},'geom',{},'status',{},...
    %                       'gfp_per_pixel',{},'red_per_pixel',{},'gfp_local_bg',{},'red_local_bg',{},...
    %                       'GFPmask',{},'nucleiN',{},'nucleiD',{},'GFPstatus',{});
    BF_cell  = fov.loadImage(1,1,'BF', 'normalize_contrast');
    segmented_cells = adTSegmentation(BF_cell,do_example_plots);
    segmented_cells(1).handpick = 0;
    segmented_cells.gfpInt.totCell          = nan;
    segmented_cells.gfpInt.totNonCell       = nan;
    segmented_cells.gfpInt.meanCell         = nan;
    segmented_cells.gfpInt.meanNonCell      = nan;
    segmented_cells.gfpInt.medCell          = nan;
    segmented_cells.gfpInt.medNonCell       = nan;
    segmented_cells.gfpInt.totPerinuc       = nan;
    segmented_cells.gfpInt.totNonPerinuc    = nan;
    segmented_cells.gfpInt.meanPerinuc      = nan;
    segmented_cells.gfpInt.meanNonPerinuc   = nan;
    segmented_cells.gfpInt.medPerinuc       = nan;
    segmented_cells.gfpInt.medNonPerinuc    = nan;
    segmented_cells(1).REDint.INcell        = nan;    
    segmented_cells(1).REDint.OUTcell       = nan;    
    segmented_cells(1).REDint.INdye         = nan;
    segmented_cells(1).REDint.OUTdye        = nan;
    if segmented_cells(1).status == Inf;   return;     end 
    
    GFP_cell0 = fov.loadImage(1,1,'GFP','none');
    segmented_cells.gfpInt.totCell       = sum(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.meanCell      = mean(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.medCell       = median(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.totNonCell    = sum(GFP_cell0(not(segmented_cells.mask)));
    segmented_cells.gfpInt.meanNonCell   = mean(GFP_cell0(not(segmented_cells.mask)));
    segmented_cells.gfpInt.medNonCell    = median(GFP_cell0(not(segmented_cells.mask)));
    
    Red_cell = fov.loadImage(1,1,'Red','none');
    Red_cell = Operations.normalize(Red_cell,'normalize_contrast');
    segmented_cells(1).REDint.INcell    = mean( Red_cell(segmented_cells.mask) );
    segmented_cells(1).REDint.OUTcell   = median( Red_cell(not(segmented_cells.mask)) );
    
    GFP_cell = Operations.normalize(GFP_cell0,'normalize_contrast');
    segmented_cells = gfpSegmentation(GFP_cell,segmented_cells,do_example_plots,spotRad);
    if ~strcmp(segmented_cells(1).GFPstatus,'perfect');    return; end
    
    if (segmented_cells.nucleiN==2) && (segmented_cells.septated==0) && (segmented_cells.mitotic==0)
        segmented_cells = septum_finder(GFP_cell,segmented_cells,false);
    end
     
    segmented_cells = redSegmentation(Red_cell,segmented_cells,do_example_plots);
    if ~strcmp(segmented_cells(1).REDstatus,'perfect');   return; end
    segmented_cells(1).REDint.OUTdye    = mean( Red_cell(segmented_cells.mask & not(segmented_cells.REDmask)));
    segmented_cells(1).REDint.INdye     = mean( Red_cell(segmented_cells.REDmask));
    segmented_cells = ratio_finder(segmented_cells);
end

function segmented_cells = growthSegmentation(fov,do_example_plots,spotRad)
    %segmented_cells = struct('mask',{},'area',{},'geom',{},'status',{},...
    %                       'gfp_per_pixel',{},'red_per_pixel',{},'gfp_local_bg',{},'red_local_bg',{},...
    %                       'GFPmask',{},'nucleiN',{},'nucleiD',{},'GFPstatus',{});
    BF_cell  = fov.loadImage(1,1,'BF', 'normalize_contrast');
    segmented_cells = adTSegmentation(BF_cell,do_example_plots);
    segmented_cells(1).handpick = 0;
    segmented_cells.gfpInt.totCell          = nan;
    segmented_cells.gfpInt.totNonCell       = nan;
    segmented_cells.gfpInt.meanCell         = nan;
    segmented_cells.gfpInt.meanNonCell      = nan;
    segmented_cells.gfpInt.medCell          = nan;
    segmented_cells.gfpInt.medNonCell       = nan;
    segmented_cells.gfpInt.totPerinuc       = nan;
    segmented_cells.gfpInt.totNonPerinuc    = nan;
    segmented_cells.gfpInt.meanPerinuc      = nan;
    segmented_cells.gfpInt.meanNonPerinuc   = nan;
    segmented_cells.gfpInt.medPerinuc       = nan;
    segmented_cells.gfpInt.medNonPerinuc    = nan;
    if segmented_cells(1).status == Inf;   return;     end

    GFP_cell0   = fov.loadImage(1,1,'GFP','none');
    segmented_cells.gfpInt.totCell       = sum(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.meanCell      = mean(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.medCell       = median(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.totNonCell    = sum(GFP_cell0(not(segmented_cells.mask)));
    segmented_cells.gfpInt.meanNonCell   = mean(GFP_cell0(not(segmented_cells.mask)));
    segmented_cells.gfpInt.medNonCell    = median(GFP_cell0(not(segmented_cells.mask)));

    Red_cell    = fov.loadImage(1,1,'Red','none');
    segmented_cells(1).red_per_pixel    = mean( Red_cell(segmented_cells.mask) );
    segmented_cells(1).red_local_bg     = double(median( Red_cell(:) )); 
    
    GFP_cell = Operations.normalize(GFP_cell0,'normalize_contrast');
    segmented_cells = gfpSegmentation(GFP_cell,segmented_cells,do_example_plots,spotRad);
    if ~strcmp(segmented_cells(1).GFPstatus,'perfect');    return; end

    if (segmented_cells.nucleiN==2) && (segmented_cells.septated==0) && (segmented_cells.mitotic==0)
        segmented_cells = septum_finder(GFP_cell,segmented_cells,false);
    end
    
    Red_cell = Operations.normalize(Red_cell,'normalize_contrast');
    if (segmented_cells.septated==0) && (segmented_cells.mitotic==0)
        segmented_cells = cut3_finder(Red_cell,segmented_cells,false);
    end
end

function segmented_cells = microSegmentation(fov,do_example_plots)
    %segmented_cells = struct('mask',{},'area',{},'geom',{},'status',{},...
    %                       'gfp_per_pixel',{},'red_per_pixel',{},'gfp_local_bg',{},'red_local_bg',{},...
    %                       'GFPmask',{},'nucleiN',{},'nucleiD',{},'GFPstatus',{});
    BF_cell  = fov.loadImage(1,1,'BF', 'normalize_contrast');
    segmented_cells = adTSegmentation(BF_cell,do_example_plots);
    segmented_cells(1).handpick = 0;
    segmented_cells.gfpInt.totCell          = nan;
    segmented_cells.gfpInt.totNonCell       = nan;
    segmented_cells.gfpInt.meanCell         = nan;
    segmented_cells.gfpInt.meanNonCell      = nan;
    segmented_cells.gfpInt.medCell          = nan;
    segmented_cells.gfpInt.medNonCell       = nan;
    segmented_cells.gfpInt.totPerinuc       = nan;
    segmented_cells.gfpInt.totNonPerinuc    = nan;
    segmented_cells.gfpInt.meanPerinuc      = nan;
    segmented_cells.gfpInt.meanNonPerinuc   = nan;
    segmented_cells.gfpInt.medPerinuc       = nan;
    segmented_cells.gfpInt.medNonPerinuc    = nan;
    if segmented_cells(1).status == Inf;   return;     end
      
    Red_cell = fov.loadImage(1,1,'Red','none');
    GFP_cell0 = fov.loadImage(1,1,'GFP','none');
    segmented_cells.gfpInt.totCell       = sum(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.meanCell      = mean(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.medCell       = median(GFP_cell0(segmented_cells.mask));
    segmented_cells.gfpInt.totNonCell    = sum(GFP_cell0(not(segmented_cells.mask)));
    segmented_cells.gfpInt.meanNonCell   = mean(GFP_cell0(not(segmented_cells.mask)));
    segmented_cells.gfpInt.medNonCell    = median(GFP_cell0(not(segmented_cells.mask)));
    
    GFP_cell = Operations.normalize(GFP_cell0,'normalize_contrast');
    segmented_cells = gfpMicroSegmentation(GFP_cell,segmented_cells,do_example_plots);
    if ~strcmp(segmented_cells(1).GFPstatus,'perfect');    return; end
    
    segmented_cells(1).gfpINcell    = sum(GFP_cell0(segmented_cells.BFonGFPmask));
    if (segmented_cells.nucleiN==2) && (segmented_cells.septated==0) && (segmented_cells.mitotic==0)
        segmented_cells = septumMicro_finder(GFP_cell,segmented_cells,false);
    end
    
    Red_cell = Operations.normalize(Red_cell,'normalize_contrast');
    if (segmented_cells.septated==0) && (segmented_cells.mitotic==0)
        segmented_cells = cut3Micro_finder(Red_cell,segmented_cells,false);
    else;   segmented_cells(1).CUT3status='perfect';
    end
end



%% SEGMENTING FUNCTIONS
function segmented_cells = adTSegmentation(BF_cell,do_example_plots)
    Cs = [0.025,0.02,0.03,0.01,0.04];
    
    geom.success = false; atts=0;
    while ((~geom.success) && atts<6)
        atts = atts+1;
        if atts<6
            cells_interiors = bwareaopen(imfill(imclearborder(adaptivethresholdSH(BF_cell,10,Cs(atts))),'holes'),100);
        else
            cells_interiors = bwareaopen(imfill(imclearborder(im2bw(BF_cell,graythresh(BF_cell))),'holes'),100);            
        end
        if sum(cells_interiors(:)) > 0
            cell_mask = Operations.getCentralCC(cells_interiors);   
            props = regionprops(cell_mask,'Area','Solidity','Centroid','PixelIdxList');
            if (props.Area > 100 && props.Area < 10000 && props.Solidity > 0.8)
                geom = Operations.mask2geom(cell_mask,60,false);
                if geom.success
                    delta  = abs(geom.width_vec(45)-geom.width_vec(15));
                    check1 = delta/geom.width_avg;  check2 = geom.l*geom.width_avg/geom.CS;
                    %if or(or( ((delta>4) && (check1>0.5)),(max(segmented_cells.geom.width_vec)>15)),check2>1.25);
                    if ~(geom.width_avg<20 && check2<=1.25) || ((delta>4) && (check1>0.5)) %new, only needed for microscopy
                        geom.success = false;
                    end        
                end
            end
        end
    end
    
    if geom.success
        segmented_cells(1).area     = props.Area;
        segmented_cells(1).centroid = props.Centroid;
        segmented_cells(1).mask     = imerode(imdilate(cell_mask,strel('disk',6)),strel('disk',6));
        segmented_cells(1).geom     = geom;
        segmented_cells(1).status   = atts;
        if do_example_plots
            figure(1);
            cell_boundary = bwboundaries(cell_mask);
            cell_boundary = cell_boundary{1};
            Operations.imshowfit([BF_cell,cell_mask]); hold on;
            plot( cell_boundary(:,2) , cell_boundary(:,1) , 'r');
        end        
    else
        disp('BF failure');
        segmented_cells(1).status = Inf;
        if do_example_plots
            figure(3); Operations.imshowfit([BF_cell,cell_mask]); hold on;
        end
    end
end

function segmented_cells = gfpSegmentation(GFP_cell,segmented_cells,do_example_plots,spotRad)
    thresh      = graythresh(GFP_cell);   atts=0;
    GFP_thresh  = im2bw(GFP_cell,thresh);
    cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
    cells_interiors = cells_interiors & segmented_cells.mask;
    
    if sum(cells_interiors(:)) > 0
        segmented_cells(1).mitotic = 0;
        segmented_cells(1).septated= 0;
        segmented_cells(1).GFPstatus = false;
        while (~strcmp(segmented_cells(1).GFPstatus,'perfect') && atts<3)
            GFPdot      = false(size(GFP_cell));
            atts  = atts+1;
            props = bwconncomp(cells_interiors);
            tVec  = zeros(props.NumObjects,1);
            for jj=1:props.NumObjects; tVec(jj) = length(props.PixelIdxList{jj});  end

            while any(tVec>50) && (thresh<0.9)
                threshSt    = thresh;
                thresh      = thresh+0.1;
                GFP_thresh  = im2bw(GFP_cell,thresh);
                cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
                cells_interiors = cells_interiors & segmented_cells.mask;
                props = bwconncomp(cells_interiors);
                tVec  = zeros(props.NumObjects,1);
                for jj=1:props.NumObjects; tVec(jj) = length(props.PixelIdxList{jj});  end
                if all(tVec<5)
                    tVec   = Inf;
                    thresh = threshSt;
                    break
                end
            end
            while any(tVec>50) && (thresh<0.99)
                threshSt    = thresh;
                thresh      = thresh+0.01;
                GFP_thresh  = im2bw(GFP_cell,thresh);
                cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
                cells_interiors = cells_interiors & segmented_cells.mask;
                props = bwconncomp(cells_interiors);
                tVec  = zeros(props.NumObjects,1);
                for jj=1:props.NumObjects; tVec(jj) = length(props.PixelIdxList{jj});  end
                if all(tVec<5)
                    tVec   = Inf;
                    thresh = threshSt;
                    break
                end
            end
            while any(tVec>50) && (thresh<0.999)
                thresh      = thresh+0.001;
                GFP_thresh  = im2bw(GFP_cell,thresh);
                cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
                cells_interiors = cells_interiors & segmented_cells.mask;
                props = bwconncomp(cells_interiors);
                tVec  = zeros(props.NumObjects,1);
                for jj=1:props.NumObjects; tVec(jj) = length(props.PixelIdxList{jj});  end
            end
            props = regionprops(cells_interiors , 'Area','Solidity', 'Centroid' , 'PixelIdxList','MajorAxisLength','MinorAxisLength');
            
            if length(props)==1     %check for hidden 2nd spots             
                bfCC = segmented_cells.centroid(2); gfpCC = props.Centroid(2); bfL = segmented_cells.geom.l;
                if  or((gfpCC <  bfCC - 0.1*bfL),(gfpCC > bfCC + 0.1*bfL))
                    if (gfpCC < bfCC - 0.1*bfL);        indX = round(bfCC):size(GFP_cell,1);
                    elseif (gfpCC > bfCC + 0.1*bfL);    indX = 1:round(bfCC);   end
                    GFP_cell2 = GFP_cell(indX,:);
                    thresh    = graythresh(GFP_cell2);
                    GFP_thresh2 = im2bw(GFP_cell2,thresh);
                    GFP_thresh2 = imerode(GFP_thresh2,strel('disk',1));
                    GFP_thresh(indX,:) = GFP_thresh2; 
                    cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
                    cells_interiors = cells_interiors & segmented_cells.mask;
                    props = regionprops(cells_interiors , 'Area','Solidity', 'Centroid' , 'PixelIdxList','MajorAxisLength','MinorAxisLength');
                end
            end
            
            
            if      length(props)==1
                area = props.Area;  solid= props.Solidity;
                nucD = props.MajorAxisLength;
                ratio= nucD/props.MinorAxisLength;
                if solid>=0.8 && area<100 && ~xor((area>50),(ratio>2))
                    segmented_cells(1).nucleiN = 1;
                    segmented_cells(1).GFPmask = props.PixelIdxList;
                    segmented_cells(1).nucleiD = nucD;
                    segmented_cells(1).GFPstatus        = 'perfect';
                    if  ratio>2; segmented_cells(1).mitotic = 1;    end
                    tCen = round(props(1).Centroid);
                    GFPdot(tCen(1),tCen(2)) = 1;
                end
            elseif  length(props)==2
                area  = max(props.Area);  solid= min(props.Solidity);
                nucD  = props(1).MajorAxisLength/2 + props(2).MajorAxisLength/2 + sqrt(sum((props(1).Centroid - props(2).Centroid).^2));
                %ratio = nucD/max(props.MinorAxisLength);
                if solid>=0.8 && area<100
                    segmented_cells(1).nucleiN = 2;
                    segmented_cells(1).GFPmask = props(1).PixelIdxList;
                    segmented_cells(1).GFPmask2= props(2).PixelIdxList;
                    segmented_cells(1).nucleiD = nucD;
                    segmented_cells(1).GFPstatus        = 'perfect';
                    if nucD/segmented_cells.geom.l <0.33
                        segmented_cells(1).mitotic = 1;
                    end
                    tCen = round(props(1).Centroid);
                    GFPdot(tCen(1),tCen(2)) = 1;
                    tCen = round(props(2).Centroid);
                    GFPdot(tCen(1),tCen(2)) = 1;                    
                end     
            else;    segmented_cells(1).septated= 10;
            end
            if any(any(GFPdot)) && not(spotRad==0)
                GFPdot = imdilate(GFPdot,strel('disk',spotRad));
                segmented_cells(1).gfpInt.totPerinuc     = sum(GFP_cell(segmented_cells.mask & GFPdot));
                segmented_cells(1).gfpInt.totNonPerinuc  = sum(GFP_cell(segmented_cells.mask & not(GFPdot)));
                segmented_cells(1).gfpInt.meanPerinuc    = mean(GFP_cell(segmented_cells.mask & GFPdot));
                segmented_cells(1).gfpInt.meanNonPerinuc = mean(GFP_cell(segmented_cells.mask & not(GFPdot)));
                segmented_cells(1).gfpInt.medPerinuc     = median(GFP_cell(segmented_cells.mask & GFPdot));
                segmented_cells(1).gfpInt.medNonPerinuc  = median(GFP_cell(segmented_cells.mask & not(GFPdot)));
            elseif any(any(GFPfdot)) && spotRad==0
                segmented_cells.gfpInt.totPerinuc       = sum(GFP_cell(cells_interiors));
                segmented_cells.gfpInt.totNonPerinuc    = sum(GFP_cell(segmented_cells.mask & not(cells_interiors)));
                segmented_cells.gfpInt.meanPerinuc      = mean(GFP_cell(cells_interiors));
                segmented_cells.gfpInt.meanNonPerinuc   = mean(GFP_cell(segmented_cells.mask & not(cells_interiors)));
                segmented_cells.gfpInt.medPerinuc       = median(GFP_cell(cells_interiors));
                segmented_cells.gfpInt.medNonPerinuc    = median(GFP_cell(segmented_cells.mask & not(cells_interiors)));                
            end
            
            if (~strcmp(segmented_cells(1).GFPstatus,'perfect') && atts<3)
                for i=1:length(props);  ycoord(i) = props(i).Centroid(2);   end
                if atts==1; midX = round(segmented_cells.centroid(2)); 
                else;       midX = round((min(ycoord)+max(ycoord))/2);
                end
                GFP_cell2 = GFP_cell;   GFP_cell2((midX-2):(midX+2),:) = 0;
                thresh    = graythresh(GFP_cell2);
                GFP_thresh = im2bw(GFP_cell2,thresh);
                cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
            end
        end
    end
    
    if sum(cells_interiors(:)) == 0
        disp('  cannot get GFP mask');
        segmented_cells(1).GFPstatus = 'thresh_fail';
        if do_example_plots;    figure(2); Operations.imshowfit(BF_cell); hold on;  end
    elseif strcmp(segmented_cells(1).GFPstatus,'perfect')
        if do_example_plots
            figure(1);
            cell_boundary = bwboundaries(cell_mask);
            cell_boundary = cell_boundary{1};
            Operations.imshowfit([BF_cell,cell_mask]); hold on;
            plot( cell_boundary(:,2) , cell_boundary(:,1) , 'r');
        end
    else
        disp(['  GFP mask not conform (area = ' num2str(props(1).Area) ' , solidity = ' num2str(props(1).Solidity) ')']);
        segmented_cells(1).GFPstatus = 'mask_fail';
        if do_example_plots;    figure(3); Operations.imshowfit([BF_cell,cell_mask]); hold on;  end
    end      
end

function segmented_cells = gfpMicroSegmentation(GFP_cell1,segmented_cells,do_example_plots)
    % very long part for finding cell mask on GFP images
    extVar      = segmented_cells.area;
    BFthr       = graythresh(GFP_cell1);
    GFP_thresh  = im2bw(GFP_cell1,BFthr);
    cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),50);
    if sum(cells_interiors(:))==0
        tMin=Inf;
    else
        props = bwconncomp(cells_interiors);
        tMin = min(abs(cellfun(@numel,props.PixelIdxList) - extVar));
    end
   
    if tMin>50
        tInt = 0:0.1:1;
        [tMin,tInd] = GFPmicro_aux1(GFP_cell1,tInt,extVar);
        if tMin>50
            tInt = tInt(max(tInd-1,1)):0.02: tInt(min(tInd+1,length(tInt)));
            [tMin,tInd] = GFPmicro_aux1(GFP_cell1,tInt,extVar);
            if tMin>50
                tInt = tInt(max(tInd-1,1)):0.004:tInt(min(tInd+1,length(tInt)));
                [tMin,tInd] = GFPmicro_aux1(GFP_cell1,tInt,extVar);
                if tMin>50
                    segmented_cells(1).GFPstatus = 'off target';
                    return
                end
            end
        end
        BFthr       = tInt(tInd);
        GFP_thresh  = im2bw(GFP_cell1,BFthr);
        cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),50);
    end
    cellONgfp = regionprops(cells_interiors,'Area','Centroid','PixelIdxList','Orientation');
    
    if length(cellONgfp)>1
        tVar = abs([cellONgfp.Area]-extVar);    %area diff and relative importance
        [~,tDelI]    = min(tVar);    tDelR = diff(abs(tVar))/max(tVar);    
        tVar = inf(1,length(cellONgfp));        %centroid diff and relative importance
        for ii=1:length(cellONgfp)
            c1 = cellONgfp(ii).Centroid./size(GFP_cell1); c2 = segmented_cells.centroid./size(segmented_cells.mask);
            tVar(ii) = sqrt(sum((c1-c2).^2));
        end
        [~,tDisI] = min(tVar);      tDisR = diff(abs(tVar))/max(tVar);
        if tDelR>tDisR; tInd = tDelI;   else;    tInd = tDisI;   end
        cellONgfp   = cellONgfp(tInd);
    end
    brightMask  = false(size(cells_interiors));
    brightMask(cellONgfp.PixelIdxList) = 1;
    cells_interiors = cells_interiors & brightMask;  
    
    
    % standard part for finding GFP mask on GFP images
    segmented_cells(1).mitotic = 0;
    segmented_cells(1).septated= 0;
    segmented_cells(1).GFPstatus = false;   lastHope = 0;
    while (~strcmp(segmented_cells(1).GFPstatus,'perfect') && lastHope<2)
        GFPdot      = false(size(GFP_cell));
        tInt = BFthr:0.1:1;
        [~,tInt] = GFPmicro_aux2(GFP_cell1,tInt,0.01);
        [~,tInt] = GFPmicro_aux2(GFP_cell1,tInt,0.001);
        [tInd,~] = GFPmicro_aux2(GFP_cell1,tInt,0);
        if tInd==Inf
            props = [];
        else
            GFPthr = tInt(tInd);
            GFP_thresh      = im2bw(GFP_cell1,GFPthr);
            GFP_thresh      = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
            GFP_thresh      = imerode(imdilate(GFP_thresh,strel('disk',1)),strel('disk',1));
            cells_interiors = GFP_thresh & brightMask;
            props = regionprops(cells_interiors , 'Area','Solidity', 'Centroid' , 'PixelIdxList','MajorAxisLength','MinorAxisLength');            
        end

        if length(props)==1     %check for hidden 2nd spots
            tDel = sqrt(sum((cellONgfp.Centroid - props.Centroid).^2));
            if  tDel/segmented_cells.geom.l>0.1
                tMask       = imdilate(cells_interiors,strel('disk',2));
                GFP_cell2   = GFP_cell1; GFP_cell2(tMask) = BFthr;
                tInt = BFthr:0.1:GFPthr;
                [~,tInt] = GFPmicro_aux2(GFP_cell2,tInt,0.01);
                [~,tInt] = GFPmicro_aux2(GFP_cell2,tInt,0.001);
                [tInd,~] = GFPmicro_aux2(GFP_cell2,tInt,0);
                if tInd~=Inf
                    GFPthr2         = tInt(tInd);
                    GFP_thresh2     = im2bw(GFP_cell2,GFPthr2);
                    GFP_thresh2     = bwareaopen(imfill(imclearborder(GFP_thresh2),'holes'),10);
                    GFP_thresh2     = imerode(imdilate(GFP_thresh2,strel('disk',1)),strel('disk',1));
                    cells_interiors = (GFP_thresh | GFP_thresh2) & brightMask;
                    props = regionprops(cells_interiors , 'Area','Solidity', 'Centroid' , 'PixelIdxList','MajorAxisLength','MinorAxisLength');
                end                
            end
        end

        if      length(props)==1
            nucD = props.MajorAxisLength;
            ratio= nucD/props.MinorAxisLength;
            segmented_cells(1).nucleiN = 1;
            segmented_cells(1).BFonGFPmask  = cellONgfp.PixelIdxList;
            segmented_cells(1).BFonGFPori   = cellONgfp.Orientation;
            segmented_cells(1).BFonGFPcen   = cellONgfp.Centroid;
            segmented_cells(1).GFPcentroid  = props.Centroid;
            segmented_cells(1).GFPmask      = props.PixelIdxList;
            segmented_cells(1).nucleiD      = nucD;
            segmented_cells(1).GFPstatus    = 'perfect';
            if  ratio>2; segmented_cells(1).mitotic = 1;    end
            tCen = round(props(1).Centroid);
            GFPdot(tCen(1),tCen(2)) = 1;         
        elseif  length(props)==2
            nucD  = props(1).MajorAxisLength/2 + props(2).MajorAxisLength/2 + sqrt(sum((props(1).Centroid - props(2).Centroid).^2));
            segmented_cells(1).nucleiN = 2;
            segmented_cells(1).BFonGFPmask  = cellONgfp.PixelIdxList;                    
            segmented_cells(1).BFonGFPori   = cellONgfp.Orientation;
            segmented_cells(1).BFonGFPcen   = cellONgfp.Centroid;
            segmented_cells(1).GFPcentroid  = props(1).Centroid;
            segmented_cells(1).GFPcentroid2 = props(2).Centroid;
            segmented_cells(1).GFPmask      = props(1).PixelIdxList;
            segmented_cells(1).GFPmask2     = props(2).PixelIdxList;
            segmented_cells(1).nucleiD      = nucD;
            segmented_cells(1).GFPstatus    = 'perfect';
            if nucD/segmented_cells.geom.l <0.33
                segmented_cells(1).mitotic = 1;
            end
            tCen = round(props(1).Centroid);
            GFPdot(tCen(1),tCen(2)) = 1;
            tCen = round(props(2).Centroid);
            GFPdot(tCen(1),tCen(2)) = 1;               
        %elseif length(props)>2;    segmented_cells(1).septated= 10;
        end
        if any(any(GFPdot))
            GFPdot = imdilate(GFPdot,strel('disk',4));
            segmented_cells(1).gfpInt.totPerinuc     = sum(GFP_cell(segmented_cells.mask & GFPdot));
            segmented_cells(1).gfpInt.totNonPerinuc  = sum(GFP_cell(segmented_cells.mask & not(GFPdot)));
            segmented_cells(1).gfpInt.meanPerinuc    = mean(GFP_cell(segmented_cells.mask & GFPdot));
            segmented_cells(1).gfpInt.meanNonPerinuc = mean(GFP_cell(segmented_cells.mask & not(GFPdot)));
            segmented_cells(1).gfpInt.medPerinuc     = median(GFP_cell(segmented_cells.mask & GFPdot));
            segmented_cells(1).gfpInt.medNonPerinuc  = median(GFP_cell(segmented_cells.mask & not(GFPdot)));
        end
            
        lastHope= lastHope +1;
        if ~strcmp(segmented_cells(1).GFPstatus,'perfect') && lastHope<2
            midC    = round(cellONgfp.Centroid);
            tMask   = false(size(GFP_cell1));  
            tMask((midC(2)-2):(midC(2)+2),(midC(1)-2):(midC(1)+2)) = 1;
            tMask   = tMask & brightMask;
            GFP_cell1(tMask) = BFthr;
        end
    end
    
    
    if strcmp(segmented_cells(1).GFPstatus,'perfect')
        if do_example_plots
            figure(1);
            cell_boundary = bwboundaries(cell_mask);
            cell_boundary = cell_boundary{1};
            Operations.imshowfit([BF_cell,cell_mask]); hold on;
            plot( cell_boundary(:,2) , cell_boundary(:,1) , 'r');
        end
    elseif sum(cells_interiors(:)) == 0
        segmented_cells(1).GFPstatus = 'thresh_fail';   disp('  cannot get GFP mask');
        if do_example_plots;    figure(2); Operations.imshowfit(BF_cell); hold on;  end
    else
        segmented_cells(1).GFPstatus = 'mask_fail';     disp('  GFP mask not conform');
        if do_example_plots;    figure(3); Operations.imshowfit([BF_cell,cell_mask]); hold on;  end
    end
end

function segmented_cells = redSegmentation(Red_cell,segmented_cells,do_example_plots)
    segmented_cells(1).REDstatus   = Inf;
    thresh     = graythresh(Red_cell);
    RED_thresh = im2bw(Red_cell,thresh);
    cells_interiors = bwareaopen(imfill(imclearborder(RED_thresh),'holes'),10);
    cells_interiors = cells_interiors & segmented_cells.mask;
    if sum(cells_interiors(:)) > 0
        dye_mask   = Operations.getCentralCC(cells_interiors);    
        props = regionprops(dye_mask,'Area','Solidity','Centroid','PixelIdxList');
        if (props.Area < 10000 && props.Solidity > 0.8)
            geom = Operations.mask2geomDye(segmented_cells.mask,dye_mask,60);
            if geom.success
                segmented_cells(1).REDcentroid = props.Centroid;
                segmented_cells(1).REDmask     = dye_mask;
                segmented_cells(1).REDgeom     = geom;
                segmented_cells(1).REDgeom.area= props.Area;
                segmented_cells(1).REDstatus   = 'perfect';
                if do_example_plots
                    figure(1);
                    cell_boundary = bwboundaries(dye_mask);
                    cell_boundary = cell_boundary{1};
                    Operations.imshowfit([Red_cell,dye_mask]); hold on;
                    plot( cell_boundary(:,2) , cell_boundary(:,1) , 'r');
                end
            end
        end
    end
    if segmented_cells.REDstatus == Inf
        disp('BF failure');
        segmented_cells(1).status = Inf;
        if do_example_plots
            figure(3); Operations.imshowfit([Red_cell,dye_mask]); hold on;
        end
    end
end


%% ADDITIONAL FUNCTIONS
function [tMin,tInd] = GFPmicro_aux1(GFP_cell1,tInt,extVar)
    tDel = inf(1,length(tInt));
    for aa=1:length(tInt)
        GFP_thresh  = im2bw(GFP_cell1,tInt(aa));
        cells_interiors = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),50);
        if sum(cells_interiors(:))==0;  continue;   end
        props = bwconncomp(cells_interiors);
        tDel(aa) = min(abs(cellfun(@numel,props.PixelIdxList) - extVar));
        if aa>2 && all(diff(tDel(aa-2:aa))>0); break; end
    end
    [tMin,tInd] = min(tDel);
end

function [tInd,tInt] = GFPmicro_aux2(GFP_cell1,tInt,tD)
    tInd = Inf;
    tAreaVec = inf(1,length(tInt));
    for aa=1:length(tInt)
        GFP_thresh  = im2bw(GFP_cell1,tInt(aa));
        GFP_thresh  = bwareaopen(imfill(imclearborder(GFP_thresh),'holes'),10);
        cells_interiors = imerode(imdilate(GFP_thresh,strel('disk',1)),strel('disk',1));
        if sum(cells_interiors(:))==0;  continue;   end
        props = regionprops(cells_interiors , 'Area','Solidity', 'MajorAxisLength','MinorAxisLength');
        tArea = [props.Area];
        tAreaVec(aa) = sum(tArea);
        if tAreaVec(aa)<100
            tSold = min([props.Solidity]);
            if tSold>0.8
                if      length(props)==1
                    tRatio = props.MajorAxisLength/props.MinorAxisLength;
                    if ~xor((tArea>50),(tRatio>2)); tInd=aa;  break;    end
                   %if      tRatio>2  && tArea>50;  tInd=aa;  break;
                   %elseif  tRatio<=2 && tArea<50;  tInd=aa;  break; 
                elseif  length(props)==2
                    if all(tArea<50);               tInd=aa;  break;    end
                end
            end
        end
    end
    
    if tD~=0
        if      tInd==Inf;  [~,tInd] = min(tAreaVec);   tInt = tInt(max(1,tInd-1)):tD:tInt(min(tInd+1,length(tInt)));
        elseif  tInd==1;    tInt = tInt(1):tD:tInt(2);
        else;               tInt = tInt(tInd-1):tD:tInt(tInd); %#ok<*NOSEM>
        end
    end
end

function segmented_cells = septum_finder(GFP_cell,segmented_cells,outputting) 
    coord = round(segmented_cells.centroid);
    nucMask = false(size(GFP_cell));
    nucMask(segmented_cells.GFPmask)=1;
    nucMask(segmented_cells.GFPmask2)=1;
    nucMask = imdilate(nucMask,[1,1,1]);    %new!
    
    posN =0;  GFPint_middle =zeros(5,1); GFPint_cellOnly =zeros(5,1);
    for i=[-2,-1,0,1,2]
        posN    = posN+1;
        midMask = false(size(GFP_cell));   
        %midMask((coord(2)-2+posN):(coord(2)+2+posN),(coord(1)-2):(coord(1)+2)) = 1;
        midMask((coord(2)-1+posN):(coord(2)+1+posN),(coord(1)-4):(coord(1)+4)) = 1;
        midMask(nucMask)=0;
        GFPint_middle(posN)    = double(median(GFP_cell(midMask)));
        GFPint_cellOnly(posN)  = double(median(GFP_cell(imerode(xor(segmented_cells.mask,or(nucMask,midMask)),[1,1,1]))));
    end
    
    GFPratio = GFPint_middle./GFPint_cellOnly;
    if or(outputting,any(GFPratio>= 1.5))
        segmented_cells(1).septated= max(GFPratio);
    end
end

function segmented_cells = septumMicro_finder(GFP_cell,segmented_cells,outputting)
    brightMask = false(size(GFP_cell));
    brightMask(segmented_cells.BFonGFPmask)=1;
    coord   = round(segmented_cells.BFonGFPcen);
    nucMask = false(size(GFP_cell));
    nucMask(segmented_cells.GFPmask)=1;
    nucMask(segmented_cells.GFPmask2)=1;
    nucMask = imdilate(nucMask,strel('disk',2));    %new!
    
    GFPint_M1 =zeros(9,1); GFPint_M2 = GFPint_M1;
    tInd = [[-1,-1,-1,0,0,0,1,1,1];[-1,0,1,-1,0,1,-1,0,1]];
    for i=1:9
        midMask1 = false(size(GFP_cell));   midMask2 = midMask1;
        midMask1((coord(2)-1+tInd(1,i)):(coord(2)+1+tInd(1,i)),(coord(1)-1+tInd(2,i)):(coord(1)+1+tInd(2,i))) = 1;   %midMask(nucMask) =0;
        midMask2((coord(2)-2+tInd(1,i)):(coord(2)+2+tInd(1,i)),(coord(1)-2+tInd(2,i)):(coord(1)+2+tInd(2,i))) = 1;   %midMask(nucMask) =0;
        GFPint_M1(i) = double(median(GFP_cell(midMask1)));
        GFPint_M2(i) = double(median(GFP_cell(midMask2)));                
    end
    outMask     = xor(brightMask,nucMask);
    outMask((coord(2)-4):(coord(2)+4),(coord(1)-4):(coord(1)+4)) = 0;
    GFPint_O    = double(median(GFP_cell(outMask)));
    
    GFPratio = [GFPint_M1/GFPint_O; GFPint_M2/GFPint_O];
    GFPdelta = [GFPint_M1; GFPint_M2] - (GFPint_O + std((GFP_cell(outMask))));
    %if or(outputting,any(GFPratio>= 1.5)); segmented_cells(1).septated= max(GFPratio);     end
    if or(outputting,any(GFPdelta>0));      segmented_cells(1).septated= max(GFPratio);     end
end

function segmented_cells = cut3_finder(Red_cell,segmented_cells,outputting) 
    nucMask = false(size(Red_cell));
    nucMask(segmented_cells.GFPmask)=1;
    if segmented_cells.nucleiN==2
        nucMask(segmented_cells.GFPmask2)=1;
    end
    nucMaskD = imdilate(nucMask,[1,1,1]);    %new!

    cut3_nuc    = double(median(Red_cell(nucMask)));
    cut3_nucD   = double(median(Red_cell(nucMaskD)));
    cut3_noNuc  = double(median(Red_cell(imerode(xor(segmented_cells.mask,nucMask),[1,1,1]))));
    cut3_noNucD = double(median(Red_cell(imerode(xor(segmented_cells.mask,nucMaskD),[1,1,1]))));
    
    cut3_ratio  = cut3_nuc/cut3_noNuc;
    cut3_ratioD = cut3_nucD/cut3_noNucD;
    if or(outputting,or(cut3_ratio>= 1.5, cut3_ratioD>=1.5))
        segmented_cells(1).mitotic = max([cut3_ratio, cut3_ratioD]);
    end
end

function segmented_cells = cut3Micro_finder(Red_cell,segmented_cells,outputting)
    extVar = regionprops(imerode(segmented_cells.mask,strel('disk',3)),'Area');
    [~,tInd] = max([extVar.Area]); extVar = extVar(tInd);
    if isempty(extVar)
        extVar = regionprops(segmented_cells.mask,'Area');
        [~,tInd] = max([extVar.Area]); extVar = extVar(tInd);
         if isempty(extVar)
             segmented_cells(1).CUT3status = 'off target';
             return
         end
    end
    for ii=[30 25 35 20 40 15 45 10 50 5 55]
        Red_thresh = adaptivethresholdSH(Red_cell,ii,0);
        cells_interiors = bwareaopen(imfill(imclearborder(~Red_thresh),'holes'),100);
        cells_interiors = imdilate(imerode(cells_interiors,strel('disk',2)),strel('disk',1));
        if sum(cells_interiors(:)) > 0
            props = regionprops(cells_interiors,'Area','Centroid','PixelIdxList','Orientation'); tArea=[props.Area];
            if any(tArea>0.75*extVar.Area) && any(tArea<1.25*extVar.Area);    break;  end
        end
    end
    if sum(cells_interiors(:)) == 0
        segmented_cells(1).CUT3status = 'off target';
        return
    end
    tVar = abs([props.Area]-segmented_cells.area);
    [~,tDelI]    = min(tVar);     tDelR = diff(abs(tVar))/max(tVar);    tVar=[];
    for ii=1:length(props)
        c1 = props(ii).Centroid./size(Red_cell);    c2 = segmented_cells.centroid./size(segmented_cells.mask);
        tVar(ii) = sqrt(sum((c1-c2).^2));
    end
    [~,tDisI] = min(tVar);         tDisR = diff(abs(tVar))/max(tVar);   %tVar=[];
    if tDelR>tDisR; tInd = tDelI;   else;    tInd = tDisI;   end
    brightMask  = false(size(cells_interiors));
    brightMask(props(tInd).PixelIdxList) = 1;

    
    %segmented_cells.BFonGFPori         %cell mask orientation on GFP    
    tOri = props(tInd).Orientation;     %cell mask orientation on Red
    if segmented_cells.BFonGFPori<0;    segmented_cells.BFonGFPori=segmented_cells.BFonGFPori+180;  end
    if tOri<0;                          tOri = tOri+180;    end
    %segmented_cells.BFonGFPcen         %cell mask centroid on GFP
    tDis = props(tInd).Centroid;        %cell mask centroid on Red
    %segmented_cells.GFPcentroid        %GFP  mask centroid on GFP   
    %GFP  mask centroid on Red: translation + rotation
    tCen  = segmented_cells.GFPcentroid + (tDis - segmented_cells.BFonGFPcen);
    alpha = pi*(tOri-segmented_cells.BFonGFPori )/180;
    tDel    = tCen - tDis;
    tDel    = [tDel(1)*cos(alpha) + tDel(2)*sin(alpha), -tDel(1)*sin(alpha) + tDel(2)*cos(alpha)];
    tCen    = round(tDel + tDis);
    nucMask = false(size(Red_cell));    nucMask(tCen(2),tCen(1))=1;
    tDilate = floor(sqrt(length(segmented_cells.GFPmask)/pi))+1;
    nucMask = imdilate(nucMask,strel('disk',tDilate));
    
    if segmented_cells.nucleiN==2
        tCen2 = segmented_cells.GFPcentroid2 + (tDis - segmented_cells.BFonGFPcen);
        alpha = pi*(tOri-tInd)/180;
        tDel    = tCen2 - tDis;
        tDel    = [tDel(1)*cos(alpha) + tDel(2)*sin(alpha), -tDel(1)*sin(alpha) + tDel(2)*cos(alpha)];
        tCen2    = round(tDel + tDis);
        nucMask2= false(size(Red_cell));    nucMask2(tCen2(2),tCen2(1))=1;
        tDilate = floor(sqrt(sum(segmented_cells.GFPmask2(:))/pi))+1;
        nucMask2 = imdilate(nucMask2,strel('disk',tDilate));
        nucMask = nucMask | nucMask2;
    end    
    nucMaskD = imdilate(nucMask,[1,1,1]);    %new!

    
    cut3_nuc    = double(median(Red_cell(nucMask)));
    cut3_nucD   = double(median(Red_cell(nucMaskD)));
    cut3_noNuc  = double(median(Red_cell(imerode(xor(brightMask,nucMask),[1,1,1]))));
    cut3_noNucD = double(median(Red_cell(imerode(xor(brightMask,nucMaskD),[1,1,1]))));    
    cut3_ratio  = cut3_nuc/cut3_noNuc;
    cut3_ratioD = cut3_nucD/cut3_noNucD;
    if or(outputting,or(cut3_ratio>= 1.5, cut3_ratioD>=1.5))
        segmented_cells(1).mitotic = max([cut3_ratio, cut3_ratioD]);
    end
    segmented_cells(1).CUT3status='perfect';
end

function segmented_cells = ratio_finder(segmented_cells)
    segmented_cells(1).REDgeom.ratio_area  = segmented_cells(1).REDarea/segmented_cells(1).area;
   %segmented_cells(1).REDgeom.ratio_width = segmented_cells(1).REDgeom.dye_width_avg/segmented_cells(1).REDgeom.bf_width_avg;
    segmented_cells(1).REDgeom.ratio_length= segmented_cells(1).REDgeom.l/segmented_cells(1).geom.l;
    segmented_cells(1).REDgeom.deltaYC     = abs(segmented_cells.REDcentroid(2) - segmented_cells.centroid(2));
end
