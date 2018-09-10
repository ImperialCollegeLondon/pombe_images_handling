function ExploNOclickPoints(fov,output3,clicks)
        if nargin<3; clicking=false; else; clicking=true; end
            
        BF_cell  = fov.loadImage(1,1,'BF', 'normalize_contrast');
        points = struct('x',{},'y',{}); counter = 0;    masks = struct('m',{},'Box',{},'im',{});
        Cs = [0,0.01,0.02,0.03,0.04];   
        for atts = 1:5
            pro2 =  bwareaopen(imfill(imclearborder(adaptivethresholdSH(BF_cell,10,Cs(atts))),'holes'),100);
            props = regionprops(pro2,'Area','Solidity','Centroid','PixelIdxList','BoundingBox','Image');
            for qq=1:length(props)
                if (props(qq).Area > 1000 && props(qq).Area < 10000 && props(qq).Solidity > 0.8)
                    counter = counter+1;
                    points(counter).x   = props(qq).Centroid(1);
                    points(counter).y   = props(qq).Centroid(2);
                    masks(counter).m    = props(qq).PixelIdxList;
                    masks(counter).Box  = props(qq).BoundingBox;
                end
            end
        end

        tFind = zeros(length(masks),1);
        for mm = 1:(length(masks)-1)
            t1 = points(mm);
            for nn = (mm+1):length(masks)
                t2 = masks(nn).Box;
                if  t1.y>t2(2) && t1.y<(t2(2)+t2(4)) && t1.x>t2(1) && t1.x<(t2(1)+t2(3))
                    tFind(nn) = 1;
                end 
            end
        end
        points = points(tFind==0);  masks = masks(tFind==0);        
        
        if clicking
            cFind = zeros(length(masks),1);
            for mm = 1:length(masks)
                t2 = masks(mm).Box;
                for nn = 1:length(clicks)
                    t1 = clicks(nn);
                    if  t1.y>t2(2) && t1.y<(t2(2)+t2(4)) && t1.x>t2(1) && t1.x<(t2(1)+t2(3))
                        cFind(mm) = 1;
                        continue
                    end 
                end
            end 
            points = points(cFind==1);  masks = masks(cFind==1);
        end
        save(output3,'points','fov','masks');     
end