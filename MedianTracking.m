%% Setup & Crop %%

inf = dicomCollection('S3_tag614.dcm');

framecount = inf{1, 9};
intr = 1;

rect = [243 158 640 510];
flip = 1;
% crop the window you want to analyze (optional)
% cropme = dicomread('M1_tag643.dcm', 'frames', intr);
% [cropped, rect] = imcrop(cropme);

numsegs = 100;
filename = 'results.gif';

%% Select Cycles %%
montage = zeros(511, framecount);
montage_col = rect(3)/2;
for t = 1:300
    cur = imcrop(dicomread('M1_tag643', 'frames', t), rect);
    col = cur(:, montage_col);
    montage(:, t) = col; 
end
imtool(uint8(montage));
pause;

%% Core Loop %%

fstart = 26;
fend = 27;
for f = fstart : +intr : fend    
    rwimg = imcrop(dicomread('S3_tag614.dcm', 'frames', f), rect);
    img = histeq(rgb2gray(rwimg));
    if(flip == 1)
        img = fliplr(img);
    end
    disp(f);

    %% Segment Walls %%  
    
    bin_epi = img > 160;
    bin_endo = img < 40;
    
    epi = imfill(bwareaopen(bin_epi, 10000), 'holes');
    endo = imfill(bwareaopen(bin_endo, 10000), 'holes');
    
    %% Detect Edges %%
    
    epiedge = edge(epi, 'canny');
    endoedge = edge(endo, 'canny');
    edges = imoverlay(epiedge, endoedge, [1 0 0]);
    
    %% Epicardium Curve %%
    
    xep = 360;
    yep = 250;
    plot(xep, yep, 'c.');
    
    st_epi = 250;
    end_epi = -30;
    ep_numpoints = 64;
    inc = (end_epi-st_epi)/ep_numpoints;
    
    epipoints = zeros(2, 1);
    count = 1;
    endoct = 1;
    for t = st_epi : inc : end_epi
        r = 0;
        searching = true;
        while(searching) 
            r = r + 1;
            x = round(r * cosd(t)) + xep;
            y = round(r * sind(t)) + yep;
        
            if(x > rect(4)-1 || x < 2) 
                break; end
            if(y > rect(3)-1 || y < 2)
                break; end
            for h = x-1 : +1 : x+1
                for k = y-1 : +1 : y+1
                    if((epiedge(h, k) > 0) && searching)
                        epipoints(1, count) = y;
                        epipoints(2, count) = x;
                        count = count + 1;
                        searching = false;
                    end
                end
            end
        end
    end
    imshow(edges); hold on,
    
    epicrv = cscvn(epipoints);
    fnplt(epicrv);
    plot(epipoints(1, : ), epipoints(2, : ), 'o');
    
    %% Endocardium Curve %%
    
    xen = xep;
    yen = yep;
    
    st_endo = 210;
    end_endo = -30;
    en_numpoints = 64;
    inc = (end_endo-st_endo)/en_numpoints;
    
    endopoints = zeros(2, 0);
    count = 1;
    for t = st_endo : inc : end_endo
        r = 0;
        maxr = r;
        searching = true;
        while(searching)
            r = r + 1;
            x = round(r * cosd(t)) + xen;
            y = round(r * sind(t)) + yen;
        
            if(x > rect(4)-1 || x < 2) 
                break; end
            if(y > rect(3)-1 || y < 2)
                break; end
            for h = x-1 : +1 : x+1
                for k = y-1 : +1 : y+1
                    if((endoedge(h, k) > 0) && searching)
                        maxr = r;
                    end
                end
            end
        end
        endopoints(1, count) = round(maxr * sind(t)) + yen;
        endopoints(2, count) = round(maxr * cosd(t)) + xen;
        count = count + 1;
    end
    
    endocrv = cscvn(endopoints);
    fnplt(endocrv);
    plot(endopoints(1, : ), endopoints(2, : ), 'o')
    
    %% Generate Ellipses %%
     
    % find ellipses
    options.isplot = false;
    options.verbose = false;
    try 
        epiellipse = ellipseFit4HC(epipoints(1, : ), epipoints(2, : ), options);
        endoellipse = ellipseFit4HC(endopoints(1, : ), endopoints(2, : ), options);
    catch ERR
       
    end
    
    % plot ellipses
    t = linspace(0, 2*pi, numsegs)';
    plot(endoellipse.X_fitfun(t), endoellipse.Y_fitfun(t),'--');
    plot(epiellipse.X_fitfun(t), epiellipse.Y_fitfun(t),'--');
%     
%     plot(endoellipse.mu_X_fit, endoellipse.nu_Y_fit, 'h');
%     plot(epiellipse.mu_X_fit, epiellipse.nu_Y_fit, 'h');
    
    %% Adjust Lines Using Ellipses %%
    
    % find avg dist. between estimation and detected points
    endo_avgdist = 0;
    for i=1 : 1 : size(endopoints, 2)
        dist = sqrt((endopoints(1, i)-endoellipse.mu_X_fit(i, 1))^2 + (endopoints(2, i)-endoellipse.nu_Y_fit(i, 1))^2);
        endo_avgdist = endo_avgdist+dist;
    end
    endo_avgdist = endo_avgdist/size(endopoints, 2);
    epi_avgdist = 0;
    for i=1 : 1 : size(epipoints, 2)
        dist = sqrt((epipoints(1, i)-epiellipse.mu_X_fit(i, 1))^2 + (epipoints(2, i)-epiellipse.nu_Y_fit(i, 1))^2);
        epi_avgdist = epi_avgdist+dist;
    end
    epi_avgdist = epi_avgdist/size(epipoints, 2);
    
    % compute standard dev.
    epi_std = 0;
    for i=1 : 1 : size(epipoints, 2)
        dist = sqrt((epipoints(1, i)-epiellipse.mu_X_fit(i, 1))^2 + (epipoints(2, i)-epiellipse.nu_Y_fit(i, 1))^2);
        epi_std = epi_std + (dist-epi_avgdist)^2;
    end
    epi_std = sqrt(epi_std/(size(epipoints, 2)-1));
    endo_std = 0;
    for i=1 : 1 : size(epipoints, 2)
        dist = sqrt((endopoints(1, i)-endoellipse.mu_X_fit(i, 1))^2 + (endopoints(2, i)-endoellipse.nu_Y_fit(i, 1))^2);
        endo_std = endo_std + (dist-endo_avgdist)^2;
    end
    endo_std = sqrt(endo_std/(size(endopoints, 2)-1));
    
    % adjust points
    endopoints_a = zeros(size(endopoints));
    for i=1 : 1 : size(endopoints, 2)
        dist = sqrt((endopoints(1, i)-endoellipse.mu_X_fit(i, 1))^2 + (endopoints(2, i)-endoellipse.nu_Y_fit(i, 1))^2);
        if abs(dist-endo_avgdist) <= endo_std
            endopoints_a(1, i) = (endopoints(1, i)+2*endoellipse.mu_X_fit(i, 1))/3;
            endopoints_a(2, i) = (endopoints(2, i)+2*endoellipse.nu_Y_fit(i, 1))/3;
        else
            endopoints_a(1, i) = endoellipse.mu_X_fit(i, 1);
            endopoints_a(2, i) = endoellipse.nu_Y_fit(i, 1);
        end        
    end
    
    epipoints_a = zeros(size(epipoints));
    for i=1 : 1 : size(epipoints, 2)
        dist = sqrt((epipoints(1, i)-epiellipse.mu_X_fit(i, 1))^2 + (epipoints(2, i)-epiellipse.nu_Y_fit(i, 1))^2);
        if abs(dist-epi_avgdist) <= epi_std
            epipoints_a(1, i) = (epipoints(1, i)+2*epiellipse.mu_X_fit(i, 1))/3;
            epipoints_a(2, i) = (epipoints(2, i)+2*epiellipse.nu_Y_fit(i, 1))/3;
        else
            epipoints_a(1, i) = epiellipse.mu_X_fit(i, 1);
            epipoints_a(2, i) = epiellipse.nu_Y_fit(i, 1);
        end
    end
    
    % create new walls and interpolate
    endoadj = cscvn(endopoints_a);
    %fnplt(endoadj);
    [endo_x, endo_y] = interpolate_cscvn(endoadj);
    
    epiadj = cscvn(epipoints_a);
    %fnplt(epiadj);
    [epi_x, epi_y] = interpolate_cscvn(epiadj);
    
    %plot(epi_x, epi_y, 'o');
    %plot(endo_x, endo_y, 'o');
    
    %% Evenly Distribute Points %%
    
    % epicardium
    max = 0;
    for i=1 : 1 : size(epi_x, 1)-1
       dist = sqrt((epi_x(i)-epi_x(i+1))^2 + (epi_y(i)-epi_y(i+1))^2);
       if dist > max
           max = dist;
       end
    end
    
    epi_spaced = zeros(2, 1);
    epi_spaced(1, 1) = epi_x(1);
    epi_spaced(2, 1) = epi_y(1);
    hasNext = true;
    ind = 1;
    sp_ind = 1;
    while hasNext
        searching = true;
        j = 1;
        while searching
            try
                dist = sqrt((epi_x(ind)-epi_x(ind+j))^2 + (epi_y(ind)-epi_y(ind+j))^2);
                if dist >= max
                    epi_spaced(1, sp_ind+1) = epi_x(ind+j);
                    epi_spaced(2, sp_ind+1) = epi_y(ind+j);
                    searching = false;
                    ind = ind+j;
                    sp_ind = sp_ind+1;
                end
                j = j+1;
            catch ME
                hasNext = false;
                break;
            end
        end
        if ind >= size(epi_x, 1)
            hasNext = false;
        end
    end
    
    plot(epi_spaced(1, : ), epi_spaced(2, : ), 'o')
    
    % endocardium
    max = 0;
    for i=1 : 1 : size(endo_x, 1)-1
       dist = sqrt((endo_x(i)-endo_x(i+1))^2 + (endo_y(i)-endo_y(i+1))^2);
       if dist > max
           max = dist;
       end
    end
    
    endo_spaced = zeros(2, 1);
    endo_spaced(1, 1) = endo_x(1);
    endo_spaced(2, 1) = endo_y(1);
    hasNext = true;
    ind = 1;
    sp_ind = 1;
    while hasNext
        searching = true;
        j = 1;
        while searching
            try
                dist = sqrt((endo_x(ind)-endo_x(ind+j))^2 + (endo_y(ind)-endo_y(ind+j))^2);
                if dist >= max
                    endo_spaced(1, sp_ind+1) = endo_x(ind+j);
                    endo_spaced(2, sp_ind+1) = endo_y(ind+j);
                    searching = false;
                    ind = ind+j;
                    sp_ind = sp_ind+1;
                end
                j = j+1;
            catch ME
                hasNext = false;
                break;
            end
        end
        if ind >= size(endo_x, 1)
            hasNext = false;
        end
    end
    
    %plot(endo_spaced(1, : ), endo_spaced(2, : ), 'o')
    
    % update interpolation
    endoadj = cscvn(endo_spaced);
    fnplt(endoadj);
    [endo_x, endo_y] = interpolate_cscvn(endoadj);
    plot(endo_x, endo_y, 'o');
    
    epiadj = cscvn(epi_spaced);
    fnplt(epiadj);
    [epi_x, epi_y] = interpolate_cscvn(epiadj);
    plot(epi_x, epi_y, 'o');
    
    %% Generate ROIs %%
    
    if f==fstart
        % adjustable #ROIs (+1)
        numROIs = 23;

        % height and width of ROIs
        hROI = 50;
        wROI = 70;
        
        %center
        h = 270;
        k = 270;
        plot(h, k, 'c.')
    
        endoscale = 1;
        % the endocardial wall is much smaller, and therefore there may be
        % too much overlap between regions. use this number to reduce its
        % size if needed
        
        % now we need to make the ROIs, evenly distributed around the walls
        endo_used = zeros(length(epi_x(:)), 1); %track whether or not we've already pointed to a point for an ROI on the epicardium
        
        epi_ROIs_x = zeros(1, numROIs+1);
        epi_ROIs_y = zeros(1, numROIs+1);

        endo_ROIs_x = zeros(1, numROIs+1);
        endo_ROIs_y = zeros(1, numROIs+1);

        % generating ellipses will allow us to pair regions better in a
        % tangential manner
        epiellipse = ellipseFit4HC(epi_x(:), epi_y(:), options);
%         endoellipse = ellipseFit4HC(endo_x(:), endo_y(:), options);
        
        t = cumsum([0; diff(epiellipse.mu_X_fit(:)).^2 + diff(epiellipse.nu_Y_fit(:)).^2]);
        
        splx = spline(t,epiellipse.mu_X_fit);
        sply = spline(t,epiellipse.nu_Y_fit);
        
        epi_tans = ppval(fnder(sply),t)./ppval(fnder(splx),t);
        
%         t = cumsum([0; diff(endoellipse.mu_X_fit(:)).^2 + diff(endoellipse.nu_Y_fit(:)).^2]);
%         
%         splx = spline(t,endoellipse.mu_X_fit);
%         sply = spline(t,endoellipse.nu_Y_fit);
%         
%         endo_tans = ppval(fnder(sply),t)./ppval(fnder(splx),t);
        
        count = 1;
        
        endo_grid = zeros(rect(4), rect(3));
        for j = 1 : length(endo_x(:))
            endo_grid(round(endo_y(j)), round(endo_x(j))) = 1;
        end
        
        for j = floor(size(epi_x, 1)/(numROIs+1)) : floor(size(epi_x, 1)/(numROIs+1)) : size(epi_x, 1)

            perp = -(epi_tans(j))^(-1);
            d = atand(perp);

            incR = 0.5;
            if (epi_y(j) < k && d < 0) || (epi_y(j) > k && d > 0)
                incR = -incR;
            end
%             beste = 1000000; %just something really high, essentially an infinite bound
%             bind = 0; %beste index
%             bdist = 1000000; % best distance
            
            distLIM = sqrt((h-epi_x(j))^2 + (k-epi_y(j))^2);

            r = 0;
            minr = 100000;
            minDist = distLIM;
            searching = true;
            while(searching)
                r = r + incR;
                x = round(r * cosd(d) + epi_x(j));
                y = round(r * sind(d) + epi_y(j));

                if(x > rect(3)-2 || x < 2) 
                    break; end
                if(y > rect(4)-2 || y < 2)
                    break; end
                for b = x-1 : x+1
                    for n = y-1 : y+1
                        dist = sqrt((b-epi_x(j))^2 + (n-epi_y(j))^2);
                        if((endo_grid(n, b) == 1) && searching) && dist < minDist
                            minr = r; 
                            minDist = dist;
                        end
                    end
                end
            end
            

            
%             for b = 1 : length(endo_x(:))
%                 dist = sqrt((endo_x(b)-epi_x(j))^2 + (endo_y(b)-epi_y(j))^2);
%                 slp = (endo_y(b)-epi_y(j))/(endo_x(b)-epi_x(j));
%                 
%                 if abs(perp-slp)<beste && endo_used(b) == 0 && dist < distLIM
%                     beste = abs(perp-slp);
%                     bind = b;
%                 end
%             end
            
            x = minr * cosd(d) + epi_x(j);
            y = minr * sind(d) + epi_y(j);
            
            % sometimes it doesn't catch the lip, so just tack it on the
            % end            
            if minDist == distLIM
                if epi_y(j) < k
                    x = endo_x(1);
                    y = endo_y(1);
                else
                    x = endo_x(length(endo_x(:))); 
                    y = endo_y(length(endo_y(:))); 
                end
            end

            line([epi_x(j) x], [epi_y(j) y])
%             endo_used(bind) = 1;

            endo_ROIs_x(count) = x;
            endo_ROIs_y(count) = y;
            epi_ROIs_x(count) = epi_x(j);
            epi_ROIs_y(count) = epi_y(j);

            endo_ROI = [endo_ROIs_x(count)-wROI*0.5*endoscale, endo_ROIs_y(count)-hROI*0.5*endoscale, wROI*endoscale, hROI*endoscale];
            epi_ROI = [epi_ROIs_x(count)-wROI*0.5, epi_ROIs_y(count)-hROI*0.5, wROI, hROI];

            rectangle('Position', endo_ROI, 'EdgeColor', 'b');
            rectangle('Position', epi_ROI, 'EdgeColor', 'b');

            count = count + 1;
        end
        
%         for j = floor(size(endo_x, 1)/(numROIs+1)) : floor(size(endo_x, 1)/(numROIs+1)) : size(endo_x, 1)            
% 
%             dx = endo_x(j) - h;
%             dy = endo_y(j) - k;
%             dist = sqrt(dx*dx + dy*dy);
%             ux = dx/dist; %unit vector components
%             uy = dy/dist;
%             theta = atan(dy/dx);
% 
%             beste = 1000000; %just something really high, essentially an infinite bound
%             bind = 0; %beste index
% 
%             for b = 1 : length(epi_x(:))
%                 dx = epi_x(b) - h;
%                 dy = epi_y(b) - k;   
%                 if (dx/ux > 0) && (dy/uy > 0) && epi_used(b)==0 && abs(dx/ux - dy/uy)<beste
%                     beste = abs(dx/ux - dy/uy);
%                     bind = b;
%                 end
%             end
%             if beste <= 30
%                 line([epi_x(bind) endo_x(j)], [epi_y(bind) endo_y(j)])
%                 epi_used(bind) = 1;
%                 
%                 endo_ROIs_x(count) = endo_x(j);
%                 endo_ROIs_y(count) = endo_y(j);
%                 epi_ROIs_x(count) = epi_x(bind);
%                 epi_ROIs_y(count) = epi_y(bind);
%                 
%                 endo_ROI = [endo_ROIs_x(count)-wROI*0.5*endoscale, endo_ROIs_y(count)-hROI*0.5*endoscale, wROI*endoscale, hROI*endoscale];
%                 epi_ROI = [epi_ROIs_x(count)-wROI*0.5, epi_ROIs_y(count)-hROI*0.5, wROI, hROI];
%                 
%                 rectangle('Position', endo_ROI, 'EdgeColor', 'b');
%                 rectangle('Position', epi_ROI, 'EdgeColor', 'b');
%                 
%                 count = count + 1;
%             end
%         end
        
        % create copies for strain calc
        
        endo_ROIs_x_OG = endo_ROIs_x;
        endo_ROIs_y_OG = endo_ROIs_y;
        
        epi_ROIs_x_OG = epi_ROIs_x;
        epi_ROIs_y_OG = epi_ROIs_y;
        
    end
    %% Trackers %%
    
    image1=rwimg;
    margin = 40;
    
    if f==fstart
        points = detectMinEigenFeatures(histeq(rgb2gray(image1)));
        pointsLoc1=points.Location;
%         plot(pointsLoc1(:,1),pointsLoc1(:,2),'r+')
        
        % we need to store the indexes of the trackers in each ROI, so that
        % they can be revisited later
        
        endo_ROIinds = zeros (numROIs+1, 1);
        epi_ROIinds = zeros (numROIs+1, 1);
        
%         new_pointsLoc1 = zeros(1,2);
%         
%         epc = 1;
%         for j = 1 : length(epi_x)  
%             for k = 1 : length(pointsLoc1)
%                 if abs(epi_x(j) - pointsLoc1(k, 1))<=margin && abs(epi_y(j) - pointsLoc1(k, 2))<=margin
%                     if epc > 1
%                         if pointsLoc1(k, :) == new_pointsLoc1(epc-1, :)
%                             break;
%                         end
%                     end
%                     plot(pointsLoc1(k,1),pointsLoc1(k,2),'y+')
%                     new_pointsLoc1(epc, :) = pointsLoc1(k, :);
%                     epc = epc+1;
% %                     break;
%                 end
%             end
%         end
% 
% %       find trackers around endocardium
% %       also puts points in radial order
%         enc = epc;
%         for j = 1 : length(endo_x)
%             for k = 1 : length(pointsLoc1)
%                 if abs(endo_x(j) - pointsLoc1(k, 1))<=margin && abs(endo_y(j) - pointsLoc1(k, 2))<=margin
%                     if enc > 1
%                         if pointsLoc1(k, :) == new_pointsLoc1(enc-1, :)
%                             break;
%                         end
%                     end
%                     plot(pointsLoc1(k,1),pointsLoc1(k,2),'y+')
%                     new_pointsLoc1(enc, :) = pointsLoc1(k, :);
%                     enc = enc+1;
% %                     break;
%                 end
%             end
%         end
%         
%         pointsLoc1 = new_pointsLoc1;
        
        for i=1:numROIs+1
            ep_trind = 1;
            en_trind = 1;
            for j=1:length(pointsLoc1)
                if pointsLoc1(j,1) >= endo_ROIs_x(i)-wROI*0.5*endoscale && pointsLoc1(j,1) <= endo_ROIs_x(i)+wROI*0.5*endoscale && pointsLoc1(j,2) >= endo_ROIs_y(i)-hROI*0.5*endoscale && pointsLoc1(j,2) <= endo_ROIs_y(i)+hROI*0.5*endoscale
                    endo_ROIinds(i, en_trind) = j;
                    plot(pointsLoc1(j,1),pointsLoc1(j,2),'g+')
                    en_trind = en_trind + 1;
                end
                
                if pointsLoc1(j,1) >= epi_ROIs_x(i)-wROI*0.5 && pointsLoc1(j,1) <= epi_ROIs_x(i)+wROI*0.5 && pointsLoc1(j,2) >= epi_ROIs_y(i)-hROI*0.5 && pointsLoc1(j,2) <= epi_ROIs_y(i)+hROI*0.5
                    epi_ROIinds(i, ep_trind) = j;
                    plot(pointsLoc1(j,1),pointsLoc1(j,2),'g+')
                    ep_trind = ep_trind + 1;
                end
            end
            line([endo_ROIs_x(i) epi_ROIs_x(i)],[endo_ROIs_y(i) epi_ROIs_y(i)])
        end
        
        %init trackers
        tracker = vision.PointTracker('MaxBidirectionalError',1,'NumPyramidLevels', 5);
        initialize(tracker, pointsLoc1, image1);
        
        % this is where we will store the median x and y displacement
        epi_meds_x = zeros(numROIs+1, fend-fstart+1);
        epi_meds_y = zeros(numROIs+1, fend-fstart+1);
        
        endo_meds_x = zeros(numROIs+1, fend-fstart+1);
        endo_meds_y = zeros(numROIs+1, fend-fstart+1);
        
    else
        [pointsLoc, validity] = tracker(image1);
        for j=1:length(pointsLoc1)
            if validity(j)
                plot(pointsLoc1(j,1),pointsLoc1(j,2),'r+')
                plot(pointsLoc(j,1),pointsLoc(j,2),'g+')
                line([pointsLoc(j,1) pointsLoc1(j,1)],[pointsLoc(j,2) pointsLoc1(j,2)])
            end
        end
        
        %record median displacements in x and y
        for j=1:numROIs+1
            
            xvals = zeros(1,1);
            yvals = zeros(1,1);
            
            for i=1:length(endo_ROIinds(j, :))
                if endo_ROIinds(j, i) == 0
                    break;
                end
                if validity(endo_ROIinds(j, i))
                    xvals(1, i) = pointsLoc(endo_ROIinds(j, i), 1) - pointsLoc1(endo_ROIinds(j, i), 1);
                    yvals(1, i) = pointsLoc(endo_ROIinds(j, i), 2) - pointsLoc1(endo_ROIinds(j, i), 2);
                    plot(pointsLoc(endo_ROIinds(j, i),1),pointsLoc(endo_ROIinds(j, i),2),'g+')
                end
            end
            
            % we can use unique here to get rid of 0s, since it is highly
            % improbable that the displacements will ever be exactly the
            % same to 3 decimals
            endo_meds_x(j, f) = mean(unique(xvals), 'all');
            endo_meds_y(j, f) = mean(unique(yvals), 'all');
            
            xvals = zeros(1,1);
            yvals = zeros(1,1);
            
            for i=1:length(epi_ROIinds(j, :))
                if epi_ROIinds(j, i) == 0
                    break;
                end
                if validity(epi_ROIinds(j, i))
                    xvals(1, i) = pointsLoc(epi_ROIinds(j, i), 1) - pointsLoc1(epi_ROIinds(j, i), 1);
                    yvals(1, i) = pointsLoc(epi_ROIinds(j, i), 2) - pointsLoc1(epi_ROIinds(j, i), 2);
                    plot(pointsLoc(epi_ROIinds(j, i),1),pointsLoc(epi_ROIinds(j, i),2),'g+')
                end
            end
            

            epi_meds_x(j, f) = median(unique(xvals), 'all');
            epi_meds_y(j, f) = median(unique(yvals), 'all');
            
            %draw ROIs
            endo_ROI = [endo_ROIs_x_OG(j)-wROI*0.5*endoscale, endo_ROIs_y_OG(j)-hROI*0.5*endoscale, wROI*endoscale, hROI*endoscale];
            epi_ROI = [epi_ROIs_x_OG(j)-wROI*0.5, epi_ROIs_y_OG(j)-hROI*0.5, wROI, hROI];
            rectangle('Position', endo_ROI, 'EdgeColor', 'b');
            rectangle('Position', epi_ROI, 'EdgeColor', 'b');
            
            % now we shift each ROI by the median displacement of its
            % trackers, then reacquire new trackers in its new spot (most
            % of which will be the same ones)

            endo_ROIs_x(j) = endo_ROIs_x(j) + endo_meds_x(j, f);
            endo_ROIs_y(j) = endo_ROIs_y(j) + endo_meds_y(j, f);
             
            epi_ROIs_x(j) = epi_ROIs_x(j) + epi_meds_x(j, f);
            epi_ROIs_y(j) = epi_ROIs_y(j) + epi_meds_y(j, f);
            
            %redraw
            endo_ROI = [endo_ROIs_x(j)-wROI*0.5*endoscale, endo_ROIs_y(j)-hROI*0.5*endoscale, wROI*endoscale, hROI*endoscale];
            epi_ROI = [epi_ROIs_x(j)-wROI*0.5, epi_ROIs_y(j)-hROI*0.5, wROI, hROI];
            rectangle('Position', endo_ROI, 'EdgeColor', 'g');
            rectangle('Position', epi_ROI, 'EdgeColor', 'g');          
        end
        
        endo_ROIinds = zeros (numROIs+1, 1);
        epi_ROIinds = zeros (numROIs+1, 1);
        
        % re-acquire trackers
        for i=1:numROIs+1
            ep_trind = 1;
            en_trind = 1;
            for j=1:length(pointsLoc1)
                if pointsLoc1(j,1) >= endo_ROIs_x(i)-wROI*0.5*endoscale && pointsLoc1(j,1) <= endo_ROIs_x(i)+wROI*0.5*endoscale && pointsLoc1(j,2) >= endo_ROIs_y(i)-hROI*0.5*endoscale && pointsLoc1(j,2) <= endo_ROIs_y(i)+hROI*0.5*endoscale
                    endo_ROIinds(i, en_trind) = j;
                     plot(pointsLoc1(j,1),pointsLoc1(j,2),'b+')
                    en_trind = en_trind + 1;
                end
                
                if pointsLoc1(j,1) >= epi_ROIs_x(i)-wROI*0.5 && pointsLoc1(j,1) <= epi_ROIs_x(i)+wROI*0.5 && pointsLoc1(j,2) >= epi_ROIs_y(i)-hROI*0.5 && pointsLoc1(j,2) <= epi_ROIs_y(i)+hROI*0.5
                    epi_ROIinds(i, ep_trind) = j;
                    plot(pointsLoc1(j,1),pointsLoc1(j,2),'b+')
                    ep_trind = ep_trind + 1;
                end
            end
        end
        
        tracker.release;
        points = detectMinEigenFeatures(rgb2gray(image1));
        pointsLoc1 = pointsLoc; %original location
        tracker = vision.PointTracker('MaxBidirectionalError',1,'NumPyramidLevels', 5); %initialize a tracker
        initialize(tracker, abs(pointsLoc1), image1);
    end
    %% Calculate Strain %%
    
    %store strain values in matrix
    %each row is the strain displacement for a region, columns are frame
        
    if f == fstart
        num_elts = fend-fstart+1;
        strain = zeros(numROIs+1, num_elts);
        strain(:, 1) = 1;
    end
    if f > fstart
        ref = f-fstart+1;
        for j=1:numROIs+1
            %isolate components
            dx = epi_ROIs_x_OG(j)-endo_ROIs_x_OG(j);
            dy = epi_ROIs_y_OG(j)-endo_ROIs_y_OG(j);
            dist  = sqrt(dx*dx + dy*dy);
            ux = dx/dist;
            uy = dy/dist;
            
            dx = (epi_ROIs_x(j)-endo_ROIs_x(j))*ux;
            dy = (epi_ROIs_y(j)-endo_ROIs_y(j))*uy;
            
            strain(j, ref) = (sqrt(dx*dx + dy*dy))/dist;
        end
    end
    
    %% Export to .GIF %%
    
    hold off
    
    print('currentframe', '-dpng');
    fig = rgb2gray(imread('currentframe.png'));
    
    if(f == fstart) % first frame being observed
        imwrite(fig, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.3);
    else
        imwrite(fig, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.3);
    end
end

%% Export .XLSX %%

lim = length(strain(1, :));
j = 1;
while j <= lim
    if strain(:, j) == zeros(numROIs+1, 1)
        strain(:, j) = [];
        j = j-1;
        lim = lim-1;
    end
    j = j+1;
end

delete('strain.xlsx');
xlswrite('strain.xlsx', strain);
