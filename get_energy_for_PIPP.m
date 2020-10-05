function energy_stat = get_energy_for_PIPP(final_pos_on,final_pos_off,varargin)


    % Specify list of optional parameters
    p = inputParser;
    p.addParameter('print_indiv_maps', false);   
    p.addParameter('print_average_map', false);   
    p.addParameter('ROI_type', 'square');  % 'square' (default), 'rectangle', 'convexhull', 'outline', 'hexagon' 
    p.addParameter('buf',1.2,@isnumeric); 
    p.addParameter('Interaction_type','Coulomb');  % 'Coulomb' (default), 'Gaussian', 'Inverse-square'
    p.addParameter('Move','ON'); % can be either 'ON' (default) or 'OFF'
    p.parse(varargin{:});
    prms = p.Results;
    
    frac = 0.0; % % of cells removed near the boundary 
    
    
    buf = prms.buf; % edge buffer (units of inter-cell distance) 
    gamma = 1.0; % max mosaic shift (units of inter-cell distance)
    exc_diameter = 0.2; % exclusion zone (units of inter-cell distance) 
    gaussfilt_std =  0.1; % 0.15 (units of inter-cell distance) % changed from 0.05 to 0.15
    move = prms.Move; 
    clear accum steps_mat;
    % figure; 
    
    if iscell(final_pos_on) % check if the varargin{1} is a cell 
        num_dsets = length(final_pos_on); 
    else 
        num_dsets = 1; 
    end
        
    bin = 0.04; % step size to move mosaic (unit of inter-cell-distance)
    for nm=1:num_dsets 
        
        if iscell(final_pos_on)
            pos_on = final_pos_on{nm}; 
            pos_off = final_pos_off{nm}; 
        else
            pos_on = final_pos_on; 
            pos_off = final_pos_off; 
        end
        
        
        % Median cell-to-cell distance 
        temp_on = pos_on;
        temp_off = pos_off;
        c2c_dist_on = nonzeros(ipdm(temp_on,'Subset','NearestNeighbor'));
        c2c_dist_off = nonzeros(ipdm(temp_off,'Subset','NearestNeighbor'));
        med_c2c_dist_on = median(c2c_dist_on);
        med_c2c_dist_off = median(c2c_dist_off);
        
        
        % Define lattice spacing in units of NN distance 
        %     Note: 
        %           The average energy map is obtained by averaging the energy 
        %           maps from different data sets (or preps). However in each 
        %           prep, the NN cell-cell distance can be different, which 
        %           means that the absolute mosaic shift distance to cover 1 NN
        %           distance can be different in different data sets. So a 
        %           constant shift applied uniformly to all data sets would 
        %           yield energy maps with X-Y axes out of scale wrt each other
        %           in units of the dimensionless inter-cell distance. To avoid
        %           this, we shift mosaics in each data by ~1.5*NN distance 
        %           corresponding to the data set.  
        if strcmp(move,'ON') 
            d_true = med_c2c_dist_off; % in stixels 
        elseif strcmp(move,'OFF') 
            d_true = med_c2c_dist_on; % in stixels 
        end
        scale_fac = d_true; % for scaling the shifts to the original lattice units 
        d_true = mean([med_c2c_dist_off med_c2c_dist_on]); 
        scale_fac = d_true;
        
        % Set ROI 
        if strcmp(prms.ROI_type,'square') % (Default)
            
            % Get roi (without buffer)
            Left = max([min(pos_on(:,1)) min(pos_off(:,1))]);
            Right = min([max(pos_on(:,1)) max(pos_off(:,1))]);
            Bottom = max([min(pos_on(:,2)) min(pos_off(:,2))]);
            Top = min([max(pos_on(:,2)) max(pos_off(:,2))]);
            ROIx = [Left Left Right Right Left]';
            ROIy = [Bottom Top Top Bottom Bottom]';
            
            buffunit = mean([med_c2c_dist_on med_c2c_dist_off]); 
            BufferX = buf*buffunit; % Define outer bound of ROI with buffer (in stixel units)
            BufferY = buf*buffunit;
            locLeft = find(ROIx==min(ROIx));
            locRight = find(ROIx==max(ROIx));
            ROIx(locLeft) = ROIx(locLeft)+BufferX;
            ROIx(locRight) = ROIx(locRight)-BufferX;
            locBottom = find(ROIy==min(ROIy));
            locTop = find(ROIy==max(ROIy));
            ROIy(locBottom) = ROIy(locBottom)+BufferY;
            ROIy(locTop) = ROIy(locTop)-BufferY;
            
            Left = min(ROIx);
            Right = max(ROIx);
            Bottom = min(ROIy);
            Top = max(ROIy);
            W = Right - Left;
            H = Top - Bottom;
            
            if W>H  % Making ROI a square 
                locLeft = find(ROIx==min(ROIx));
                locRight = find(ROIx==max(ROIx));
                ROIx(locLeft) = ROIx(locLeft) + (W-H)/2;
                ROIx(locRight) = ROIx(locRight) - (W-H)/2; 
            elseif W<H 
                locBottom = find(ROIy==min(ROIy));
                locTop = find(ROIy==max(ROIy)); 
                ROIy(locBottom) = ROIy(locBottom) + (H-W)/2;
                ROIy(locTop) = ROIy(locTop) - (H-W)/2;
            end
            
        elseif strcmp(prms.ROI_type,'rectangle')
            
            % Get roi (without buffer)
            Left = max([min(pos_on(:,1)) min(pos_off(:,1))]);
            Right = min([max(pos_on(:,1)) max(pos_off(:,1))]);
            Bottom = max([min(pos_on(:,2)) min(pos_off(:,2))]);
            Top = min([max(pos_on(:,2)) max(pos_off(:,2))]);
            ROIx = [Left Left Right Right Left]';
            ROIy = [Bottom Top Top Bottom Bottom]';
            
            buffunit = mean([med_c2c_dist_on med_c2c_dist_off]); 
            BufferX = buf*buffunit; % Define outer bound of ROI with buffer (in stixel units)
            BufferY = buf*buffunit;
            locLeft = find(ROIx==min(ROIx));
            locRight = find(ROIx==max(ROIx));
            ROIx(locLeft) = ROIx(locLeft)+BufferX;
            ROIx(locRight) = ROIx(locRight)-BufferX;
            locBottom = find(ROIy==min(ROIy));
            locTop = find(ROIy==max(ROIy));
            ROIy(locBottom) = ROIy(locBottom)+BufferY;
            ROIy(locTop) = ROIy(locTop)-BufferY;
            
            Left = min(ROIx);
            Right = max(ROIx);
            Bottom = min(ROIy);
            Top = max(ROIy);
            W = Right - Left;
            H = Top - Bottom;
            
        elseif strcmp(prms.ROI_type,'convexhull')
            
            % Step 1: For each point detect k=6 nearest neighbors 
            % Step 2: Calculate median of the distance: dM6
            % Step 3: Calculate std of the inter-cell distance for all points : stdMAll
            % Step 4: If dM6 is >= 3 * stdMAll, then reject the point
            % Step 5: Determine "convexhull polygon" that encloses all included points 
            % Step 6: Repeat the process for ON and OFF mosaics; find intersection
            %         polygon 
            % Step 7: Scale the polygon down to include a certain buffer 

            thresh = 7; % units of std
            pos_store{1} = pos_on;
            pos_store{2} = pos_off;
            
            [Polyshape,~] = deal(cell(1,length(pos_store)));
            parfor i=1:length(pos_store)
                pos_temp = pos_store{i};
                [Idx,D] = knnsearch(pos_temp,pos_temp,'k',7);
                Idx(:,1) = [];
                D(:,1) = [];
                dM6 = mean(D,2);
                c2c_dist = nonzeros(ipdm(pos_temp,'Subset','NearestNeighbor'));
                std_dMAll = std(c2c_dist);
                mean_dMMAll = mean(c2c_dist);
                conf = dM6 < std_dMAll*thresh+ mean_dMMAll; % condition for inclusion
                
                pos_temp = pos_temp(conf,:); % after removing outliers
                
                DT = delaunayTriangulation(pos_temp);
                Ch = convexHull(DT);
                Poly = [DT.Points(Ch,1) DT.Points(Ch,2)];
                Polyshape{i} = polyshape(Poly(:,1),Poly(:,2));
            end
            Polyshape_intersect = intersect(Polyshape{1},Polyshape{2});
            Polyshape_intersect.Vertices(find(isnan(Polyshape_intersect.Vertices),1,'first'):end,:) = []; 
            
            [centroid_x,centroid_y] = centroid(Polyshape_intersect);
            cents = [centroid_x,centroid_y];
            tempvec = repmat(cents,size(Polyshape_intersect.Vertices(:,1),1),1) - Polyshape_intersect.Vertices;
            
            poly_radius = median(hypot(tempvec(:,1),tempvec(:,2))); % median radius of the polygon
            Buffer = buf*d_true;
            scaledown = (poly_radius - Buffer)/poly_radius;
            if scaledown<0 
                Polyshape_final = Polyshape_intersect; 
            else 
                Polyshape_final = scale(Polyshape_intersect,scaledown,cents); % this needs to be changed
            end

            ROIx = [Polyshape_final.Vertices(:,1); Polyshape_final.Vertices(1,1)];
            ROIy = [Polyshape_final.Vertices(:,2); Polyshape_final.Vertices(1,2)];
            
%             %------------- Previous version -----------------------------
%             % Calculate convex hull using Delaunay triangulation
%             pos = [pos_on; pos_off];
%             DT = delaunayTriangulation(pos);
%             Ch = convexHull(DT);
%             Poly = [DT.Points(Ch,1) DT.Points(Ch,2)];
%             Polyshape = polyshape(Poly(:,1),Poly(:,2));
%             Polysect = Polyshape; 
%             [centroid_x,centroid_y] = centroid(Polysect);
%             cents = [centroid_x,centroid_y]; 
%             
%             tempvec = repmat(cents,size(Polysect.Vertices,1),1) - [Polysect.Vertices(:,1) Polysect.Vertices(:,2)];
%             ddist = prctile((hypot(tempvec(:,1),tempvec(:,2))),95); % 95% of the distances from the center to edges
%             scaleup = ddist;
%             
%             Buffer = buf*d_true; 
%             scaledown = (scaleup-Buffer)/scaleup;
% %             if scaledown<0.60 
% %                 scaledown = 0.60; 
% %             end
%             
%             Polysect = scale(Polysect,scaledown,[centroid_x centroid_y]); % this needs to be changed
%             ROIx = [Polysect.Vertices(:,1); Polysect.Vertices(1,1)];
%             ROIy = [Polysect.Vertices(:,2); Polysect.Vertices(1,2)];
            
        elseif strcmp(prms.ROI_type,'outline')
            
            % Step 1: For each point detect k=6 nearest neighbors 
            % Step 2: Calculate median of the distance: dM6
            % Step 3: Calculate std of the inter-cell distance for all points : stdMAll
            % Step 4: If dM6 is >= 3 * stdMAll, then reject the point
            % Step 5: Determine "outline polygon" that encloses all included points 
            % Step 6: Repeat the process for ON and OFF mosaics; find intersection
            %         polygon 
            % Step 7: Scale the polygon down to include a certain buffer
            
            thresh = 7; % units of std
            pos_store{1} = pos_on;
            pos_store{2} = pos_off;
            
            [Polyshape,~] = deal(cell(1,length(pos_store)));
            parfor i=1:length(pos_store)
                pos_temp = pos_store{i};
                [Idx,D] = knnsearch(pos_temp,pos_temp,'k',7);
                Idx(:,1) = [];
                D(:,1) = [];
                dM6 = mean(D,2);
                c2c_dist = nonzeros(ipdm(pos_temp,'Subset','NearestNeighbor'));
                std_dMAll = std(c2c_dist);
                mean_dMMAll = mean(c2c_dist);
                conf = dM6 < std_dMAll*thresh+ mean_dMMAll; % condition for inclusion
                
                pos_temp = pos_temp(conf,:); % after removing outliers
                
                k = boundary(pos_temp(:,1),pos_temp(:,2));
                Polyshape{i} = polyshape(pos_temp(k,1),pos_temp(k,2));  
            end
            Polyshape_intersect = intersect(Polyshape{1},Polyshape{2});
            Polyshape_intersect.Vertices(find(isnan(Polyshape_intersect.Vertices),1,'first'):end,:) = []; 
            
            [centroid_x,centroid_y] = centroid(Polyshape_intersect);
            cents = [centroid_x,centroid_y];
            tempvec = repmat(cents,size(Polyshape_intersect.Vertices(:,1),1),1) - Polyshape_intersect.Vertices;
            
            poly_radius = median(hypot(tempvec(:,1),tempvec(:,2))); % median radius of the polygon
            Buffer = buf*d_true;
            scaledown = (poly_radius - Buffer)/poly_radius;
            if scaledown<0 
                Polyshape_final = Polyshape_intersect; 
            else 
                Polyshape_final = scale(Polyshape_intersect,scaledown,cents); % this needs to be changed
            end
            
            ROIx = [Polyshape_final.Vertices(:,1); Polyshape_final.Vertices(1,1)];
            ROIy = [Polyshape_final.Vertices(:,2); Polyshape_final.Vertices(1,2)];
            
%             % ------------- Previous version ----------------------------           
%             pos = [pos_on; pos_off];
%             k = boundary(pos(:,1),pos(:,2));
%             pgon = polyshape(pos(k,1),pos(k,2));             
%             [centroid_x,centroid_y] = centroid(pgon);
%             cents = [centroid_x,centroid_y]; 
%             tempvec = repmat(cents,length(k),1) - [pos(k,1) pos(k,2)];
%             ddist = prctile((hypot(tempvec(:,1),tempvec(:,2))),95); % 95% of the distances from the center to edges
%             scaleup = ddist;
%             
%             Buffer = buf*d_true; 
%             scaledown = (scaleup-Buffer)/scaleup;
%             if scaledown<0.60 
%                 scaledown = 0.60; 
%             end
%             
%             Polysect = scale(pgon,scaledown,[centroid_x centroid_y]); % this needs to be changed
%             ROIx = [Polysect.Vertices(:,1); Polysect.Vertices(1,1)];
%             ROIy = [Polysect.Vertices(:,2); Polysect.Vertices(1,2)];

            
        elseif strcmp(prms.ROI_type,'hexagon')
            
            pos = [pos_on; pos_off];
            cents = mean(pos,1);
            k = boundary(pos(:,1),pos(:,2));
            tempvec = repmat(cents,length(k),1) - [pos(k,1) pos(k,2)];
            ddist = prctile((hypot(tempvec(:,1),tempvec(:,2))),90); % 90% of the distances from the center to edges
            
            scaleup = ddist; % this is used to scale up the hexagon
            n_sides = 6;
            t=(1/(n_sides*2):1/n_sides:1)'*2*pi;
            x=sin(t+(pi/2));
            y=cos(t+(pi/2));
            x=scaleup*[x; x(1)]; x = x+cents(1); 
            y=scaleup*[y; y(1)]; y = y+cents(2); 
            pgon = polyshape(x,y); 
            [centroid_x,centroid_y] = centroid(pgon);
            
            Buffer = buf*d_true;
            scaledown = (scaleup-Buffer)/scaleup;
            if scaledown<0.60
                scaledown = 0.60;
            end
            
            Polysect = scale(pgon,scaledown,[centroid_x centroid_y]); % this needs to be changed
            ROIx = [Polysect.Vertices(:,1); Polysect.Vertices(1,1)];
            ROIy = [Polysect.Vertices(:,2); Polysect.Vertices(1,2)];
            
        elseif strcmp(prms.ROI_type,'circle')
            
            pos = [pos_on; pos_off];
            cents = [median(pos(:,1)) median(pos(:,2))]; 
            DT = delaunayTriangulation(pos(:,1),pos(:,2)); 
            [C,~] = convexHull(DT); 
            xch = DT.Points(C,1); 
            ych = DT.Points(C,2); 
            xch = xch(:); ych = ych(:); 
            A = [-2*xch -2*ych ones(length(xch),1)];
            xch = A\-(xch.^2+ych.^2);
            xo=xch(1);
            yo=xch(2);
            R = sqrt(  xo.^2 + yo.^2  - xch(3));
            Rscaleddown = R - buf*d_true;
            [ROIx,ROIy] = scircle1(cents(1),cents(2),Rscaleddown);
        end
            
        
%         % Remove cells +/- 1.2 icd around the boundary
%         % -----------------------------------------------------------------
%         tempd = mean([med_c2c_dist_off  med_c2c_dist_on]).*1.2;  
%         ROIx_inner = ROIx; ROIx_outer = ROIx; 
%         ROIy_inner = ROIy; ROIy_outer = ROIy; 
%         
%         locLeft = find(ROIx==min(ROIx));
%         locRight = find(ROIx==max(ROIx));
%         locBottom = find(ROIy==min(ROIy));
%         locTop = find(ROIy==max(ROIy));
%             
%         
%         ROIx_inner(locLeft) = ROIx_inner(locLeft)+tempd;
%         ROIx_inner(locRight) = ROIx_inner(locRight)-tempd;
%         ROIy_inner(locBottom) = ROIy_inner(locBottom)+tempd;
%         ROIy_inner(locTop) = ROIy_inner(locTop)-tempd;
%         
%         ROIx_outer(locLeft) = ROIx_outer(locLeft)-tempd;
%         ROIx_outer(locRight) = ROIx_outer(locRight)+tempd;
%         ROIy_outer(locBottom) = ROIy_outer(locBottom)-tempd;
%         ROIy_outer(locTop) = ROIy_outer(locTop)+tempd;
%         
%         
%         pos_on_excised = pos_on; 
%         pos_off_excised = pos_off; 
%         out_on = ~inpolygon(pos_on_excised(:,1), pos_on_excised(:,2), ROIx, ROIy); 
%         out_off = ~inpolygon(pos_off_excised(:,1), pos_off_excised(:,2), ROIx, ROIy); 
%         
%         figure; set(gcf,'position',[440   551   603   247]); 
%         subplot(1,2,1); 
%         plot(pos_on_excised(:,1),pos_on_excised(:,2),'.r',pos_off_excised(:,1),pos_off_excised(:,2),'.b'); hold on; 
%         plot(ROIx,ROIy,'-k'); 
%         title({'nON=102, nOFF=102';'No cells removed at ROI boundary'});
%         set(gca,'xlim',[-5 1400],'ylim',[-5 700]); 
%         
%         index_on = find(out_on); 
%         index_off = find(out_off); 
%         index_on__ = index_on(randperm(sum(out_on), round(frac*sum(out_on))));
%         index_off__ = index_off(randperm(sum(out_off), round(frac*sum(out_off))));
%         
%         pos_on_excised(index_on__,:) = []; 
%         pos_off_excised(index_off__,:) = [];
%         
%         out_on_ = ~inpolygon(pos_on_excised(:,1), pos_on_excised(:,2), ROIx, ROIy); 
%         out_off_ = ~inpolygon(pos_off_excised(:,1), pos_off_excised(:,2), ROIx, ROIy);
%         
%         
%         subplot(1,2,2); 
%         plot(pos_on_excised(:,1),pos_on_excised(:,2),'.r',pos_off_excised(:,1),pos_off_excised(:,2),'.b'); hold on; 
%         plot(pos_on(index_on__,1),pos_on(index_on__,2),'or',pos_off(index_off__,1),pos_off(index_off__,2),'ob'); hold on; 
%         plot(ROIx,ROIy,'-k'); 
%         title({'nON=102, nOFF=102';'50% of ON, 50% of OFF cells removed at ROI boundary'});
%         set(gca,'xlim',[-5 1400],'ylim',[-5 700]); 
%         
%         
%         % Median cell-to-cell distance 
%         temp_on = pos_on_excised;
%         temp_off = pos_off_excised;
%         c2c_dist_on = nonzeros(ipdm(temp_on,'Subset','NearestNeighbor'));
%         c2c_dist_off = nonzeros(ipdm(temp_off,'Subset','NearestNeighbor'));
%         med_c2c_dist_on = median(c2c_dist_on);
%         med_c2c_dist_off = median(c2c_dist_off);
%         if strcmp(move,'ON') 
%             d_true = med_c2c_dist_off; % in stixels 
%         elseif strcmp(move,'OFF') 
%             d_true = med_c2c_dist_on; % in stixels 
%         end
%         scale_fac = d_true; % for scaling the shifts to the original lattice units 
% 
%         
%         
%         
%         % -----------------------------------------------------------------
%         
        

        % Set up lattice shift parameters, ROI, and grid matrix 
        x_range = scale_fac.*[-gamma:bin:gamma]; % shift steps (units of stixel)
        y_range = x_range;
        num_steps = length(x_range);

        for x_step = 1:num_steps
            for y_step = 1:num_steps 
                steps_mat(x_step,y_step,1) = x_range(x_step); % Z dimension: X,Y,r,theta
                steps_mat(x_step,y_step,2) = y_range(y_step); 
                [th,rd] = cart2pol(x_range(x_step),y_range(y_step)); 
                steps_mat(x_step,y_step,3) = th;
                steps_mat(x_step,y_step,4) = rd;
            end
        end

        % Minimum distance 2 cells can be within each other 
        min_sep = exc_diameter*d_true; 

        % Calculate energy 
        [ON_pos_temp,OFF_pos_temp] = deal([]); 
        [energy,energy_nonnorm] = deal(zeros(num_steps,num_steps)); 
        [countsin_on, countsin_off, countsdiff_on] = deal(zeros(num_steps,num_steps));  
        ON_pos_initial = temp_on;
        OFF_pos_initial = temp_off; 


    %         figure(100); clf; 
    %         plot(ON_pos_initial(:,1),ON_pos_initial(:,2),'or'); hold on;
    %         plot(OFF_pos_initial(:,1),OFF_pos_initial(:,2),'ob'); hold on;
    %         plot(bound(:,1),bound(:,2),'-k'); 
    %         plot(ROIx,ROIy,'--k'); 
    %         pause();


        %figure; 
        
        in_true = inpolygon(pos_on(:,1),pos_on(:,2), ROIx, ROIy); 
        std_gaussian = mean([median(nonzeros(ipdm(pos_on,'Subset','NearestNeighbor')))   median(nonzeros(ipdm(pos_off,'Subset','NearestNeighbor')))])/2; % for Gaussian interaction 
        
        for x_step = 1:num_steps
            for y_step = 1:num_steps 

                if strcmp(move,'ON')  
                    % freeze on cell mosaic
                    ON_pos_temp = [ON_pos_initial(:,1)+x_range(x_step) ON_pos_initial(:,2)+y_range(y_step)];
                    % shift off cell mosaic 
                    OFF_pos_temp = OFF_pos_initial;
                elseif strcmp(move,'OFF') 
                    % freeze on cell mosaic
                    ON_pos_temp = ON_pos_initial;
                    % shift off cell mosaic 
                    OFF_pos_temp = [OFF_pos_initial(:,1)+x_range(x_step) OFF_pos_initial(:,2)+y_range(y_step)];
                end

                % Get the cells inside ROI 
                in_on_updated = inpolygon(ON_pos_temp(:,1), ON_pos_temp(:,2), ROIx, ROIy); 
                in_off_updated = inpolygon(OFF_pos_temp(:,1), OFF_pos_temp(:,2), ROIx, ROIy); 
                countsin_on(x_step,y_step) = sum(in_on_updated);
                countsin_off(x_step,y_step) = sum(in_off_updated); 
                
                % Get the change in number of cells 
                countsdiff_on(x_step,y_step) = sum(in_true) - sum(in_on_updated) ; 
                
                % Get heterotypic pair-wise distance 
                pairwise_dists = ipdm(ON_pos_temp(in_on_updated,:), OFF_pos_temp(in_off_updated,:)); % row: ON, col: OFF            
                

                % Calculate energy
                if strcmp(prms.Interaction_type,'Coulomb')
                    % Set energy of particles closer than "s" soma diameter to energy at soma diameter separation.
                    pairwise_dists(pairwise_dists<min_sep) = min_sep;
                    energy_temp = (1./pairwise_dists);
                    energy(x_step,y_step) = sum(energy_temp(:))./length(energy_temp(:));
                    %energy(x_step,y_step) = sum(energy_temp(:));
                    
                elseif strcmp(prms.Interaction_type,'Inverse-square')
                    % Set energy of particles closer than "s" soma diameter to energy at soma diameter separation.
                    pairwise_dists(pairwise_dists<min_sep) = min_sep;
                    energy_temp = (1./pairwise_dists.^2);
                    energy(x_step,y_step) = sum(energy_temp(:))./length(energy_temp(:));
                    
                elseif strcmp(prms.Interaction_type,'Gaussian')
                    energy_temp = (1/(2*pi*std_gaussian^2)).*exp(- pairwise_dists.^2./(2*std_gaussian^2));
                    energy(x_step,y_step) =sum(energy_temp(:))./length(energy_temp(:));
                    
                end
                
                energy_nonnorm(x_step,y_step) = sum(energy_temp(:)); 
                
    %             cla(gca); 
    %             plot(ON_pos_temp(:,1),ON_pos_temp(:,2),'o','color',clr(1,:)); hold on; 
    %             plot(OFF_pos_temp(:,1),OFF_pos_temp(:,2),'o','color',clr(2,:)); hold on; 
    %             plot(ON_pos_temp(in_on_updated,1),ON_pos_temp(in_on_updated,2),'^m');
    %             plot(OFF_pos_temp(in_off_updated,1),OFF_pos_temp(in_off_updated,2),'^m');
    %             plot(ROIx,ROIy,'--k'); axis equal; 
    %             set(gca,'xlim',[min(bound(:,1))-s max(bound(:,1))+s],'ylim',[min(bound(:,2))-s max(bound(:,2))+s]);
    %             pause(0.01);

            end
        end
        %fprintf('Data %d out of %d \n',nm,num_dsets); 
        
        accum.energy_map(:,:,nm) = energy; 
        accum.energy_nonnorm(:,:,nm) = energy_nonnorm; 
        
        energy_smooth = imgaussfilt(energy, [ceil(gaussfilt_std/bin) ceil(gaussfilt_std/bin)]); % gaussfilt_std is the std of the Gaussian filter
        energy_smooth_zsc = (energy_smooth - mean(energy_smooth(:)))./robust_std(energy_smooth(:)); 
        
        accum.energy_map_smooth(:,:,nm) = energy_smooth; 
        accum.energy_map_smooth_zsc(:,:,nm) = energy_smooth_zsc; 
        
        %****
        
        energy_nonnorm_zsc = (energy_nonnorm - mean(energy_nonnorm(:)))./robust_std(energy_nonnorm(:)); 
        
        % Gather count in ROI
        accum.countsin_on{nm} = countsin_on; 
        accum.countsin_off{nm} = countsin_off; 

        % Gather COMs and ROIs
        accum.ON_pos_initial{nm} = ON_pos_initial; 
        accum.OFF_pos_initial{nm} = OFF_pos_initial; 
        accum.ROIx{nm} = ROIx; 
        accum.ROIy{nm} = ROIy; 
        accum.d_true{nm} = d_true; 
        accum.x_range{nm} = x_range;
    end


    % Gather matrix of X,Y,rho,theta corresponding to each shift (units of unitless NN distance)
    sc_x_range = -gamma:bin:gamma;
    sc_y_range = sc_x_range;
    global_sc_steps_mat = zeros(length(sc_x_range),length(sc_x_range),4); 
    for x_step = 1:num_steps
        for y_step = 1:num_steps 
            global_sc_steps_mat(x_step,y_step,1) = sc_x_range(x_step); % Z dimension: X,Y,r,theta
            global_sc_steps_mat(x_step,y_step,2) = sc_x_range(y_step); 
            [th,rd] = cart2pol(sc_x_range(x_step),sc_x_range(y_step)); 
            global_sc_steps_mat(x_step,y_step,3) = th;
            global_sc_steps_mat(x_step,y_step,4) = rd;
        end
    end

    % Average smoothed Z scored energy maps
    energy_map = zeros(size(accum.energy_map,1), size(accum.energy_map,2)); 
    for nm=1:size(accum.energy_map,3) 
        temp = accum.energy_map(:,:,nm); 
        temp = imgaussfilt(temp, [ceil(gaussfilt_std/bin) ceil(gaussfilt_std/bin)]); % smooth  
        if robust_std(temp(:))~=0
            temp = (temp-mean(temp(:)))./robust_std(temp(:)); % z-score
        end
        energy_map = energy_map + temp; 
    end
    mean_smoothed_energy_map = energy_map./size(accum.energy_map,3); 
    mean_smooothed_energy_map_zsc = (mean_smoothed_energy_map-mean(mean_smoothed_energy_map(:)))./robust_std(mean_smoothed_energy_map(:));
   
    
    % Set axes units for image 
    x_un = sc_x_range;
    y_un = fliplr(sc_y_range);
    tickmarks = [-1 0 1];
    [~,lcx] = intersect(x_un, tickmarks); 
    vlx = sc_x_range(lcx); 
    [~,lcy] = intersect(y_un, tickmarks); 
    if ~isequal(lcy,sort(lcy,'ascend'))
        [lcy,invrt] = sort(lcy,'ascend'); 
        vly = tickmarks(invrt); 
    end
    
    
    % Get radial energy dependence for average map 
    r_bin_edge = linspace(0,max(sc_x_range),21);
    r_bin_cent = (r_bin_edge(1:end-1) + r_bin_edge(2:end))./2;
    r_max = r_bin_edge(end);
    [lin_energy_avgmap, lin_energy_std_avgmap] = deal(zeros(1,length(r_bin_cent))); 
    for rb=1:length(r_bin_edge)-1
        linind = find(global_sc_steps_mat(:,:,4)>r_bin_edge(rb) & global_sc_steps_mat(:,:,4)<=r_bin_edge(rb+1));
        lin_energy_avgmap(rb) = median(mean_smooothed_energy_map_zsc(linind));
        lin_energy_std_avgmap(rb) = std(mean_smooothed_energy_map_zsc(linind)); 
    end
    rads = global_sc_steps_mat(:,:,4); 
    [sorted_rads, sorted_inds] = sort(rads(:),'ascend');
    sorted_energies = mean_smooothed_energy_map_zsc(sorted_inds); 

    
    %% Plot individual & average energy maps and radial energy dependence  
    
    store_smooth_zsc_energy = cell(1,num_dsets);
    r_bin_edge = linspace(0,max(sc_x_range),21);
    r_bin_cent = (r_bin_edge(1:end-1) + r_bin_edge(2:end))./2;
    [lin_energy, lin_energy_std, lin_rad_bin] = deal(cell(1,num_dsets));
    for nm=1:num_dsets
        % Get 2D energy map  
        localenergy = accum.energy_map(:,:,nm); 
        localenergy_smooth = imgaussfilt(localenergy, [ceil(gaussfilt_std/bin) ceil(gaussfilt_std/bin)]); % gaussian filter with sigma = 1/4*inter-cell dist
        if robust_std(localenergy_smooth(:))~=0        
            localenergy_smooth_zsc = (localenergy_smooth-mean(localenergy_smooth(:)))./robust_std(localenergy_smooth(:)); 
        end
        store_smooth_zsc_energy{nm} = localenergy_smooth_zsc; 
        
        
        % Get radial energy dependence  
        lin_energy{nm} = zeros(1,length(r_bin_cent)); 
        lin_energy_std{nm} = zeros(1,length(r_bin_cent)); 
        for rb=1:length(r_bin_edge)-1
            linind = find(global_sc_steps_mat(:,:,4)>r_bin_edge(rb) & global_sc_steps_mat(:,:,4)<=r_bin_edge(rb+1));
            lin_energy{nm}(rb) = median(store_smooth_zsc_energy{nm}(linind)); 
            lin_energy_std{nm}(rb) = robust_std(store_smooth_zsc_energy{nm}(linind)); 
        end
        lin_rad_bin{nm} = r_bin_cent; 
    end
    
    
    if prms.print_indiv_maps
        % Figure of energy maps
        hf= figure(); clf(hf); %set(hf,'position',[-1536 385 1413 310]); 
        nc = ceil(sqrt(num_dsets));
        nr = ceil(num_dsets/nc); 
        
        for nm=1:num_dsets
            subplot(nr,nc,nm);
            localenergy_smooth_zsc = store_smooth_zsc_energy{nm}; 
            h = pcolor(global_sc_steps_mat(:,:,1),global_sc_steps_mat(:,:,2), localenergy_smooth_zsc); hold on; 
            plot([sc_x_range(1) sc_x_range(end)],[0 0],'-w',[0 0],[sc_y_range(1) sc_y_range(end)],'-w');
            colorbar; 
            set(h,'EdgeColor','None'); axis square; 
            title(['Run ',num2str(nm)]); 
            %set(gca,'fontsize',12);
        end
        %ylabel('Mosaic shift (unit of d)');
    
    
        % Figure of radial energy dependence 
        hf= figure(); clf(hf); 
        nc = ceil(sqrt(num_dsets));
        nr = ceil(num_dsets/nc); 
        for nm=1:num_dsets
            subplot(nr,nc,nm);
            % scatter(sorted_rads, sorted_energies,15,...
            %     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3); hold on; 
            errorbar(lin_rad_bin{nm}, lin_energy{nm},lin_energy_std{nm}./2, 'o-b','linewidth',2); hold on; 
            plot([0 lin_rad_bin{nm}(end)],[0 0],'--k','linewidth',1); hold on; 
            if nm==(nr-1)*nc+1
                xlabel('Radial dist (unit of inter-cell dist)');
                ylabel('Z-scored energy');
                set(gca,'fontsize',12);
            end
            title(['Run ',num2str(nm)]); 
        end
    end
    
    if prms.print_average_map
        % Figure of average energy map 
        figure;
        imagesc(mean_smooothed_energy_map_zsc); hold on; 
        plot([1 length(sc_x_range)],repmat((length(sc_x_range)+1)/2,1,2),'-k','linewidth',1);
        plot(repmat((length(sc_y_range)+1)/2,1,2),[1 length(sc_y_range)],'-k','linewidth',1);
        caxis([-1 1].*max(mean_smooothed_energy_map_zsc(:))); 
        hc = colorbar; ylabel(hc,'Z-score'); axis square;
        set(gca,'xtick',lcx,'xticklabel',vlx,'ytick',lcy,'yticklabel',vly); 


        % Plot radial energy 
        figure;
        errorbar(r_bin_cent, lin_energy_avgmap,lin_energy_std_avgmap./2,'o-b','linewidth',2); hold on; 
        plot([0 r_bin_cent(end)],[0 0],'--k','linewidth',2); 
        xlabel('Mosaic shift (unit of inter-cell dist)');
        ylabel('Z-score');
        title('Energy profile');
    end

   
    % Prepare return stats 
    energy_stat.ROIx = ROIx; 
    energy_stat.ROIy = ROIy; 
    energy_stat.pos_on = pos_on;
    energy_stat.pos_off = pos_off;
    energy_stat.sc_x_range = sc_x_range; 
    energy_stat.sc_y_range = sc_y_range; 
    energy_stat.global_sc_steps_mat = global_sc_steps_mat; 
    energy_stat.lin_energy = lin_energy{nm};
    energy_stat.lin_rad_bin = lin_rad_bin{nm};
    energy_stat.lin_energy_std = lin_energy_std{nm};
    
    energy_stat.store_smooth_zsc_energy = store_smooth_zsc_energy{nm};
    energy_stat.countsdiff_on = countsdiff_on;
    energy_stat.countsin_on = countsin_on;
    energy_stat.countsin_off = countsin_off;
    
    energy_stat.energy = energy;
    energy_stat.energy_smooth = energy_smooth;
    energy_stat.energy_smooth_zsc = energy_smooth_zsc;
    
    energy_stat.energy_nonnorm = energy_nonnorm; 
    energy_stat.energy_nonnorm_zsc = energy_nonnorm_zsc; 
    
    
    
end


