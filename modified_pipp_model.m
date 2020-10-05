% Set density so that number of ON and OFF cells are each ~50
% [1.47:0.36=80:20, 
density_on = 0.91; 
density_off = 0.91;  


X_bound = [0 7.4 7.4 0 0]; 
Y_bound = [0 0 7.4 7.4 0]; 
Left = min(X_bound); Right = max(X_bound); 
Bottom = min(Y_bound); Top = max(Y_bound); 
W = Right - Left; 
H = Top - Bottom; 
gridsize = [W H]; 
grid_area = prod(gridsize);


% Assign # of cells 
fprintf('-----------------------')
n_on = round(density_on * grid_area);
n_off = round(density_off * grid_area); 
[n_on n_off]
fprintf('-----------------------\n')


%% Set up interaction functions & exclusion zones

% Set parameter values [0.75 0.50 0.50]:repulsion/anti-aligned. [0.68 0.45 0.55]:attraction/aligned
int_type = 'anti-aligned'; % 'aligned' / 'anti-aligned' / 'independent' 
exc_homotypic_sf = 0.80; % exclusion zone scale factor (scaling down the exclusion zone to a fraction of the median inter-cell distance)
het_rep_bound_onoff = 0.60; % in relative units of established mosaic
het_attr_bound_onoff = 0.35; % in relative units of established mosaic
jitter = 0.15; % makes bounds flexible
%num_mosaics = 1; 
num_mosaics = 100;  

pos_on_aa = zeros(n_on,2,num_mosaics);
pos_off_aa = zeros(n_off,2,num_mosaics);



% Set up interaction functions (homotypic)
h_on_on_fun = @(x,delta_onon) ( (x>delta_onon) .* (1) + (x<=delta_onon) .* (0) );  % Hard exclusion 
h_off_off_fun = @(x,delta_offoff) ( (x>delta_offoff) .* (1) + (x<=delta_offoff) .* (0) );  % Hard exclusion 

% Set up interaction functions (heterotypic)
h_onoff_repul_fun = @(x,delta_onoff) ( (x>delta_onoff) .* (1) + (x<=delta_onoff) .* (0) );  % Repulsion 
h_onoff_attr_fun = @(x,delta_onoff) ( (x<delta_onoff) .* (1) + (x>=delta_onoff) .* (0) );  % Attraction 



%% Simulate mosaic 

pos_store = cell(num_mosaics,1); 
for nsm=1:num_mosaics
    
    % Step 1: Set up ON mosaic with homotypic interaction only
    exc_mean = sqrt(sqrt(4/3) ./ (n_on ./ prod(gridsize))) * exc_homotypic_sf;
    exc_std = exc_mean * jitter; % making exclusion zone flexible by 10% of the exclusion mean
    
    pos_on = [];
    pos_on(1,:) = rand(1,2).*gridsize; % first element 
    for nc =2:n_on
        temp_pos = rand(1,2).*gridsize;
        dists = pdist2(temp_pos, pos_on);
        delta_onon = normrnd(exc_mean, exc_std);
        while ~(sum(h_on_on_fun(dists,delta_onon))==length(dists))
            temp_pos = rand(1,2).*gridsize;
            dists = pdist2(temp_pos, pos_on);
            delta_onon = normrnd(exc_mean, exc_std);
        end
        pos_on(nc,:) = temp_pos;
        %fprintf('Cell %d of %d done \n',nc,n_on);
    end
    
    
    % Step 2: Set up OFF mosaic with self-homotypic and repulsive heterotypic
    % interaction with ON mosaic
    
    exc_mean = sqrt(sqrt(4/3) ./ (n_off ./ prod(gridsize))) * exc_homotypic_sf;
    exc_std = exc_mean * jitter; % making exclusion zone flexible by 10% of the exclusion mean
    
    if strcmp(int_type, 'anti-aligned')
        
        nnd_on = ipdm(pos_on,'Subset','NearestNeighbor');
        nnd_on = median(nonzeros(nnd_on));
        mean_het_rep_bound_onoff = nnd_on*het_rep_bound_onoff;
        std_het_rep_bound_onoff = mean_het_rep_bound_onoff * jitter;
        
        pos_off = zeros(1,2);
 
        for nc =1:n_off
            temp_pos = rand(1,2).*gridsize;
            if nc==1 % only heterotypic
                dists1 = pdist2(temp_pos, pos_on);
                delta_onoff = normrnd(mean_het_rep_bound_onoff, std_het_rep_bound_onoff);
                while ~(sum(h_onoff_repul_fun(dists1,delta_onoff))==length(dists1))
                    temp_pos = rand(1,2).*gridsize;
                    dists1 = pdist2(temp_pos, pos_on);
                    delta_onoff = normrnd(mean_het_rep_bound_onoff, std_het_rep_bound_onoff);
                end
                pos_off(nc,:) = temp_pos;
                
            else % heterotypic & homotypic
                
                % heterotypic
                dists1 = pdist2(temp_pos, pos_on);
                delta_onoff = normrnd(mean_het_rep_bound_onoff, std_het_rep_bound_onoff);
                % homotypic
                dists2 = pdist2(temp_pos, pos_off);
                delta_offoff = normrnd(exc_mean, exc_std);
                
                % (1) all heterotypic probabilities should be > p_thresh, (2) new
                % position should be outside the exlusion zone of each homotypic
                % particle
                while ~( (sum(h_onoff_repul_fun(dists1,delta_onoff))==length(dists1)) && (sum(h_off_off_fun(dists2,delta_offoff))==length(dists2)) )
                    temp_pos = rand(1,2).*gridsize;
                    dists1 = pdist2(temp_pos, pos_on);
                    delta_onoff = normrnd(mean_het_rep_bound_onoff, std_het_rep_bound_onoff);
                    dists2 = pdist2(temp_pos, pos_off);
                    delta_offoff = normrnd(exc_mean, exc_std);
                end
                pos_off(nc,:) = temp_pos;
            end
            %fprintf('Cell %d of %d done \n',nc,n_off);
        end 
        % figure; cla; plot(pos_on_bt(:,1),pos_on_bt(:,2),'or',pos_off_bt(:,1),pos_off_bt(:,2),'ob'); legend('ON bt','OFF bt'); hold on
        
    elseif strcmp(int_type,'aligned')
        
        nnd_on = ipdm(pos_on,'Subset','NearestNeighbor');
        nnd_on = median(nonzeros(nnd_on));
        mean_het_attr_bound_onoff = nnd_on*het_attr_bound_onoff;
        std_het_attr_bound_onoff = mean_het_attr_bound_onoff * jitter;
        
        pos_off = zeros(1,2);
        for nc =1:n_off
            temp_pos = rand(1,2).*gridsize;
            if nc==1 % only heterotypic
                dists1 = pdist2(temp_pos, pos_on);
                delta_onoff = normrnd(mean_het_attr_bound_onoff, std_het_attr_bound_onoff);
                while ~(sum(h_onoff_attr_fun(dists1,delta_onoff))>=1)  % inter-particle distances < threshold distance
                    temp_pos = rand(1,2).*gridsize;
                    dists1 = pdist2(temp_pos, pos_on);
                    delta_onoff = normrnd(mean_het_attr_bound_onoff, std_het_attr_bound_onoff);
                end
                pos_off(nc,:) = temp_pos;
                
            else % heterotypic & homotypic
                
                % heterotypic
                dists1 = pdist2(temp_pos, pos_on);
                delta_onoff = normrnd(mean_het_attr_bound_onoff, std_het_attr_bound_onoff);
                % homotypic
                dists2 = pdist2(temp_pos, pos_off);
                delta_offoff = normrnd(exc_mean, exc_std);
                
                % (1) all heterotypic probabilities should be > p_thresh, (2) new
                % position should be outside the exlusion zone of each homotypic
                % particle
                while ~( (sum(h_onoff_attr_fun(dists1,delta_onoff))>=1) && (sum(h_off_off_fun(dists2,delta_offoff))==length(dists2)) )
                    temp_pos = rand(1,2).*gridsize;
                    dists1 = pdist2(temp_pos, pos_on);
                    delta_onoff = normrnd(mean_het_attr_bound_onoff, std_het_attr_bound_onoff);
                    dists2 = pdist2(temp_pos, pos_off);
                    delta_offoff = normrnd(exc_mean, exc_std);
                end
                pos_off(nc,:) = temp_pos;
            end
            %fprintf('Cell %d of %d done \n',nc,n_off);
        end
        
    elseif strcmp(int_type,'independent')
        
        pos_on = [];
        pos_on(1,:) = rand(1,2).*gridsize; % first element
        for nc =2:n_on
            temp_pos = rand(1,2).*gridsize;
            dists = pdist2(temp_pos, pos_on);
            delta_onon = normrnd(exc_mean, exc_std);
            while ~(sum(h_on_on_fun(dists,delta_onon))==length(dists))
                temp_pos = rand(1,2).*gridsize;
                dists = pdist2(temp_pos, pos_on);
                delta_onon = normrnd(exc_mean, exc_std);
            end
            pos_on(nc,:) = temp_pos;
        end
        pos_off = [];
        pos_off(1,:) = rand(1,2).*gridsize; % first element
        for nc =2:n_off
            temp_pos = rand(1,2).*gridsize;
            dists = pdist2(temp_pos, pos_off);
            delta_offoff = normrnd(exc_mean, exc_std);
            while ~(sum(h_off_off_fun(dists,delta_offoff))==length(dists))
                temp_pos = rand(1,2).*gridsize;
                dists = pdist2(temp_pos, pos_off);
                delta_offoff = normrnd(exc_mean, exc_std);
            end
            pos_off(nc,:) = temp_pos;
        end
        
        
    end
    
    %figure; cla; plot(pos_on(:,1),pos_on(:,2),'or',pos_off(:,1),pos_off(:,2),'ob'); legend('ON','OFF'); hold on
    
    % Store data
    pos_store{nsm}.on_centers = pos_on; 
    pos_store{nsm}.off_centers = pos_off;
    fprintf('Mosaics %d of %d done \n',nsm,num_mosaics);
    
end
 


%% Calculate energy for the simulated mosaics 


% ROI: 
%   1. 'square' (default)
%   2. 'rectangle' 
%   3. 'outline' 
%   4. 'convexhull'
%   5. 'hexagon' (for 30 micron array)
%   6. 'circle' 
roi_type = 'convexhull'; 

% Interaction energy: 
%   1. 'Coulomb' (default) ~1/r 
%   2. 'Gaussian'
%   3. 'Inverse-square' 
%   4. 'Inverse-cube' 
interaction_type = 'Inverse-square'; 

verbose = true; 
buf2use = 1.2;
Rb = 21; % number of bins+1 for radial energy    

[elin_bin_gather,elin_median_gather,elin_std_gather] = deal(zeros(1,Rb-1)); 
[Emat_smooth_zsc, Emat, Emat_smooth_gather_onoff] = deal(cell(1,length(pos_store))); 


for k=1:length(pos_store)
    
    pos_on = pos_store{k}.on_centers;
    pos_off = pos_store{k}.off_centers;

    energy_stat_ON  = get_energy_for_PIPP(pos_on,pos_off,'Move','ON','ROI_type',roi_type,'print_indiv_maps',false,...
        'print_average_map',false,'buf',buf2use,'Interaction_type',interaction_type); 
    energy_stat_OFF  = get_energy_for_PIPP(pos_on,pos_off,'Move','OFF','ROI_type',roi_type,'print_indiv_maps',false,...
        'print_average_map',false,'buf',buf2use,'Interaction_type',interaction_type); 
    
    % Calculate linear energy
    [sc_x_range_gather_onoff,sc_y_range_gather_onoff,elin_bin_gather_onoff,elin_median_gather_onoff,elin_std_gather_onoff,...
        lcx_gather_onoff,lcy_gather_onoff,vlx_gather_onoff,vly_gather_onoff] = deal([]); 
    [emat_gather_onoff,emat_smooth_gather_onoff,emat_smooth_zsc_gather_onoff] = deal([]); 
    for nr = 1:2
        if nr==1
            energy_stat = energy_stat_ON;
        elseif nr==2
            energy_stat = energy_stat_OFF;
        end
        
        energymat = energy_stat.energy_smooth_zsc;
        energymat(isnan(energymat))=0;
        sc_x_range = energy_stat.sc_x_range;
        sc_y_range = energy_stat.sc_y_range;
        global_sc_steps_mat = energy_stat.global_sc_steps_mat;
        
        r_bin_edge = linspace(0,max(sc_x_range),Rb);
        r_bin_cent = (r_bin_edge(1:end-1) + r_bin_edge(2:end))./2;
        [lin_energy, lin_energy_std] = deal(zeros(1,length(r_bin_cent)));
        
        for rb=1:length(r_bin_edge)-1
            linind = find(global_sc_steps_mat(:,:,4)>r_bin_edge(rb) & global_sc_steps_mat(:,:,4)<=r_bin_edge(rb+1));
            lin_energy(rb) = median(energymat(linind));
            lin_energy_std(rb) = robust_std(energymat(linind));
        end
        lin_rad_bin= r_bin_cent;
        nrbins = length(lin_rad_bin);
        
        elin_bin = lin_rad_bin;
        elin_median = lin_energy;
        elin_std = lin_energy_std;
        
        % Set axes units for 2D energy map
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
        
        emat_gather_onoff(:,:,nr) = energy_stat.energy; 
        emat_smooth_gather_onoff(:,:,nr) = energy_stat.energy_smooth; 
        emat_smooth_zsc_gather_onoff(:,:,nr) = energymat;
        sc_x_range_gather_onoff(1,:,nr) = sc_x_range;
        sc_y_range_gather_onoff(1,:,nr) = sc_y_range;
        elin_bin_gather_onoff(1,:,nr) = elin_bin;
        elin_median_gather_onoff(1,:,nr) = elin_median;
        elin_std_gather_onoff(1,:,nr) = elin_std;
        lcx_gather_onoff(1,:,nr) = lcx';
        lcy_gather_onoff(1,:,nr) = lcy';
        vlx_gather_onoff(1,:,nr) = vlx;
        vly_gather_onoff(1,:,nr) = vly;
    end
    Emat_smooth_zsc{k}.onoff = emat_smooth_zsc_gather_onoff; 
    Emat{k}.onoff = emat_gather_onoff;
    Emat_smooth_gather_onoff{k}.onoff = emat_smooth_gather_onoff; 
    elin_bin_gather(k,:) = mean(elin_bin_gather_onoff(1,:,:),3);
    elin_median_gather(k,:) = mean(elin_median_gather_onoff(1,:,nr),3);
    elin_std_gather(k,:) = sqrt((elin_std_gather_onoff(1,:,1).^2 + elin_std_gather_onoff(1,:,2).^2)/2);
    
end


%% Plot energy figures 

nn=1;
nnd_on = ipdm(pos_on,'Subset','NearestNeighbor'); 
nnd_on = median(nonzeros(nnd_on));
nnd_off = ipdm(pos_off,'Subset','NearestNeighbor'); 
nnd_off = median(nonzeros(nnd_off));

% Fig 1e: Plot RF center positions (48-52 : ON-OFF)
clx = {'r','b'};
FS = 20;
hf = figure('Name','Simulated: 100-100'); clf(hf); hold on;
plot(pos_store{nn}.on_centers(:,1),pos_store{nn}.on_centers(:,2),'o','color',clx{1},...
    'markersize',7,'markerfacecolor',clx{1},'markeredgecolor','none'); hold on;
theta = 1:360; 
rad = nnd_on/2; 
for cc=1:length(pos_store{nn}.on_centers)
    plot(rad.*cosd(theta) + pos_store{nn}.on_centers(cc,1),rad.*sind(theta) + pos_store{nn}.on_centers(cc,2),...
        '-','color',clx{1},'linewidth',2); 
end
plot(pos_store{nn}.off_centers(:,1),pos_store{nn}.off_centers(:,2),'o','color',clx{1},...
    'markersize',7,'markerfacecolor',clx{2},'markeredgecolor','none');
theta = 1:360; 
rad = nnd_off/2; 
for cc=1:length(pos_store{nn}.off_centers)
    plot(rad.*cosd(theta) + pos_store{nn}.off_centers(cc,1),rad.*sind(theta) + pos_store{nn}.off_centers(cc,2),...
        '-','color',clx{2},'linewidth',2); 
end
axis square;
box on;


% Fig 1: Energy map 
FS=20;
hf= figure(1); clf(hf);
EE = emat{nn}.onoff(:,:,1);
imagesc(EE); hold on;
plot([0 length(sc_x_range)+0.5],repmat((length(sc_x_range)+1)/2,1,2),'--k','linewidth',2);
plot(repmat((length(sc_y_range)+1)/2,1,2),[0 length(sc_y_range)+0.5],'--k','linewidth',2);
caxis([-1 1].*max(EE(:)));
hc = colorbar; ylabel(hc,'Z-score'); axis square;
colormap('parula'); %'redblue'
set(hc,'fontsize',FS,'linewidth',2);
xlabel('Mosaic shift');
set(gca,'fontsize',FS,'linewidth',2,'tickdir','out');
set(gca,'xtick',lcx,'xticklabel',vlx,'ytick',lcy,'yticklabel',vly);
x1=get(gca,'position');
x=get(hc,'Position');
x(3)=0.045;
set(hc,'Position',x);
set(gca,'position',x1);


% Fig 2: Radial energy curve for the 48-52 case
FS=40; alph = 0.15;
clx = [0 0 0]; % black
elin_median = elin_median_gather(nn,:);
xx_ebar = [elin_bin fliplr(elin_bin)];
yy_ebar = [elin_median - elin_std  fliplr([elin_median + elin_std])];
hf= figure(2); clf(hf); set(hf,'position',[440   378   501   420]);
plot(elin_bin, elin_median,'-','color',clx,'linewidth',3); hold on;
h=fill(xx_ebar, yy_ebar,'r');
h.FaceColor = clx; h.FaceAlpha = alph; h.EdgeColor='none';
plot([0 elin_bin(end)],[0 0],'--k','linewidth',2); hold on;
xlabel('Radial shift (r)');
ylabel('Z-score');
set(gca,'fontsize',FS,'linewidth',2,'box','off', 'tickdir','out','xlim',[0 1]);



