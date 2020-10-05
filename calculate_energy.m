%% Get positions 

tempvar = load('/Documents/posmat_rat.mat'); 
num_mosaics = size(tempvar.posmat.on_bt_off_bt,1);
[pos_mosaic1,pos_mosaic2] = deal(cell(num_mosaics,1)); 
for ii=1:num_mosaics
    pos_mosaic1{ii} = tempvar.posmat.on_bt_off_bt{ii,1}; 
    pos_mosaic2{ii} = tempvar.posmat.on_bt_off_bt{ii,2}; 
end


%% Calculate energy for mosaic pairs

% ROI: 
%   1. 'square' (default)
%   2. 'rectangle' 
%   3. 'outline' 
%   4. 'convexhull'
%   5. 'hexagon' (for 30 micron array)
%   6. 'circle' 
roi_type = 'convexhull'; 

% Interaction energy: 
%   1. 'Coulomb' (default)
%   2. 'Gaussian'
%   3. 'Inverse-square' 
%   4. 'Inverse-cube' 
interaction_type = 'Inverse-square'; 

verbose = true; 
buf2use = 1.0;
Rb = 21; % number of bins+1 for radial energy    


[elin_bin_gather,elin_median_gather,elin_std_gather] = deal(zeros(length(allcombstable),Rb-1)); 

for k=1:num_mosaics 
    pos_1 = pos_mosaic1{k};
    pos_2 = pos_mosaic2{k};
    
    energy_stat_1  = get_energy_for_PIPP(pos_1,pos_2,'Move','ON','ROI_type',roi_type,'print_indiv_maps',false,...
        'print_average_map',false,'buf',buf2use,'Interaction_type',interaction_type); 
    energy_stat_2  = get_energy_for_PIPP(pos_1,pos_2,'Move','OFF','ROI_type',roi_type,'print_indiv_maps',false,...
        'print_average_map',false,'buf',buf2use,'Interaction_type',interaction_type); 
    
    % Calculate linear energy
    [sc_x_range_gather_12,sc_y_range_gather_12,elin_bin_gather_12,elin_median_gather_12,elin_std_gather_12,...
        lcx_gather_12,lcy_gather_12,vlx_gather_12,vly_gather_12] = deal([]); 
    energymat_gather_12 = []; 
    for nr = 1:2
        if nr==1
            energy_stat = energy_stat_1;
        elseif nr==2
            energy_stat = energy_stat_2;
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
        
        energymat_gather(:,:,nr,k) = energymat;
        sc_x_range_gather_12(1,:,nr) = sc_x_range;
        sc_y_range_gather_12(1,:,nr) = sc_y_range;
        elin_bin_gather_12(1,:,nr) = elin_bin;
        elin_median_gather_12(1,:,nr) = elin_median;
        elin_std_gather_12(1,:,nr) = elin_std;
        lcx_gather_12(1,:,nr) = lcx';
        lcy_gather_12(1,:,nr) = lcy';
        vlx_gather_12(1,:,nr) = vlx;
        vly_gather_12(1,:,nr) = vly;
        
    end
    elin_bin_gather(k,:) = mean(elin_bin_gather_12(1,:,:),3);
    elin_median_gather(k,:) = mean(elin_median_gather_12(1,:,nr),3);
    elin_std_gather(k,:) = sqrt((elin_std_gather_12(1,:,1).^2 + elin_std_gather_12(1,:,2).^2)/2);
end


%% Plot 2D IMCE map and radial average energy 


% 2D IMCE map (shifting ON mosaic) 
hf= figure; clf(hf);
IMCE2d = energymat_gather(:,:,1,1); 
imagesc(IMCE2d); hold on;
plot([0 length(sc_x_range)+0.5],repmat((length(sc_x_range)+1)/2,1,2),'--k','linewidth',2);
plot(repmat((length(sc_y_range)+1)/2,1,2),[0 length(sc_y_range)+0.5],'--k','linewidth',2);
caxis([-1 1].*max(IMCE2d(:)));
hc = colorbar; ylabel(hc,'Z-score'); axis square;
colormap('parula'); 
xlabel('Mosaic shift');
set(gca,'xtick',lcx,'xticklabel',vlx,'ytick',lcy,'yticklabel',vly); 


% Radial average IMCE curve 
FS=40; alph = 0.15;
clx = [0 0 0]; % black
mean_elin_median = mean(elin_median_gather,1);
xx_ebar = [elin_bin fliplr(elin_bin)];
yy_ebar = [mean_elin_median - elin_std  fliplr([mean_elin_median + elin_std])];
hf= figure(); clf(hf); 
plot(elin_bin, mean_elin_median,'-','color',clx,'linewidth',2); hold on;
h=fill(xx_ebar, yy_ebar,'k');
h.FaceColor = clx; h.FaceAlpha = alph; h.EdgeColor='none';
plot([0 elin_bin(end)],[0 0],'--k','linewidth',2); hold on;
xlabel('Radial shift (r)');
ylabel('Z-score');
set(gca,'box','off','xlim',[0 1]);




