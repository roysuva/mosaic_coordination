%% Load data 

positions = load('/Documents/posmat_rat.mat');
dsub = 1:size(positions.posmat.on_bt_off_bt,1);

ctype1 = 1; % ON=1, OFF=2
ctype2 = 2; % ON=1, OFF=2
ctype_sub1 = 2; % 1=brisk sustained, 2=brisk transient
ctype_sub2 = 2; % 1=brisk sustained, 2=brisk transient

%% Get combinations 

% #########################################################################
% ctype1 : ON, OFF
% ctype2 : ON, OFF
% ctype_sub1 : 1=brisk sustained, 2=brisk transient
% ctype_sub2 : 1=brisk sustained, 2=brisk transient
% #########################################################################


% Set up combinations  
onoffcombs = [];
for i=1:length(dsub)
    onoffcombs = [onoffcombs; dsub(i) ctype1 ctype_sub1; dsub(i) ctype2 ctype_sub2]; 
end
n_onoffcombs = size(onoffcombs,1); 

allcombs = nchoosek(1:n_onoffcombs, 2); 
allcombstable = zeros(size(allcombs,1),size(onoffcombs,2)*2); 
% col 1: data, col 2: ON/OFF(1/2), col 3: cell-type 
% col 4: data, col 5: ON/OFF(1/2), col 6: cell-type  
locate_ONOFF_samerec = zeros(length(dsub),1);
cnt = 1; 
for ii=1:size(allcombs,1)
    allcombstable(ii,:) = [onoffcombs(allcombs(ii,1),:) onoffcombs(allcombs(ii,2),:)]; 
    if allcombstable(ii,1)==allcombstable(ii,4)
        locate_ONOFF_samerec(cnt) = ii;
        cnt = cnt+1; 
    end
end


%% Generate pairwise combinations, store outlines of RF 

[store_pos_1,store_pos_2,store_sd_1,store_sd_2,store_angle_1,store_angle_2,...
    store_X_1,store_Y_1,store_X_2,store_Y_2] = deal(cell(length(allcombstable),1)); 
if ~exist('htype','var')
    htype = {'ON','OFF'}; 
end


for k=1:size(allcombstable,1)

    n = allcombstable(k,1); % data 1
    nty = allcombstable(k,2); % ON/OFF type 1
    ctype_sub1 = allcombstable(k,3); % cell type 
    
    pos1 = positions.posmat.on_bt_off_bt{n,nty}; 
    
    m = allcombstable(k,4); % data 2
    mty = allcombstable(k,5); % ON/OFF type 2
    ctype_sub2 = allcombstable(k,6); % cell type 
    
    pos2 = positions.posmat.on_bt_off_bt{m,mty}; 
    
    store_pos_1{k} = pos1; 
    store_pos_2{k} = pos2; 
end




%% 2.  Calculate energy for all combinations 


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
buf2use = 1.2;
Rb = 21; % number of bins+1 for radial energy    


[elin_bin_gather,elin_median_gather,elin_std_gather] = deal(zeros(length(allcombstable),Rb-1)); 


parfor k=1:length(allcombstable)
    % need to register COMs for mixed mosaics 
    if ~ismember(k,locate_ONOFF_samerec)
        pos_on = store_pos_1{k};
        pos_off = store_pos_2{k}; 
        dist = centroid(pos_on) - centroid(pos_off);
        pos_off = pos_off + dist; 
    else
        pos_on = store_pos_1{k};
        pos_off = store_pos_2{k}; 
    end
    store_pos_on_COMregistered{k} = pos_on;
    store_pos_off_COMregistered{k} = pos_off;

    energy_stat_1  = get_energy_for_PIPP(pos_on,pos_off,'Move','ON','ROI_type',roi_type,'print_indiv_maps',false,...
        'print_average_map',false,'buf',buf2use,'Interaction_type',interaction_type); 
    energy_stat_2  = get_energy_for_PIPP(pos_on,pos_off,'Move','OFF','ROI_type',roi_type,'print_indiv_maps',false,...
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
        
        energymat_gather_12(:,:,nr) = energymat;
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



%% Distribution of area under curve: True vs Fake mosaic pairs 

% Generate true null (mixed pairs) distribution and empirical distribution
% : weighted area under curve 
% ******** Use only ConvexHull data ********

gather_elin_bin_avg = elin_bin_gather; 
gather_elin_median_avg = elin_median_gather; 
gather_elin_std_avg = elin_std_gather; 
dx = diff(gather_elin_bin_avg(1,1:2)); 

empindx = locate_ONOFF_samerec; 
nullindx = setxor(1:size(allcombstable,1),empindx);
gather_elin_median_avg_fit = zeros(size(gather_elin_bin_avg)); 
AreaGather = zeros(1,size(allcombstable,1)); 
parfor k=1:size(allcombstable,1)
    rr = gather_elin_bin_avg(k,:); 
    yy = gather_elin_median_avg(k,:);
    th_rad = linspace(0,3*pi/2,length(rr)); 
    fun = @(x,xdata) x(1)*cos(xdata + x(2)).*exp(-xdata.*x(3)./(2*pi)) + x(4); 
    
    x0 = [-2,0,1,0];
    lb = [-5,0,0.1,-3]; ub=[5,0,3,3];   
    A = []; b = []; Aeq = []; beq = [];
    
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxIterations',1000,'OptimalityTolerance',1e-8,...
        'FunctionTolerance',1e-8,'display','iter'); %'UseParallel',true); 
    x = lsqcurvefit(fun,x0,th_rad,yy,lb,ub,options);
    gather_elin_median_avg_fit(k,:) = fun(x,th_rad);
    
    if x(1)*cos(th_rad(1))>0
        AreaGather(k) = - sum( (1/2).*(abs(yy(2:end))-abs(yy(1:end-1))).*dx + abs(yy(1:end-1)).*dx ); 
    else
        AreaGather(k) = + sum( (1/2).*(abs(yy(2:end))-abs(yy(1:end-1))).*dx + abs(yy(1:end-1)).*dx ); 
    end   
end
AreaEmp = AreaGather(empindx); 
AreaNull = AreaGather(nullindx); 



fittype = 'normal';
nBoot = 1000; % bootstrap with replacement 
alpha_ci = 0.05; % 95% CI 
opt = statset('UseParallel',true);
[bootarea, bootsam] = bootstrp(nBoot,@mean,AreaNull,'Options', opt);
bootbinedge = linspace(min(bootarea),max(bootarea),50); 
db = diff(bootbinedge(1:2)); 
bootbinedge = linspace(bootbinedge(1)-db,bootbinedge(end)+db,50); 
bootbincent = (bootbinedge(1:end-1)+bootbinedge(2:end))./2; 
[bootfreq,~] = histcounts(bootarea,bootbinedge); 
bootfreqnorm = bootfreq./(sum(bootfreq)*db); 
bootCi = std(bootarea)*abs(norminv(alpha_ci/2)); % from one population only because the other is only 4 samples (otherwise use ttest2)  

Xll = max([bootbinedge(end)+db mean(AreaEmp)*1.1]); 

Xfit = linspace(bootbinedge(1)-db,Xll,100);
if strcmpi(fittype,'lognormal')
    pd1 = fitdist(bootarea, 'lognormal');
    Yfit = pdf(pd1,Xfit);
    Yfit_normalized = Yfit.*diff(Xfit(1:2)); 
    norm_mean = exp(pd1.mu + 0.5*pd1.sigma^2); 
    norm_std = sqrt(exp(2*(pd1.mu + pd1.sigma^2)) - exp(2*pd1.mu + pd1.sigma^2)); 
elseif strcmpi(fittype,'normal')
    pd1 = fitdist(bootarea, 'normal');
    Yfit = pdf(pd1,Xfit);
    Yfit_normalized = Yfit.*diff(Xfit(1:2)); 
    norm_mean = pd1.mu; 
    norm_std = pd1.sigma; 
end

bootCi_upper = mean(bootarea) + std(bootarea)*abs(norminv(alpha_ci/2));
bootCi_lower = mean(bootarea) - std(bootarea)*abs(norminv(alpha_ci/2));

%[h,p,ci,stats] = ttest(bootarea, mean(AreaEmp)) 
[h,p,ci,stats] = ttest(AreaNull, mean(AreaEmp)) 


% Compute Cohen's d: Assessing Effect size based on difference of means and size of sample 
cohen_d = (mean(AreaEmp) - mean(bootarea))/sqrt( ((numel(AreaEmp)-1)*std(AreaEmp) +...
    (numel(bootarea)-1)*std(bootarea))/(numel(AreaEmp)+numel(bootarea)));

%% Define colors for RF fill, and cell types whose mosaics will be analyzed 

new_purple = [136/256 46/256 114/256];
new_blue = [25/256 101/256 176/256];
new_orange = [230/256 85/256 24/256];

light_purple = [171/256 122/256 150/256];
light_blue = [32/256 164/256 255/256];
light_orange = [237/256 152/256 108/256];

new_green = [0.13333 0.5451 0.13333]; % forest green

rf_fill_alph = 0.2; 
rf_contour_wid = 3; 
clx_off_rf_edge = [0 0 0]; 
clx_off_com = [0.5 0.5 0.5]; 

matlab_parula1 = [0    0.4470    0.7410];


%% Figure showing null distribution of IMCE for mosaic pairs from different retinas and IMCE for mosaic pairs from same retina 

FS = 33; MS=15; alph = 0.15; 
clx_sim = [0 0 0];
clx_emp = new_blue; 


[~,ind_1] = min(abs(bootCi_upper-Xfit));
ind_2 = length(Xfit); 
fillx_sim = [Xfit(ind_1) Xfit(ind_2) Xfit(ind_2:-1:ind_1) ...
    Xfit(ind_1)];
filly_sim = [0 0 Yfit_normalized(ind_2:-1:ind_1) Yfit_normalized(ind_1)];
xll = [Xfit(1) Xfit(end)]; 
yll = [0 0.055]; % 0.042

barx = bootbincent; dbx = diff(barx(1:2)); rbarx= []; 
bary = bootfreqnorm.*diff(Xfit(1:2)); rbary = []; 
for i=1:length(barx)
    rbarx = [rbarx barx(i)-dbx/2 barx(i)-dbx/2]; 
    rbary = [rbary bary(i) bary(i)]; 
end
rbarx = [rbarx barx(end)+dbx/2 barx(end)+dbx/2]; 
rbary = [0 rbary 0]; 

hf100 = figure(100); clf(hf100); hold on;
set(hf100,'position',[473   391   434   362]); 
plot(rbarx, rbary, '-','color',clx_sim,'linewidth',1); 
plot(Xfit, Yfit_normalized, '-','linewidth',2,'color',[clx_sim 0.5]); 
h1 = fill(fillx_sim,filly_sim,clx_sim); 
set(h1,'facealpha',0.3,'edgecolor','none'); 
plot(mean(AreaEmp),zeros(1,1),'o','markersize',MS,'markeredgecolor','none','markerfacecolor',clx_emp); 
set(gca,'xlim',xll,'ylim',yll,'fontsize',FS,'linewidth',2,'tickdir','out'); 
xlabel('Coordination index')
ylabel('Pdf');


    



