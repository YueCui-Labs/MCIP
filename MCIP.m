function [indiv_atlas, lambdas] = MCIP(function_ts,ana_conn, topo_cost, label_weight, smooth_weight, ref_atlas, neighbor_mat)
% MCIP: Multimodal Connectivity-based Individualized Parcellation
% Li, Chengyi, 2022.1.18

% ----------inputs-----------
% function_ts_cell: functional time series cell, #MESH x ts (elem matrix)
% ana_conn: anatomical connection matrix, #MESH x #tareget
% topo_cost: #parc_num x #MESH
% label_weight: how much you trust reference atlas
% smooth_weight
% ref_atlas: reference atlas label, #MESH x 1
% neighbor_mat: search light, only find NN in this neighbor, #MESH x #MESH
% ----------outputs-----------
% indiv_atlas: individual atlas label, #MESH x 1


[MESH,time_point] = size(function_ts);
profile_num = size(ana_conn,2);

parc_num = max(ref_atlas);

% log transform
ana_conn = log2(ana_conn+1); % #MESH x #wholebrain

% prepare normal distribution resampling
%resam = randn(1,MESH);
%resam = (resam - min(resam)) / ( max(resam) - min(resam)); % (0,1)
%resam_ordered = sort(resam);
quantile = linspace(0,1,MESH+2);
resam_ordered = norminv(quantile(2:end-1), 0.5, 0.125);
%resam_ordered = (resam_ordered - min(resam_ordered)) / ( max(resam_ordered) - min(resam_ordered)); % (0,1)

% reliability
mid_point = round(time_point / 2);
function_ts_1 = function_ts(:,1:mid_point);
function_ts_2 = function_ts(:,mid_point+1:end);

%% iterative clustering

ind_new=ref_atlas; ind_old=zeros(MESH,1);

%size(ind_new)
%size(ind_old)
t=0;

t1 = []; t2 = [];

lambdas = [];
while sum(ind_new ~= ind_old) / MESH > 0.005

    tic; % t1

    ind_old = ind_new;
    
    %% func. conn.
    mfunc_ts = zeros(parc_num,time_point);
    for p=1:parc_num
        mfunc_ts(p,:) = mean(function_ts(ind_new==p,:), 1);
    end
    
%% reliability
    mfunc_ts_1 = mfunc_ts(:,1:mid_point); 
    mfunc_ts_2 = mfunc_ts(:,mid_point+1:end); 
    func_conn_1 = corr(function_ts_1', mfunc_ts_1'); % #tp x #MESH , #tp x #parc --> #MESH x #parc
    func_conn_2 = corr(function_ts_2', mfunc_ts_2'); % #tp x #MESH , #tp x #parc --> #MESH x #parc
    
    reli= zeros(MESH, 1);
    for i=1:MESH
        reli(i) = corr(func_conn_1(i,:)', func_conn_2(i,:)');
    end
    % reli = rescale(reli, scale, 1-scale);
    lambdas = [lambdas, reli];
    reli(reli<0) = 0;
    
%%
%    func_dist = 1 - (func_conn_1 + func_conn_2)/2;
     func_dist = pdist2(function_ts, mfunc_ts, 'correlation');

    % distribution resample: #MESH x #parc
    for p=1:parc_num
        func_dist_p = func_dist(:,p);
        [~,func_dist_ind] = sort(func_dist_p);
        func_dist(func_dist_ind, p) = resam_ordered;
    end

    %% anat. conn.
    mana_conn = zeros(parc_num,profile_num);
    for p=1:parc_num
        mana_conn(p,:) = mean(ana_conn(ind_new==p,:), 1);
    end

    ana_dist = pdist2(ana_conn, mana_conn, 'correlation');

    for p=1:parc_num
        ana_dist_p = ana_dist(:,p);
        [~, ana_dist_ind] = sort(ana_dist_p);
        ana_dist(ana_dist_ind, p) = resam_ordered;
    end

    
    %% fusing
    mid_dist = reli .* func_dist + (1- reli) .* ana_dist;

    t1(end+1) = toc;  % t1

    tic; % t2
    %% gco algorithm
    h = GCO_Create(MESH,parc_num);

%    smaller_weight = min(label_weight, smooth_weight);
    scaling = 100;
    data_cost = (mid_dist'  + label_weight * topo_cost) * scaling;
    data_cost = data_cost - min(data_cost(:));
    GCO_SetDataCost(h, double(data_cost));
    
    GCO_SetNeighbors(h, double(neighbor_mat));

    Smooth_cost=ones(parc_num, parc_num);
    Smooth_cost(1:(parc_num+1):end)=0;% sets the diagonal to 0
    GCO_SetSmoothCost(h, double(smooth_weight * Smooth_cost * scaling)); % increase the smooth cost, to get smoother results
    
    GCO_Expansion(h);% maximum number of iterations should be a finite number
    
    ind_new=GCO_GetLabeling(h);

    GCO_Delete(h);
    t2(end+1) = toc; % t2


    % stop
    t=t+1;
    disp(t)
    if t>10
%        disp('ITERATION NOT CONVERAGE!!');
        break
    end

end
indiv_atlas = ind_new;


%mean(t1)
%std(t1)

%mean(t2)
%std(t2)
%
mean(t1+t2)
std(t1+t2)

end
