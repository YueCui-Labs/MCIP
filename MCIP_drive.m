function MCIP_drive(hemi, function_ts_list_file, ana_conn_file ,ref_atlas_file, surf_file, neighbor_file, topo_file, label_w, smooth_w, out_dir, method)
%MCIP drive program of multi-modality atlas individualization
% Li, Chengyi, 2022.1.18


% read atlas
atlas_st = gifti(ref_atlas_file);
all_vert=atlas_st.cdata;
all_vert(all_vert<0)=0;
atlas_mask = all_vert>0;
ref_atlas = all_vert(atlas_mask);

lookup = unique(ref_atlas); % 105
uni_atlas = ref_atlas;
for i=1:length(lookup)
    uni_atlas(ref_atlas==lookup(i))=i;
end
roi_num = sum(atlas_mask);

if hemi == 'L'
    he=1;
else
    he=2;
end
    
% read function time series file list
function_ts_list = textscan(fopen(function_ts_list_file), '%s');
function_ts_list =function_ts_list{1};
func_num = length(function_ts_list);

function_ts=[];
for ff=1:func_num
    function_ts_st=cifti_read(function_ts_list{ff});
    cdata=function_ts_st.cdata; % 91282 x ts
    start=function_ts_st.diminfo{1,1}.models{1,he}.start;
    count=function_ts_st.diminfo{1,1}.models{1,he}.count;
    vertlist=function_ts_st.diminfo{1,1}.models{1,he}.vertlist;
    
    all_vert_ts = zeros(size(all_vert,1), size(cdata,2));
    tmpts = cdata(start:start+count-1,:);
    tmpts = (tmpts - mean(tmpts,2)) ./ (std(tmpts,[],2) + eps);
    all_vert_ts(1+vertlist,:) = tmpts;
    
    function_ts = [function_ts, all_vert_ts(atlas_mask, :)];
end

%% TODO: stability test
%function_ts = function_ts(:,1:240);

% load anatomical connection matrix
if ana_conn_file(end-2:end) == 'dot'
    conn_sp = spconvert(load(ana_conn_file));
    ana_conn = conn_sp; % #MESH x #volume
    
    save([ana_conn_file(1:end-3),'mat'],'conn_sp');
    delete(ana_conn_file);

elseif ana_conn_file(end-2:end) == 'mat'
    ana_conn = load(ana_conn_file);
    ana_conn = ana_conn.conn_sp;
end
ana_conn = full(ana_conn);
ana_conn = ana_conn(atlas_mask, :); % 29696/29716 x #MESH

% load topo distance matrix
topo_dis = gifti(topo_file);
topo_dis = topo_dis.cdata';
topo_dis = single(topo_dis(:, atlas_mask));  % #parc x #MESH

% read surface
surf = gifti(surf_file);
face=surf.faces;

% construct neighbor matrix
if ~exist(neighbor_file)
    neighbor_mat = zeros(size(all_vert,1), size(all_vert,1));
    for e=1:size(face,1)
        for vi=1:3
            for vj=vi+1:3
                neighbor_mat(face(e,vi),face(e,vj))=1; %vi - vj
                neighbor_mat(face(e,vj),face(e,vi))=1; %vj - vi
            end
        end
    end
    neighbor_mat = neighbor_mat(atlas_mask,atlas_mask);

    neighbor_mat = sparse(neighbor_mat);
    save(neighbor_file,'neighbor_mat');
else
    neighbor_mat = load(neighbor_file);
    neighbor_mat = neighbor_mat.neighbor_mat;
end

%% MCIP
disp('MCIP start');
[individual_atlas, lambdas] = MCIP(function_ts, ana_conn, topo_dis, label_w, smooth_w, uni_atlas, neighbor_mat);

save(fullfile(out_dir,[method, '_', hemi, '_lambdas.mat']), 'lambdas');

%% save func file
for i=1:length(lookup)
    uni_atlas(individual_atlas==i)=lookup(i);
end

all_vert(atlas_mask) = uni_atlas;
sv=gifti(all_vert);

func_gii_file=fullfile(out_dir,[method, '_', hemi, '.32k_fs_LR.func.gii']);
save(sv,func_gii_file);

end

