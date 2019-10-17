%% initialize environment
close all;
clearvars
clc

%% configs

configs.numRegions = 360; %Number of brain regions
configs.mask_ut = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask
configs.numDiffMetrics = 10; %
configs.numFCs = 420;
configs.numSubjects = configs.numFCs / configs.numDiffMetrics;
configs.numEdges = nnz(configs.mask_ut);
configs.max_numPCs = 2*configs.numSubjects; % maximum number of PCs == data dimension
configs.score = 'frobenius';
configs.non_negative = false;

%% changing default fontsize
fontsize = 20;
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',fontsize-2);

set(0,'DefaultTextFontname','Times New Roman');
set(0,'DefaultTextFontSize',fontsize);

%% load data

connectivity_matrix = zeros(configs.numRegions, configs.numRegions, 420);
cd('../Connectivity_data')
D = dir;
i = 1;
for k = 4:length(D)
    if D(k).isdir
        currD = D(k).name;
        cd(currD);
        files = dir('*_mean.csv');
        for file = files'
            if i == 289 || i == 301
                fprintf('%s\n', file.name)
            end
            connectivity_matrix(:,:,i) = csvread(file.name); %indexed by subject then by metric
            i = i + 1;
        end
        cd('..')
    end
end
cd('../Scripts')
fprintf('Finished Loading Files\n')
%% Check connectedness of SC metrics
fully_connected = 0;
non_connected = [];
nc_counter = 1;
nc_index = [];
for i = 1:420
    if i ~= 289 %PAT_s9310_Po_mean.csv is not symmetric
        [comps, comp_sizes] = get_components(connectivity_matrix(:,:,i)); %get comps
    end
    
    if length(comp_sizes) == 1 %if fully connected
        fully_connected = fully_connected + 1;
    else %if not, keep track of which ones are not
        non_connected(:,:,nc_counter) = connectivity_matrix(:,:,i);
        nc_index = [nc_index, i];
        nc_counter = nc_counter + 1;
    end
end


% Subject 9 does not have fully connected SC matrices. Also,
% PAT_s9312_Da_mean is not fully connected but other metrics are fully
% connected. Nnz of 5667 while other metrics have 6392. (corrupted data?) 

%% Communicability (test one one matrix)

A = connectivity_matrix(:,:,1);
D = diag(sum(A));
beta = 0.1;
communicability = exp(beta*D^(-1/2)*A*D^(-1/2));

%Beta 0.1 scales our connectivity values down to more manageable numbers

%% Communicability distance
A_ii = repmat(diag(communicability), 1, configs.numRegions);
A_jj = A_ii';
distance = A_ii + A_jj - (2*communicability);
distance(distance < 0 ) = 0;
distance = sqrt(distance);

%% Compile all communicability into edges x subject form
beta = 0.1;
communicability_matrix = NaN(configs.numRegions^2, configs.numFCs);
for sc = 1:configs.numFCs
    temp_A = connectivity_matrix(:,:,sc);
    temp_A(temp_A == 0) = eps;
    temp_D = diag(sum(temp_A));
    temp_comm = exp(beta*temp_D^(-1/2)*temp_A*temp_D^(-1/2));
    communicability_matrix(:,sc) = reshape(temp_comm, [configs.numRegions^2,1]);
end

%% Sorted list of MnF group avg edges across subjects
mnf_communicability = communicability_matrix(:,6:10:end);
group_avg = mean(mnf_communicability,2);
[sorted, index] = sort(group_avg,'descend');

%% Pairwise identifiability for diffusion metrics
metrics = {'Da', 'Dr', 'Fa', 'Md', 'Mll', 'Mnf', 'Mw', 'OD', 'Po', 'Vic'};
pairwise_mat = nan(configs.numEdges, configs.numSubjects * 2);
test_index = 1:42;
retest_index = 43:84;
idiff_mat = nan(10, 10, 10);
perms = nchoosek(1:10, 2); %permutations of metrics
thresholds = 1; 

for k = 1:length(perms) %Indexing throughout meshgrid (10 2) metrics
    for thresh_index = 1 %index for 10 different thresholds, 10% to 100%
        threshold = thresholds(thresh_index);
        to_keep = floor(threshold*configs.numEdges); %number of edges to keep
        mask_mnf = index(1:to_keep); %indices of kept edges
        connectivity_matrix = communicability_matrix(mask_mnf,:); %apply mask to orig matrix
        metric1_index = perms(k,1):10:configs.numFCs;
        metric2_index = perms(k,2):10:configs.numFCs;
        pairwise_mat = connectivity_matrix(:,[metric1_index, metric2_index]); %get metrics of interest
        
        fprintf('Computing Idiff for metric %d and %d at threshold %.2f\n', perms(k,1), perms(k,2), threshold)
        metric_1_mat = pairwise_mat(:,test_index); %"test"
        metric_2_mat = pairwise_mat(:,retest_index); %"retest"
        Ident_mat = pdist2(metric_1_mat',metric_2_mat','correlation');
        mask_diag = logical(eye(size(Ident_mat)));
        Idiff = -mean(Ident_mat(mask_diag)) + mean(Ident_mat(~mask_diag));
%       Idiff = -mean((Ident_mat(mask_diag) - mean([mean(Ident_mat,2) mean(Ident_mat)'],2) ) ./ ((std(Ident_mat,0,2) + std(Ident_mat,0)')./2) ); %compute idiff
        idiff_mat(perms(k,1), perms(k,2), thresh_index) = Idiff;

    end
end

%% Plotting
fig = imagesc(idiff_mat(:,:,1));
axis square
title('Identifiability across 42 subjects, all edges, pearson corr')
set(gca,'xtick',1:10,'xticklabel',metrics)
set(gca,'ytick',1:10,'yticklabel',metrics)
colorbar
saveas(fig, '../Images/ident_mat_100pct_communicability_corr.fig')


















