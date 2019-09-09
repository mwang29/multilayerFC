%% initialize environment
close all;
clearvars
clc

%% changing default fontsize
fontsize = 20;
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',fontsize-2);

set(0,'DefaultTextFontname','Times New Roman');
set(0,'DefaultTextFontSize',fontsize);

%% configs

configs.numRegions = 360; %Number of brain regions
configs.mask_ut = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask
configs.numDiffMetrics = 10; %
configs.numFCs = 420;
configs.numSubjects = configs.numFCs / configs.numDiffMetrics;
configs.numEdges = nnz(configs.mask_ut);
% configs.numVisits = 2; % 2 visits per subject (test-retest)
configs.max_numPCs = 2*configs.numSubjects; % maximum number of PCs == data dimension
configs.score = 'frobenius';
configs.non_negative = false;
%% load data

connectivity_matrix = zeros(configs.numEdges, 420);
cd('../Connectivity_data')
D = dir;
col = 1;
for k = 4:length(D)
    if D(k).isdir
        currD = D(k).name;
        cd(currD);
        files = dir('*_mean.csv');
        for file = files'
            temp_mat = csvread(file.name);
            upper_triu = temp_mat(configs.mask_ut);
            original_matrix(:,col) = upper_triu;
            col = col + 1;
        end
        cd('..')
    end
end
cd('../Scripts')


%% Sorted list of MnF group avg edges across subjects
mnf_connectivity = original_matrix(:,6:10:end);
group_avg = mean(mnf_connectivity,2);
[sorted, index] = sort(group_avg,'descend');


%% Pairwise identifiability for diffusion metrics
metrics = {'Da', 'Dr', 'Fa', 'Md', 'Mll', 'Mnf', 'Mw', 'OD', 'Po', 'Vic'};
pairwise_mat = nan(configs.numEdges, configs.numSubjects * 2);
test_index = 1:42;
retest_index = 43:84;
idiff_mat = nan(10, 10, 10);
perms = nchoosek(1:10, 2);%permutations of metrics
thresholds = 0.05:0.05:0.5; 


for k = 1:length(perms) %Indexing throughout meshgrid (10 2) metrics
    for thresh_index = 1:10 %index for 10 different thresholds, 10% to 100%
        threshold = thresholds(thresh_index);
        to_keep = floor(threshold*configs.numEdges); %number of edges to keep
        mask_mnf = index(1:to_keep); %indices of kept edges
        connectivity_matrix = original_matrix(mask_mnf,:); %apply mask to orig matrix
        metric1_index = perms(k,1):10:configs.numFCs;
        metric2_index = perms(k,2):10:configs.numFCs;
        pairwise_mat = connectivity_matrix(:,[metric1_index, metric2_index]); %get metrics of interest
        
        fprintf('Computing Idiff for metric %d and %d at threshold %.2f\n', perms(k,1), perms(k,2), threshold)
        metric_1_mat = pairwise_mat(:,test_index); %"test"
        metric_2_mat = pairwise_mat(:,retest_index); %"retest"
        Ident_mat = pdist2(metric_1_mat',metric_2_mat','spearman');
        mask_diag = logical(eye(size(Ident_mat)));
        Idiff = -mean((Ident_mat(mask_diag) - mean([mean(Ident_mat,2) mean(Ident_mat)'],2) ) ./ ((std(Ident_mat,0,2) + std(Ident_mat,0)')./2) ); %compute idiff
        idiff_mat(perms(k,1), perms(k,2), thresh_index) = Idiff;
    end
end

%% Plot idiff curves
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(thresholds, reshape(idiff_mat(1,6,:),10,1))
hold on
plot(thresholds, reshape(idiff_mat(1,4,:),10,1))
plot(thresholds, reshape(idiff_mat(3,6,:),10,1))
legend('Da-Mnf', 'Da-Md', 'Fa-Mnf', 'Location', 'southeast')
title('Idiff thresholded based on group average')
xlabel('Threshold')
ylabel('Idiff')
saveas(fig, '../Images/group_avg_threshold.png')
%% Plot 10x10 max I_diff
fig = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(idiff_mat(:,:,1))
axis square
title('Max Identifiability across 42 subjects')
set(gca,'xtick',1:10,'xticklabel',metrics)
set(gca,'ytick',1:10,'yticklabel',metrics)
colorbar
saveas(fig, '../Images/ident_mat_5pct_threshold.png')

%% Hierarchical clustering by sorting pairs using idiff 
idiff_slice = idiff_mat(:,:,1); %slice at defined threshold
sorted_iDiff = sort(idiff_slice(:)); %sort idiff at slice
for place = 1:length(perms) %get indices based on decreasing idiff
    [index_iDiff(place,1), index_iDiff(place,2)] = find(idiff_slice == sorted_iDiff(place));
end

%% TODO: Calculate aggregate graph entropy (individual subject or all 42 subjects?)


%% TODO: Combine layers based on hierachical clustering (max idiff)

%% How many shared edges in all 42 FCs?
metric = 'Dr';
indices = find(ismember(metrics, metric)):10:configs.numFCs; %index of conn. matrix for given metric
matching_edges = find(connectivity_matrix(:,indices(1)));

for i = indices(2:end)
    next_edges = find(connectivity_matrix(:,i));
    matching_edges = intersect(matching_edges, next_edges);
end

avg_num_edges = nnz(connectivity_matrix(:,indices))/configs.numSubjects;
fprintf('Average edges per subject for %s: %.2f\n', metric, avg_num_edges)
fprintf('Total matching edges across all %d subjects for %s: %d\n', configs.numSubjects, metric, length(matching_edges))
fprintf('%.2f%% of edges are common across subjects\n', 100*length(matching_edges)/avg_num_edges)








