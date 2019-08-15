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
cd('Connectivity_data')
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
            connectivity_matrix(:,col) = upper_triu;
            col = col + 1;
        end
        cd('..')
    end
end
cd('..')

%% Distance matrix for each patient (n_edges x n_metrics x n_subjects)
metrics = {'Da', 'Dr', 'Fa', 'Md', 'Mll', 'Mnf', 'Mw', 'OD', 'Po', 'Vic'};
ut_mask = triu(true(configs.numDiffMetrics),1);
diag_mask = logical(eye(configs.numDiffMetrics));
cols = 1:10:configs.numFCs;
for subject = 1:configs.numSubjects
    range = cols(subject):cols(subject)+ 9;
    subject_mat = zscore(connectivity_matrix(:,range));
    dist_mat(:,:,subject) = pdist2(subject_mat',subject_mat', 'euclidean');
    subplot(6,7,subject)
    imagesc(dist_mat(:,:,subject))
    axis square
    colorbar
end
dist_mat_avg = squeeze(mean(dist_mat,3));
fig = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(dist_mat_avg)
axis square
title('Average of SC metric distances across 42 subjects')
set(gca,'xtick',[1:10],'xticklabel',metrics)
set(gca,'ytick',[1:10],'yticklabel',metrics)
colorbar
saveas(fig, 'dist_mat_avg.png')

%% Pairwise identifiability for closely related metrics
% pairs = [1, 4; 2, 4; 3, 10; 9, 10];
pairs = nchoosek(1:10,2);
metrics = {'Da', 'Dr', 'Fa', 'Md', 'Mll', 'Mnf', 'Mw', 'OD', 'Po', 'Vic'};
pairwise_mat = nan(configs.numEdges, configs.numSubjects * 2);
test_index = 1:42;
retest_index = 43:84;
max_idiff_mat = nan(10);
for metric_1 = 1:10
    for metric_2 = 1:10
        metric1_index = metric_1:10:configs.numFCs;
        metric2_index = metric_2:10:configs.numFCs;
        pairwise_mat = connectivity_matrix(:,[metric1_index, metric2_index]);
        %pairwise_mat(pairwise_mat~=0) = 1; %binarize the matrix before ident
        configs.binary = 0;
        for type = {'frobenius'}
            configs.score = type{1};
            fprintf('Computing for metric %d and %d\n', metric_1, metric_2)
            [Idiff_orig,Ident_mat_orig,Idiff_recon,Idiff_opt,recon_matrix_opt,Ident_mat_recon_opt,PCA_comps_range,m_star, latent] = f_PCA_identifiability(pairwise_mat, test_index, retest_index, configs);
            max_idiff_mat(metric_1, metric_2) = max(Idiff_recon);
            r2 = cumsum(latent) ./ sum(latent);
    %         figure,
    %         plot(r2)
    %         fig = figure('units','normalized','outerposition',[0 0 1 1]); % Plot Figures
    %         %suptitle(['Opt Identifiability Recon at ' int2str(m_star) ' PCA comps'])
    %         subplot(1,3,1);
    %         imagesc(Ident_mat_orig); 
    %         axis square; xlabel('subjects (Test)');ylabel('subjects (Retest)'); 
    %         title(sprintf('original data, Idiff = %0.2f',Idiff_orig)); colorbar; % caxis([0.2 1]);
    %         set(gca,'XTick',[],'YTick',[]);
    % 
    %         subplot(1,3,2);
    %         plot(PCA_comps_range,Idiff_orig*ones(length(PCA_comps_range),1),'--r','LineWidth',2); hold on;
    %         plot(PCA_comps_range,Idiff_recon,'-ob','LineWidth',2,'MarkerFaceColor','b','MarkerSize',8); 
    %         plot(m_star,Idiff_opt,'-sk','LineWidth',2,'MarkerFaceColor','k','MarkerSize',12); 
    % 
    %         xlabel('number of PCA components'); ylabel('Idiff (%)');
    %         axis square;
    %         legend({'original data','reconstruction','optimal reconstruction'},'Location','SouthEast');
    %         title(sprintf('%s | %s | %s', metrics{temp_pair(1)}, metrics{temp_pair(2)}, type{1}))
    %         
    %         sprintf('Maximal Identifiability (%0.2f) found at %d PCA comps',Idiff_opt,m_star)
    % 
    %         subplot(1,3,3); 
    %         imagesc(Ident_mat_recon_opt); axis square; xlabel('Subjects Test');ylabel('Subjects Retest'); 
    %         title(sprintf('optimal reconstruction, Idiff = %0.2f',max(Idiff_recon))); colorbar; % caxis([0.2 1]);
    %         set(gca,'XTick',[],'YTick',[]);
    %         saveas(fig, sprintf('Images/1-6_binarized_%s.png', type{1}))
        end

    end
end
%% Save figure
fig = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(max_idiff_mat)
axis square
title('Max Identifiability across 42 subjects')
set(gca,'xtick',[1:10],'xticklabel',metrics)
set(gca,'ytick',[1:10],'yticklabel',metrics)
colorbar
saveas(fig, 'Images/dist_mat_avg.png')
%% Shared edges in all 42 FCs
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








