%% batch example of the PCA-decomposition method, proposed in (Amico & Goñi, Nature Scientific Reports 2018),
%% to maximize identifiability in human functional connectomes.
%% 
%% Enrico Amico & Joaquín Goñi, Purdue University
%% version 1.1. May 10, 2018
%
%% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:
%% 1) Enrico Amico and Joaquín Goñi. The quest for identifiability in human functional connectomes. Scientific
%% Reports, (in press), 2018. 

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

%% Load FC (functional connectivity) resting-state data. This is a sample data.
%  dimensions are brain regions x brain regions x subjects
%  The FCs of this sample are organized in this order: Subj1_test, Subj1_retest....SubjN_test, SubjN_retest
% (see Amico & Goñi, SciRep 2018 for details)

load FC_test_retest.mat;

%% Configuration
configs.numRegions = 164; % aparc2009 parcellation. Number of brain regions
configs.mask_ut = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask
configs.numFCs = size(FC,3); %number of connectomes (2 FCs per subject in this case)
configs.numEdges = nnz(configs.mask_ut);
configs.numVisits = 2; % 2 visits per subject (test-retest)
configs.max_numPCs = configs.numFCs; % maximum number of PCs == data dimension

%% Index vectors for test and retest (visit 1 and visit 2; first half and second half)
Test_index = 1:configs.numVisits:configs.numFCs; % change this and next line if you have a different FCs ordering
Retest_index = 2:configs.numVisits:configs.numFCs;


%% create 2D input matrix from original FCs (n FC edges x nSubj, see also Fig.1 of Amico & Goñi, SciRep 2018)
orig_matrix = zeros(nnz(configs.mask_ut),configs.numFCs);
disp('Creating orig_matrix...')
for i=1:configs.numFCs
    aux = FC(:,:,i);
    orig_matrix(:,i) = aux(configs.mask_ut);
end
disp('Done.');

[Idiff_orig,Ident_mat_orig,Idiff_recon,Idiff_opt,recon_matrix_opt,Ident_mat_recon_opt,PCA_comps_range,m_star, latent]  = f_PCA_identifiability_old(orig_matrix,Test_index,Retest_index,configs);
r2 = cumsum(latent) ./ sum(latent);
figure,
plot(r2)

figure('units','normalized','outerposition',[0 0 1 1]), 
%suptitle(['Opt Identifiability Recon at ' int2str(m_star) ' PCA comps'])
subplot(1,3,1);
imagesc(Ident_mat_orig); 
axis square; xlabel('subjects (Test)');ylabel('subjects (Retest)'); 
title(sprintf('original data, Idiff = %0.2f',Idiff_orig)); colorbar; caxis([0.2 1]);
set(gca,'XTick',[],'YTick',[]);

subplot(1,3,2);
plot(PCA_comps_range,Idiff_orig*ones(length(PCA_comps_range),1),'--r','LineWidth',2); hold on;
plot(PCA_comps_range,Idiff_recon(2:end),'-ob','LineWidth',2,'MarkerFaceColor','b','MarkerSize',8); 
plot(m_star,Idiff_opt,'-sk','LineWidth',2,'MarkerFaceColor','k','MarkerSize',12); 

xlabel('number of PCA components'); ylabel('Idiff (%)');
axis square;
legend({'original data','reconstruction','optimal reconstruction'},'Location','SouthEast');
title('Idiff assessment based on PCA decomposition')
sprintf('Maximal Identifiability (%0.2f) found at %d PCA comps',Idiff_opt,m_star)

subplot(1,3,3); 
imagesc(Ident_mat_recon_opt); axis square; xlabel('Subjects Test');ylabel('Subjects Retest'); 
title(sprintf('optimal reconstruction, Idiff = %0.2f',max(Idiff_recon))); colorbar; caxis([0.2 1]);
set(gca,'XTick',[],'YTick',[]);


