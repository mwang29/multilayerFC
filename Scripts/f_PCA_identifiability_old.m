function [Idiff_orig,Ident_mat_orig,Idiff_recon,Idiff_opt,recon_matrix_opt,Ident_mat_recon_opt,PCA_comps_range,m_star, latent] = f_PCA_identifiability_old(orig_matrix,Test_index,Retest_index,configs)

%% Compute Identifiability matrix, original FCs
orig_matrix_test = orig_matrix(:,Test_index);
orig_matrix_retest = orig_matrix(:,Retest_index);   
Ident_mat_orig = corr(orig_matrix_retest,orig_matrix_test); 
%% Idiff computation, original FCs
%% Idiff = <Iself> - <Iothers> (Amico & Goñi, SciRep 2018; i.e. <Main diagonal> - <Triu> of the Identifiability matrix)     
mask_diag = logical(eye(size(Ident_mat_orig)));
Iself_orig = mean(Ident_mat_orig(mask_diag));
Iothers_orig = mean(Ident_mat_orig(~mask_diag));
Idiff_orig = (Iself_orig - Iothers_orig) * 100; % Identifiability score for original connectivity data       

Idiff_recon = zeros(1,configs.max_numPCs); %Initialization for Identifiability score for reconstructed connectivty data (indexed by number of components).


%% Differential Identifiability (Idiff) evaluation of PCA decomposition into FC-modes (see also Fig.2 of Amico & Goñi, SciRep 2018)
PCA_comps_range = 2:configs.max_numPCs; % Note that at least 2 PCA components are needed for Idiff evaluation
disp('FC reconstruction with:')
[FC_modes, projected_FC_modes, latent] = pca(orig_matrix,'NumComponents',configs.max_numPCs);     
for i= PCA_comps_range
    recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)'; % PCA reconstructed demeaned data    
    recon_matrix = bsxfun(@plus, recon_matrix,mean(orig_matrix)); % plug the mean back
    recon_matrix_test = recon_matrix(:,Test_index);
    recon_matrix_retest = recon_matrix(:,Retest_index);
    %% Compute Identifiability matrix, reconstructed FCs
    Ident_mat_recon = corr(recon_matrix_retest,recon_matrix_test);
    %% Idiff computation, reconstructed FCs
    %% Idiff = <Iself> - <Iothers> (Amico & Goñi, SciRep 2018; i.e. <Main diagonal> - <Triu> of the Identifiability matrix)       
    Iself_recon = mean(Ident_mat_recon(mask_diag));
    Iothers_recon = mean(Ident_mat_recon(~mask_diag));
    Idiff_recon(i) = (Iself_recon - Iothers_recon) * 100;
end

%% Identifiability matrix at optimal reconstruction (see also Fig 3C1 in Amico & Goñi, SciRep 2018)
Idiff_opt = max(Idiff_recon);
m_star = PCA_comps_range(Idiff_recon(2:end)==Idiff_opt); % obtain the optimal number of PCA components for reconstruction. First point is intentionally skipped since range of components starts on 2.
recon_matrix_opt = projected_FC_modes(:,1:m_star) * FC_modes(:,1:m_star)'; % PCA reconstructed demeaned data    
recon_matrix_opt = bsxfun(@plus, recon_matrix_opt,mean(orig_matrix)); % plug the mean back
recon_matrix_opt_test = recon_matrix_opt(:,Test_index);
recon_matrix_opt_retest = recon_matrix_opt(:,Retest_index);
%% Compute Recon Identifiability matrix at optimal point
Ident_mat_recon_opt = corr(recon_matrix_opt_retest,recon_matrix_opt_test);

