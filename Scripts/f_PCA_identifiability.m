function [Idiff_orig,Ident_mat_orig,Idiff_recon,Idiff_opt,recon_matrix_opt,Ident_mat_recon_opt,PCA_comps_range,m_star,latent,num_neg] = ...
    f_PCA_identifiability(orig_matrix,Test_index,Retest_index,configs)

%% Compute Identifiability matrix, original FCs
orig_matrix_test = orig_matrix(:,Test_index);
orig_matrix_retest = orig_matrix(:,Retest_index);

%% Compute Distance/Similarity Matrix for pcacov
dist_mat = pdist2(orig_matrix_test',orig_matrix_retest', 'euclidean');
similarity_mat = ones(42,42) - dist_mat;

%% Plot Distance/Similarity Matrix
% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% imagesc(similarity_mat)
% axis square
% title('Similarity matrix (subject x subject) between metrics')
% xlabel('Subject'), ylabel('Subject')
% colorbar
% saveas(fig, '../Images/similarity_between_metrics.png')

%% Non-negative and binarization
if configs.non_negative
    num_neg = nan(configs.max_numPCs-1);
end

if configs.binary == 1
    orig_matrix_test(orig_matrix_test~=0) = 1;
    orig_matrix_retest(orig_matrix_retest~=0) = 1;
end
%% Compute Identifiability matrix, original FCs
if strcmp(configs.score,'correlation')     
    Ident_mat_orig = corr(orig_matrix_retest,orig_matrix_test); 
elseif strcmp(configs.score,'frobenius')
    Ident_mat_orig = pdist2(orig_matrix_test',orig_matrix_retest','euclidean')./configs.numEdges;
elseif strcmp(configs.score,'geodesic')
    Ident_mat_orig = f_GD_group(orig_matrix_test,orig_matrix_retest,configs.numRegions);
elseif strcmp(configs.score,'minkowski')
    Ident_mat_orig = pdist2(orig_matrix_test',orig_matrix_retest',configs.score,configs.DistParam);
elseif strcmp(configs.score,'corrDist')
    Ident_mat_orig = pdist2(orig_matrix_test',orig_matrix_retest','correlation');
elseif strcmp(configs.score,'hamming')
    Ident_mat_orig = pdist2(orig_matrix_test',orig_matrix_retest','hamming');
elseif strcmp(configs.score,'jaccard')
    Ident_mat_orig = pdist2(orig_matrix_test',orig_matrix_retest','jaccard');
elseif strcmp(configs.score,'spearman')
    Ident_mat_orig = pdist2(orig_matrix_test',orig_matrix_retest','spearman');
end
        
%% Idiff computation, original FCs
%% Idiff = <Iself> - <Iothers> (Amico & Goñi, SciRep 2018; i.e. <Main diagonal> - <Triu> of the Identifiability matrix)     
mask_diag = logical(eye(size(Ident_mat_orig)));
if strcmp(configs.score,'correlation')  
    Idiff_orig = mean((Ident_mat_orig(mask_diag) - mean([mean(Ident_mat_orig,2) mean(Ident_mat_orig)'],2) ) ./ ((std(Ident_mat_orig,0,2) + std(Ident_mat_orig,0)')./2) );
else %if strcmp(configs.score,'frobenius') || strcmp(configs.score,'geodesic') || strcmp(configs.score,'minkowski')
    Idiff_orig = -mean((Ident_mat_orig(mask_diag) - mean([mean(Ident_mat_orig,2) mean(Ident_mat_orig)'],2) ) ./ ((std(Ident_mat_orig,0,2) + std(Ident_mat_orig,0)')./2) );
end

%% Differential Identifiability (Idiff) evaluation of PCA decomposition into FC-modes (see also Fig.2 of Amico & Goñi, SciRep 2018)
Idiff_recon = zeros(1,configs.max_numPCs); %Initialization for Identifiability score for reconstructed connectivty data (indexed by number of components).
PCA_comps_range = 1:configs.max_numPCs;
disp('FC reconstruction with:')
[COEFF, SCORE, latent] = pcacov(similarity_mat);  
for i = PCA_comps_range
    recon_matrix = SCORE(:,1:i)*COEFF(:,1:i)'; % PCA reconstructed demeaned data    
    recon_matrix = bsxfun(@plus,recon_matrix,mean(orig_matrix)); % plug the mean back
    if configs.non_negative == 1
        num_neg(i) = nnz(recon_matrix<-0.05);
        recon_matrix(recon_matrix<0) = 0;
    end
    recon_matrix_test = recon_matrix(:,Test_index);
    recon_matrix_retest = recon_matrix(:,Retest_index);
    if configs.binary == 1
        recon_matrix_test(recon_matrix_test>0.01|recon_matrix_test<-0.01) = 1;
        recon_matrix_test(recon_matrix_test<=0.01&recon_matrix_test>=-0.01) = 0;
        recon_matrix_retest(recon_matrix_retest>0.01|recon_matrix_retest<-0.01) = 1;
        recon_matrix_retest(recon_matrix_test<=0.01&recon_matrix_retest>=-0.01) = 0;
    end
    %% Compute Identifiability matrix, reconstructed FCs
    if strcmp(configs.score,'correlation')  
        Ident_mat_recon = corr(recon_matrix_retest,recon_matrix_test);
    elseif strcmp(configs.score,'frobenius')
        Ident_mat_recon = pdist2(recon_matrix_test',recon_matrix_retest','euclidean')./configs.numEdges;
    elseif strcmp(configs.score,'geodesic')
        Ident_mat_recon = f_GD_group(recon_matrix_test,recon_matrix_retest,configs.numRegions);
    elseif strcmp(configs.score,'minkowski')
        Ident_mat_recon = pdist2(recon_matrix_test',recon_matrix_retest',configs.score,configs.DistParam);
    elseif strcmp(configs.score,'corrDist')
        Ident_mat_recon = pdist2(recon_matrix_test',recon_matrix_retest','seuclidean');
    elseif strcmp(configs.score,'hamming')
        Ident_mat_recon = pdist2(recon_matrix_test',recon_matrix_retest','hamming');
    elseif strcmp(configs.score,'jaccard')
        Ident_mat_recon = pdist2(recon_matrix_test',recon_matrix_retest','jaccard');
    elseif strcmp(configs.score,'spearman')
        Ident_mat_recon = pdist2(recon_matrix_test',recon_matrix_retest','spearman');
    else        
        error('configs.score %s is not a valid method',configs.score)
    end
    
    %% Idiff computation, reconstructed FCs
    % Idiff = <Iself> - <Iothers> (Amico & Goñi, SciRep 2018; i.e. <Main diagonal> - <Triu> of the Identifiability matrix)       
    if strcmp(configs.score,'correlation')  
        Idiff_recon(i) = mean((Ident_mat_recon(mask_diag) - mean([mean(Ident_mat_recon,2) mean(Ident_mat_recon)'],2) ) ./ ((std(Ident_mat_recon,0,2) + std(Ident_mat_recon,0)')./2) );
    else %if strcmp(configs.score,'frobenius') || strcmp(configs.score,'geodesic') || strcmp(configs.score,'minkowski')
        Idiff_recon(i) = -mean((Ident_mat_recon(mask_diag) - mean([mean(Ident_mat_recon,2) mean(Ident_mat_recon)'],2) ) ./ ((std(Ident_mat_recon,0,2) + std(Ident_mat_recon,0)')./2) );
    end
end

%% Identifiability matrix at optimal reconstruction (see also Fig 3C1 in Amico & Goñi, SciRep 2018)
Idiff_opt = max(Idiff_recon);
m_star = PCA_comps_range(Idiff_recon==Idiff_opt); % obtain the optimal number of PCA components for reconstruction.
recon_matrix_opt = SCORE(:,1:m_star)*COEFF(:,1:m_star)'; % PCA reconstructed demeaned data    
recon_matrix_opt = bsxfun(@plus,recon_matrix_opt,mean(orig_matrix)); % plug the mean back
if configs.non_negative == 1
    recon_matrix_opt(recon_matrix_opt<0) = 0;
end
recon_matrix_opt_test = recon_matrix_opt(:,Test_index);
recon_matrix_opt_retest = recon_matrix_opt(:,Retest_index);
if configs.binary == 1
    recon_matrix_opt_test(recon_matrix_opt_test>0.01|recon_matrix_opt_test<-0.01) = 1;
    recon_matrix_opt_test(recon_matrix_opt_test<=0.01&recon_matrix_opt_test>=-0.01) = 0;
    recon_matrix_opt_retest(recon_matrix_opt_retest>0.01|recon_matrix_opt_retest<-0.01) = 1;
    recon_matrix_opt_retest(recon_matrix_opt_test<=0.01&recon_matrix_opt_retest>=-0.01) = 0;
end
%% Compute Recon Identifiability matrix at optimal point
if strcmp(configs.score,'correlation')
    Ident_mat_recon_opt = corr(recon_matrix_opt_test,recon_matrix_opt_retest);
elseif strcmp(configs.score,'frobenius')
    Ident_mat_recon_opt = pdist2(recon_matrix_opt_test',recon_matrix_opt_retest','euclidean')./configs.numEdges;
elseif strcmp(configs.score,'geodesic')
    Ident_mat_recon_opt = f_GD_group(recon_matrix_opt_test,recon_matrix_opt_retest,configs.numRegions);
elseif strcmp(configs.score,'minkowski')
    Ident_mat_recon_opt = pdist2(recon_matrix_opt_test',recon_matrix_opt_retest',configs.score,configs.DistParam);
elseif strcmp(configs.score,'corrDist')
    Ident_mat_recon_opt = pdist2(recon_matrix_opt_test',recon_matrix_opt_retest','seuclidean');
elseif strcmp(configs.score,'hamming')
    Ident_mat_recon_opt = pdist2(recon_matrix_opt_test',recon_matrix_opt_retest','hamming');
elseif strcmp(configs.score,'jaccard')
    Ident_mat_recon_opt = pdist2(recon_matrix_opt_test',recon_matrix_opt_retest','jaccard');
elseif strcmp(configs.score,'spearman')
    Ident_mat_recon_opt = pdist2(recon_matrix_opt_test',recon_matrix_opt_retest','spearman');
end

