% Computes the information (I) needed to travel through shortest-paths across pairs of nodes {source,target}.
% adj is the adjacency matrix of a weighted undirected connected graph
% Note that I(source,target) may be different that I(target,source) and
% hence matrix I is expected to be assymetric
function [I] = get_information_shortest_paths_wei_und(adj,SPL,B,K,flagComputeStepBack)

N=length(K);

I = zeros(N,N);

for i=1:N-1
    
    for j=i+1:N
        
        %% generalization, for weighted networks, of Rosvall et al. Networks and Cities: An Information Perspective

        path = retrieve_shortest_path(i,j,SPL,B);
        Kpath = K(path);

        w_step_ff = nan(1,length(path)-1);
        w_step_bk = nan(1,length(path)-1);
        for z=1:length(path)-1 %particle is ABLE to recognize where it came from and thus never go backwards
            w_step_ff(z) = adj(path(z),path(z+1));
            w_step_bk(z) = adj(path(z+1),path(z));
        end
        
        if flagComputeStepBack %particle does NOT recognize where it came from and thus can go backwards
            prob_sp_ff = (w_step_ff(1)./Kpath(1)).*prod(w_step_ff(2:end)./(Kpath(2:end-1)),2); %prob of shortest-path forward source->target
            prob_sp_bk = (w_step_bk(end)./Kpath(end)).*prod(w_step_bk(end-1:-1:1)./(Kpath(end-1:-1:2)),2); %prob of shortest-path backwards source<-target      
        else
            prob_sp_ff = (w_step_ff(1)./Kpath(1)).*prod(w_step_ff(2:end)./(Kpath(2:end-1)-w_step_ff(1:end-1)),2); %prob of shortest-path forward source->target
            prob_sp_bk = (w_step_bk(end)./Kpath(end)).*prod(w_step_bk(end-1:-1:1)./(Kpath(end-1:-1:2)-w_step_bk(end:-1:2)),2); %prob of shortest-path backwards source<-target    
        end   
        
        I(i,j) = -log(sum(prob_sp_ff));
        I(j,i) = -log(sum(prob_sp_bk));
    end
end
