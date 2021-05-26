clear all
close all
clc
addpath('utils/')
models_name = {'center_prior_3d+4', 'flat_prior', 'flat_prior_3d+4', ...
                'judd_model_1_prior_3d+4', 'judd_model_2_prior_3d+4', 'judd_model_3_prior_3d+4', 'judd_model_4_prior_3d+4', 'judd_model_5_prior_3d+4', ...
                'mlnet_prior_3d+4_all', 'mlnet_prior_3d+4', 'noisy_prior_3d+4'};

%for i=1:134
i = 1;
    path = char(strcat('../matrix/images/info_per_subj_img_', num2str(i), '.mat')); fprintf('%s\n',path);
    load(path)
    
    thr = 4
    s = {};
    for su = 1:length(info_per_subj)
        if (info_per_subj(su).target_found==1 && size(info_per_subj(su).fixations_matrix_reduced,1)>thr)
            s = [s info_per_subj(su).fixations_matrix_reduced];
        end
    end
    
%%
  
U   = s{1};
AU  = sqrt(sum((U(1:end-1,:) - U(2:end,:))'.^2)');

V   = s{2};
AV  = sqrt(sum((V(1:end-1,:) - V(2:end,:))'.^2)');

M = abs(repmat(AU,[1 length(AU)]) - repmat(AV',[length(AU) 1]));

M = [1 2 3; 4 5 6];
A = zeros(prod(size(M)));
s = [];
t = [];
weights = [];
for i=1:size(M,1)
    for j=1:size(M,2)
        disp([i j])
        if (j < size(M,2))  
            A((i-1)*size(M,2)+j,(i-1)*size(M,2)+j+1)    = M(i,j+1);
            s       = [s (i-1)*size(M,2)+j];
            t       = [t (i-1)*size(M,2)+j+1];
            weights = [weights M(i,j+1)];
        end        
        if (i < size(M,1))  
            A((i-1)*size(M,2)+j,i*size(M,2)+j)          = M(i+1,j);
            s       = [s (i-1)*size(M,2)+j];
            t       = [t i*size(M,2)+j];
            weights = [weights M(i+1,j)];
        end        
        if (i < size(M,1) && j < size(M,2))  
            A((i-1)*size(M,2)+j,i*size(M,2)+j+1)        = M(i+1,j+1);
            s       = [s (i-1)*size(M,2)+j];
            t       = [t i*size(M,2)+j+1];
            weights = [weights M(i+1,j+1)];
        end        
    end
end
A
G = graph(s,t,weights);
plot(G,'EdgeLabel',G.Edges.Weight)

[path,d] = shortestpath(G,1,prod(size(M))); % Chequear si es dirigido.


%%
    cols = rainbow_colors(length(s));
%     ind = 2:3;
    figure;
        hold on
            for ii = 1:length(s)
%                 plot(s{ii}(ind,1)+(2*rand-1)*0.1,s{ii}(ind,2)+(2*rand-1)*0.1,'.-','Color',cols(ii,:))
                plot(s{ii}(end-1:end,1)+(2*rand-1)*0.1,s{ii}(end-1:end,2)+(2*rand-1)*0.1,'.-','Color',cols(ii,:))
                plot(s{ii}(end,1)+(2*rand-1)*0.1,s{ii}(end,2)+(2*rand-1)*0.1,'o','Color',cols(ii,:))
            end
        hold off
        set(gca,'XLim',[0.5 20.5],'YLim',[0.5 20.5])
        set(gca,'XTick',0:21,'YTick',0:21)
        grid on
