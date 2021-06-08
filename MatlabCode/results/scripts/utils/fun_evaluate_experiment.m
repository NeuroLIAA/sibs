function [correct,MaxFix,NFix,Pc] = fun_evaluate_experiment(Nfix_img_model, N)
    maxFix = [2 4 8 12];
    propTr = [13.4 14.9 29.9 41.8];
    nTrial = round(propTr * 134/100);
    filtro = [repmat(maxFix(1),nTrial(1),1);...
                repmat(maxFix(2),nTrial(2),1);...
                repmat(maxFix(3),nTrial(3),1);...
                repmat(maxFix(4),nTrial(4),1)];
            
    correct = nan(length(Nfix_img_model),N);
    MaxFix  = nan(length(Nfix_img_model),N);
    NFix    = nan(length(Nfix_img_model),N);
    Pc      = nan(length(maxFix),N);
    for i=1:N
        MaxFix(:,i)         = shuffle(filtro);
        correct(:,i)        = (Nfix_img_model <= MaxFix(:,i));
        NFix(find(correct(:,i)),i)= Nfix_img_model(find(correct(:,i)));
        
        for j = 1:length(maxFix)
            Pc(j,i) = nanmean(correct(MaxFix(:,i)==maxFix(j),i)); 
        end
    end
    
end
