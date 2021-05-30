function info_per_subj = reduce_scanpaths(info_per_subj, delta, image_size)
    
    % INPUT: 
    %   - info_per_subj = struct that contains 
    %   - delta         = grid size (each grid's cell size is delta x delta)
    %   - image_size    = image size expected as (width, hight)
    % OUTPUT:
    %   - info_per_subj = struct containing all previous information plus
    %   time and space reduction
    % NOTE: The grid keeps the increase direction of fixation position
    % NOTE: (x,y)=(0,0) is in the top-left corner
    % NOTE: The fixations that fall outside the image were pulled inside
    % WARNING: This function does not save the info_per_subj struct anymore

    for tr = 1:length(info_per_subj)
        fixations = [info_per_subj(tr).x' info_per_subj(tr).y'];
        if not(isempty(fixations))
            timeFix = [info_per_subj(tr).t' info_per_subj(tr).t'+info_per_subj(tr).dur' ...
                info_per_subj(tr).dur'] ;
            fixations_matrix = nan(size(fixations));
            for i = 1:size(fixations,1)
                fixations_matrix(i,:) = reduce_fixation(fixations(i,:), delta, image_size);
            end

            mascara = [1 sum(diff(fixations_matrix)'~=0)];
            fixations_matrix_reduced    = fixations_matrix(1,:);
            timefix_matrix_reduced      = timeFix(1,:);
            for i = 2:size(fixations_matrix,1)
                if mascara(i)==0
                    timefix_matrix_reduced(end,:) = [timefix_matrix_reduced(end,1), ...
                                                        timeFix(i,2),...
                                                        timeFix(i,2)-timefix_matrix_reduced(end,1)];
                else
                    fixations_matrix_reduced    = [fixations_matrix_reduced;fixations_matrix(i,:)];
                    timefix_matrix_reduced      = [timefix_matrix_reduced;  timeFix(i,:)];
                end
            end

            % pull fixations outside the grid inside it
            % Warning, be
            indtmp = fixations_matrix_reduced(:,1)<1; 
            if ~isempty(indtmp); fixations_matrix_reduced(indtmp,1)=1; end
            
            indtmp = fixations_matrix_reduced(:,2)<1;
            if ~isempty(indtmp); fixations_matrix_reduced(indtmp,2)=1; end
            
            indtmp = fixations_matrix_reduced(:,1)>floor(image_size(2)/delta);
            if ~isempty(indtmp); fixations_matrix_reduced(indtmp,1)=floor(image_size(2)/delta); end
            
            indtmp = fixations_matrix_reduced(:,2)>floor(image_size(1)/delta);
            if ~isempty(indtmp); fixations_matrix_reduced(indtmp,2)=floor(image_size(1)/delta); end
            
            x_grid = fixations_matrix_reduced(:,2);
            y_grid = fixations_matrix_reduced(:,1);
            
            info_per_subj(tr).x_grid = x_grid;
            info_per_subj(tr).y_grid = y_grid;
            info_per_subj(tr).fixations_matrix_reduced = fixations_matrix_reduced;
            info_per_subj(tr).timefix_matrix_reduced = timefix_matrix_reduced;
            info_per_subj(tr).mascara_matrix_reduced = mascara;
            info_per_subj(tr).delta = delta;
        end
        %keyboard
        % save(path, 'info_per_subj');
    
        %     if doplot
        %         figure(1); clf
        %             ax1 = axes('position',[0.1 0.1 0.7 0.8]);
        %                 imagesc(max(prior(:))-prior)
        %                 xlim([0.5 size(prior,2)+0.5])
        %                 ylim([0.5 size(prior,1)+0.5])
        %         %         colormap bone
        %                 colormap(ax1,'bone')
        %             ax2 = axes('position',[0.1 0.1 0.7 0.8]);
        %                 hold on
        %                     plot(fixations_matrix_reduced(:,1),fixations_matrix_reduced(:,2),'k-')
        %                     scatter(fixations_matrix_reduced(:,1),fixations_matrix_reduced(:,2),...
        %                             timefix_matrix_reduced(:,3)/7,timefix_matrix_reduced(:,1)-timefix_matrix_reduced(1,1),'filled')
        %                             xlim([0.5 size(prior,2)+0.5])
        %                             ylim([0.5 size(prior,1)+0.5])
        %                 hold off
        %                 set(gca,'Color','none','visible','off');
        %         %         colormap jet
        %                 colormap(ax2,'jet');
        %                 hc = colorbar;
        %                 set(hc,'Position',[0.85 0.1 0.05 0.8]);
        %                 ylabel(hc,'Time Fix (ms)');
        %         % mappedFix = reduce_fixation(fixation, dims, top_left, delta)
        %     end
    end
end

