function [dist,path,M_assignment] = shortestPath(M)
% function [M_assignment] = shortestPath(M) %added by cstreet because dist variable not assigned when graphshortestpath function missing
global imwritePath fontSize debugMode

szM = size(M); %number of euclidean saccade distances passed to this variable
AM = zeros(szM(1)*szM(2)); %create vector the same length as the number of items in M
M_assignment = zeros(szM); %predefine variable

% Create adjacency matrix
for i = 1:szM(1) %for all the rows of M
    for j = 1:szM(2) %for all the columns of M
        
        currentNode = (i-1)*szM(2) + j;
        M_assignment(i,j) = currentNode;
        i_inc = 1;
        j_inc = 1;
        
        %if in the last row but not the last column
        if i == szM(1) && j < szM(2)
            i_inc = 0;
            node(i,j).adjacentNr = [currentNode, currentNode + 1]; % To the right
            node(i,j).weight = [M(i,j+j_inc)];
            AM(currentNode,currentNode + 1) = 1;

        %if in the last column but not the last row
        elseif i < szM(1) && j == szM(2)
             j_inc = 0;
             node(i,j).adjacentNr = [currentNode, currentNode + szM(2)];% Below 
             node(i,j).weight = [M(i+i_inc,j)];
             AM(currentNode, currentNode + szM(2)) = 1;
             
        %if in the last row and the last column
        elseif i == szM(1) && j == szM(2)
            i_inc = 0;
            j_inc = 0;
            node(i,j).adjacentNr  = [currentNode currentNode];
            node(i,j).weight = [0];
            
        else %otherwise
            node(i,j).adjacentNr = [currentNode, currentNode + 1;... % To the right
                                currentNode, currentNode + szM(2);...% Below  
                                currentNode, currentNode + szM(2) + 1];% Below  to the right
            node(i,j).weight = [M(i,j+j_inc),M(i+i_inc,j),M(i+i_inc, j+j_inc)];
            AM(currentNode,currentNode + 1) = 1;
            AM(currentNode, currentNode + szM(2)) = 1;
            AM(currentNode, currentNode + szM(2) + 1) = 1;
                            
        end
            
       
      
    end
end

G = cat(1,node.adjacentNr);
W = cat(2,node.weight)';

Graph = sparse(G(:,1),G(:,2),W);

if debugMode == 1
    h = view(biograph(Graph,[],'ShowWeights','on'));
end

% g = biograph.bggui(h);
% f = figure();
% copyobj(g.biograph.hgAxes,f);
% pause(2)
% delete([imwritePath,'Graph.pdf'])
% set(gca,'fontsize',fontSize)
% set(gcf,'color','none'),set(gca,'color','none')
% export_fig([imwritePath,'Graph.pdf'],'-pdf')
% 
% close(f);
% close(g.hgFigure);

%**BIOINFORMATIC TOOLBOX FUNCTION
[dist, path, pred] = graphshortestpath(Graph, 1, szM(1)*szM(2));

if debugMode == 1

    set(h.Nodes(path),'Color',[1 0.4 0.4])
    edges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
    set(edges,'LineColor',[1 0 0])
    set(edges,'LineWidth',1.5)
    pause(2)
    delete(h)
end

% g = biograph.bggui(h);
% f = figure();
% copyobj(g.biograph.hgAxes,f);
% pause(2)
% delete([imwritePath,'GraphShortestPath.pdf'])
% set(gca,'fontsize',fontSize)
% set(gcf,'color','none'),set(gca,'color','none')
% export_fig([imwritePath,'GraphShortestPath.pdf'],'-pdf')
% close(f);
% close(g.hgFigure);
% spy(AM)


% return
% AM = cat(2,node.adjacent)'
% W =  cat(2,node.weight)';
% 
% idx = unique(AM,'rows')