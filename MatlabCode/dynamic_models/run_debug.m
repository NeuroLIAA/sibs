clear all

incfg.dinamic_model   = 'geisler';
incfg.iniimg 	      = 10;
incfg.endimg          = 10;
incfg.delta           = 32;
incfg.a               = 3;            % integers (?)   
incfg.b               = 4;            % integers (?)
incfg.static_model    = 'deepgaze';

incfg.debug = 1;
incfg.norm_cdf_tolerance = 0.0001;
% incfg.norm_cdf_tolerance = 0;
incfg.parfor = 0;
out = main_debug(incfg);

%% how to plot probabilistic maps

figure()
%imshow(out.p(:,:,1), [], 'InitialMagnification', 800, 'Colormap', parula)
imshow(out.detectability_map(:,:,1), [], 'InitialMagnification', 800, 'Colormap', parula)

%% Plot them togheter

f = figure(1);
    set(f, 'Color','w')
    f.Position = [0 0 1500 600];
    %f.OuterPosition = [0, 0, 1100, 420];
    Ax = gca;
    outerpos = Ax.OuterPosition;
    ti = Ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    Ax.Position = [left bottom ax_width ax_height];
    
    PosX = linspace(-0.13, 0, 4);
    PosY = [-0.04, -0.05];
    SizeX = 0.092;
    SizeY = 0.092;
    
    for i=1:4
        ax = subplot(2,4,i);
        fig = imshow(out.p(:,:,i), [], 'Colormap', parula);
        PosVec = ax.Position;
        OutPos = ax.OuterPosition;
        hold on
        plot(out.scanpath(i,2), out.scanpath(i,1),'r+', 'MarkerSize',15,'LineWidth', 3)
        ax.Position = PosVec+[PosX(i) PosY(1) SizeX SizeY];
        %ax.OuterPosition = OutPos+ [0, 0, 0, 0];
        %if i==1
        %    ylabel('Peron')
        %end

        ax = subplot(2,4,4+i);
        fig = imshow(out.detectability_map(:,:,i), [], 'Colormap', parula);
        [~, idx_y] = max(max(out.detectability_map(:,:,i)));
        [~, idx_x] = max(out.detectability_map(:,idx_y,i));
        PosVec = ax.Position;
        hold on
        plot(idx_y, idx_x,'r+', 'MarkerSize',15,'LineWidth', 3)
        ax.Position = PosVec+[PosX(i) PosY(2) SizeX SizeY];
    end

%%

