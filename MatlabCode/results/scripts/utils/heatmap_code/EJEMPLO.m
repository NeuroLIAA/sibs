set(0,'DefaultFigureWindowStyle','docked')

%% EJEMPLO 1

%%% Genero fake data de fijaciones

fix     = [434 317; 766 528; 625 325; 986 170];
[Nf, d] = size(fix);
L       = 50;

gaze_x = []; gaze_y = [];
for i=1:Nf    
    gaze_x = [gaze_x round(fix(i,1) + L*randn(1,100))];
    gaze_y = [gaze_y round(fix(i,2) + L*randn(1,100))];
end
    
%%% genera heatmap

im_file  = '/Users/gsolovey/Dropbox/work/2016/2016 - proyectos/eyetech/experimentos/basico/images/coke.jpg';

settings.filterSize = 300; % tama?o del filtro gaussiano en px.
settings.filterSD   = 30;  % SD del filtro gaussiano en px.
settings.max_nhist  = 50;  % controla el colormap. ajustarlo para que no sature el heatmap.
settings.alpha      = 0.8; % entre 0 y 1. cuanto mayor, mas se ven las fijaciones.

fun_heatmap(im_file, gaze_x, gaze_y, settings)


%% EJEMPLO 2: igual que el 1 solo cambia la simulacion de las fijaciones

im_file     = '/Users/gsolovey/Dropbox/work/2016/2016 - proyectos/eyetech/experimentos/basico/images/coke.jpg';
im          = imread(im_file);
[Ly, Lx, ~] = size(im);

%%% genera fijaciones
T = 20;          % tiempo de exploraci?n en seg
f = 30;          % frecuencia en Hz
dt = 1/f;        % en seg
N = floor(T/dt); % numero de samples

w_saccades = 3;  % prob por unidad de tiempo de una sacada
p_saccades = w_saccades * dt;

gaze_x = -ones(1,N);
gaze_y = -ones(1,N);

gaze_x(1) = 200;
gaze_y(1) = 200;

for i=2:N
    
    while( gaze_x(i)<1 || gaze_x(i)>Lx || gaze_y(i)<1 || gaze_y(i)>Ly )
        
        if rand < p_saccades
            
            L         = 500; % largo en px
            gaze_x(i) = gaze_x(i-1) + L * cos(2*3.14*rand);
            gaze_y(i) = gaze_y(i-1) + L * sin(2*3.14*rand);
            
        else
            
            L         = 5; % largo en px
            gaze_x(i) = gaze_x(i-1) + L * cos(2*3.14*rand);
            gaze_y(i) = gaze_y(i-1) + L * sin(2*3.14*rand);
            
        end
        
    end
end
    
%%% genera heatmap

im_file  = '/Users/gsolovey/Dropbox/work/2016/2016 - proyectos/eyetech/experimentos/basico/images/coke.jpg';

settings.filterSize = 300; % tama?o del filtro gaussiano en px.
settings.filterSD   = 30;  % SD del filtro gaussiano en px.
settings.max_nhist  = 50;  % controla el colormap. ajustarlo para que no sature el heatmap.
settings.alpha      = 0.9; % entre 0 y 1. cuanto mayor, mas se ven las fijaciones.

fun_heatmap(im_file, gaze_x, gaze_y, settings)

