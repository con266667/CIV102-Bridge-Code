clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis (0 to 1201)

%% 1. SFD, BMD under train loading

x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train=[1 1 1 1 1 1]*P/6; 
n_train = 2400;

SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations

for i = 1:n_train
    locations = x_train - 1200 + fix(i);
    locations(0 > locations | locations > 1201) = 0;
    loads = P_train;
    loads(locations == 0) = 0;
    
    
    % Reaction forces (Rearranged moment and Fy equations)
    B_y = sum(locations .* loads) / L;
    A_y = sum(loads) - B_y;

    % Create a vector with forces at locations
    w = zeros(n+1, 1);
    w(1) = A_y;
    w(locations(locations>0)) = w(locations(locations>0)) -(P/6);
    w(L) = B_y;

    % Create SFD and BMD vectors with integration
    SFDi(i,:) = cumsum(w');
    BMDi(i,:) = cumtrapz(SFDi(i,:));
end

% Get max SFD, BMD and plot
SFD = max(abs(SFDi));
BMD = max(BMDi);

% plot(SFD)

% max_V = max(SFD);
max_V = 240;
% max_M = max(BMD);
max_M = 6.75 * 10^4;

param = [
    0, 100, 1, 80, 73.73, 1, 5, 1, 1;
    400, 100, 1, 80, 73.73, 1, 5, 1, 1;
    800, 100, 1, 80, 73.73, 1, 5, 1, 1;
    L, 100, 1, 80, 73.73, 1, 5, 1, 1;
];

deck_width = interp1(param(:,1), param(:,2), x);
deck_layers = interp1(param(:,1), param(:,3), x);
bottom_width = interp1(param(:,1), param(:,4), x);
bottom_layers = interp1(param(:,1), param(:,9), x);
wall_height = interp1(param(:,1), param(:,5), x);
wall_layers = interp1(param(:,1), param(:,6), x);
tab_width = interp1(param(:,1), param(:,7), x);
tab_layers = interp1(param(:,1), param(:,8), x);


deck_area = deck_width .* deck_layers .* 1.27;
deck_dist = bottom_layers * 1.27 + wall_height + deck_layers * 1.27 * 0.5;
walls_area = wall_height .* wall_layers * 1.27 * 2;
walls_dist = bottom_layers * 1.27 + wall_height * 0.5;
tab_area = tab_width .* tab_layers * 1.27 * 2;
tab_dist = bottom_layers * 1.27 + wall_height - tab_layers * 1.27 * 0.5;
bottom_area = bottom_width .* bottom_layers * 1.27;
bottom_dist = bottom_layers * 1.27 * 0.5;

ybar = deck_area .* deck_dist + walls_area .* walls_dist + tab_area .* tab_dist + bottom_area .* bottom_dist;
ybar = ybar / (deck_area + walls_area + tab_area + bottom_area);

ytop = bottom_layers * 1.27 + wall_height + deck_layers * 1.27 - ybar;
ybot = ybar;

% I0
I = deck_width .* (deck_layers * 1.27) .^ 3;
I = I + 2 * wall_layers .* 1.27 .* (wall_height .^ 3);
I = I + 2 * tab_width .* (tab_layers * 1.27) .^ 3;
I = I + bottom_width .* (bottom_layers * 1.27) .^ 3;
I = I / 12;

% Ad^2
I = I + deck_area .* (deck_dist - ybar) .^ 2;
I = I + walls_area .* (walls_dist - ybar) .^ 2;
I = I + tab_area .* (tab_dist - ybar) .^ 2;
I = I + bottom_area .* (bottom_dist - ybar) .^ 2;

Qcent = bottom_area .* (ybot - bottom_layers * 1.27 * 0.5);
Qcent = Qcent + (ybot - (bottom_layers * 1.27)) .^ 2 .* wall_layers * 1.27;

Qglue = deck_area .* (ytop - deck_layers * 0.5 * 1.27);

S_top = max_M .* ytop ./ I;
S_bottom = max_M .* ybot ./ I;

T_cent = Qcent * max_V ./ (I * 2 .* wall_layers * 1.27);
T_glue = Qglue * max_V ./ (I * 2 .* (tab_width + wall_layers * 1.27));

E  = 4000;
mu = 0.2;
S_tens  = 30;
S_comp  = 6;
T_max   = 4;
T_gmax  = 2;

S_buck1 = (4 * pi^2 * E) / (12 * (1 - mu^2));
S_buck1 = ((deck_layers * 1.27) ./ (bottom_width - wall_layers * 1.27)) .^ 2 * S_buck1;

S_buck2 = (0.425 * pi^2 * E) / (12 * (1 - mu^2));
S_buck2 = ((deck_layers * 1.27) ./ ((deck_width - bottom_width) * 0.5 + (wall_layers * 1.27) * 0.5)) .^ 2 * S_buck2;


S_buck3 = (6 * pi^2 * E) / (12 * (1 - mu^2));
S_buck3 = ((wall_layers * 1.27) ./ (ytop - deck_layers * 1.27)) .^ 2 * S_buck3;

S_glue = max_M * (ytop - deck_layers * 1.27) ./ I;

a = S_buck3 ./ S_glue * P;
a(1);

a = S_buck2 ./ S_top * P;
a(1);

% T_buck  =


