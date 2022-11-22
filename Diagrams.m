clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis (1 to 1201)

%% 1. SFD, BMD under train loading

x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train=[1 1 1 1 1 1]*P/6; 
n_train_start = 1;
n_train = 345;

% n_train_start = 164;
% n_train = 164;

SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations

i_max_sfd = 0;
max_sfd_val = 0;
i_max_bmd = 0;
max_bmd_val = 0;

for i = n_train_start:n_train
    locations = x_train - 52 + fix(i);
    locations(1 > locations | locations > 1201) = 0;
    loads = P_train;
    loads(locations == 0) = 0;
    
    % Reaction forces (Rearranged moment and Fy equations)
    B_y = sum((locations - 1) .* loads) / L;
    A_y = sum(loads) - B_y;

    % Create a vector with forces at locations
    w = zeros(n+1, 1);
    w(1) = A_y;
    w(L + 1) = B_y;
    w(locations(locations>0)) = w(locations(locations>0)) - (P/6);

    % Create SFD and BMD vectors with integration
    SFDi(i,:) = cumsum(w');
%     if max(SFDi(i,:)) > max_sfd_val
%         max_sfd_val = max(SFDi(i,:));
%         i_max_sfd = i;
%     end
    BMDi(i,:) = cumsum(SFDi(i,:));
    if max(BMDi(i,:)) > max_bmd_val
        max_bmd_val = max(BMDi(i,:));
        i_max_bmd = i;
    end
end

% Get max SFD, BMD and plot
SFD = max(abs(SFDi));
BMD = max(BMDi);

max(BMD)

tiledlayout(2, 1);

nexttile;
sfdplot = plot(SFD, 'LineWidth', 3);
title("Shear Force Diagram");
nexttile;
bmdplot = plot(BMD, 'LineWidth', 3);
title("Bending Moment Diagram");

% plot(SFD, 'LineWidth', 3);
% plot(BMD, 'LineWidth', 3);


% plot(SFD)

% max_V = max(SFD);
max_V = 240;
% max_M = max(BMD);
max_M = 6.75 * 10^4;

diaphram_max = 400;

param = [
    0, 100, 1, 80, 73.73, 1, 5, 1, 1;
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

T_buck = (5 * pi^2 * E) / (12 * (1 - mu^2));
t = 1.27 * wall_layers;
b = wall_height;
a = diaphram_max;
T_buck =((t./b).^2 + (t./a).^2) * T_buck;


%% 6. FOS
FOS_tens    = S_tens ./ S_bottom;
FOS_comp    = S_comp ./ S_top;
FOS_shear   = T_max ./ T_cent;
FOS_glue    = T_gmax ./ T_glue;
FOS_buck1   = S_buck1 ./ S_comp;
FOS_buck2   = S_buck2 ./ S_comp;
FOS_buck3   = S_buck3 ./ S_comp;
FOS_buckV   = T_buck ./ T_cent;

%% 7. Min FOS and the failure load Pfail
% minFOS = min(min(FOS_tens), min(FOS_comp), min(FOS_shear), min(FOS_glue), min(FOS_buck1), min(FOS_buck2), min(FOS_buck3), min(FOS_buckV));

%% 8. Vfail and Mfail
Mf_tens  = max_M * FOS_tens;
Mf_comp  = max_M * FOS_comp;
Vf_shear = max_V * FOS_shear;
Vf_glue  = max_V * FOS_glue;
Mf_buck1 = max_M * FOS_buck1;
Mf_buck2 = max_M * FOS_buck2;
Mf_buck3 = max_M * FOS_buck3;
Vf_buckV = max_V * FOS_buckV;

%% 9. Output plots of Vfail and Mfail
subplot(2,3,1)
hold on; grid on; grid minor;
plot(x, Vf_shear, 'r')
plot(x, max(abs(SFDi)), 'k');
legend('Matboard Shear Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Shear Force (N)')

subplot(2,3,2)
hold on; grid on; grid minor;
plot(x, Vf_glue, 'r')
plot(x, max(abs(SFDi)), 'k');
legend('Glue Shear Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Shear Force (N)')

subplot(2,3,3)
hold on; grid on; grid minor;
plot(x, Vf_buckV, 'r')
plot(x, max(abs(SFDi)), 'k');
legend('Matboard Shear Buckling Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Shear Force (N)')

subplot(2,3,4)
hold on; grid on; grid minor;
plot(x, Mf_tens, 'r')
plot(x, Mf_comp, 'b')
plot(x, max(abs(BMDi)), 'k');
legend('Matboard Tension Failure', 'Matboard Compression Faliure') 
xlabel('Distance along bridge (mm)') 
ylabel('Bending Moment (Nmm)')


subplot(2,3,5)
hold on; grid on; grid minor;
plot(x, Mf_buck1, 'r')
plot(x, Mf_buck2, 'b')
plot(x, max(abs(BMDi)), 'k');
legend('Matboard Buckling Failure, Top Flange - Mid', 'Matboard Buckling Failure, Top Flange - Sides') 
xlabel('Distance along bridge (mm)') 
ylabel('Bending Moment (Nmm)')


subplot(2,3,6)
hold on; grid on; grid minor;
plot(x, Mf_buck3, 'r')
plot(x, max(abs(BMDi)), 'k');
legend('Matboard Buckling Failure, Webs') 
xlabel('Distance along bridge (mm)') 
ylabel('Bending Moment (Nmm)')
