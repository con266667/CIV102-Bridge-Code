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
    if max(SFDi(i,:)) > max_sfd_val
        max_sfd_val = max(SFDi(i,:));
        i_max_sfd = i;
    end
    BMDi(i,:) = cumsum(SFDi(i,:));
    if max(BMDi(i,:)) > max_bmd_val
        max_bmd_val = max(BMDi(i,:));
        i_max_bmd = i;
    end
end

% Get max SFD, BMD and plot
SFD = max(abs(SFDi));
BMD = max(BMDi);

tiledlayout(2, 1);

nexttile;
sfdplot = plot(SFD, 'LineWidth', 3);
title("Shear Force Diagram");
nexttile;
bmdplot = plot(BMD, 'LineWidth', 3);
title("Bending Moment Diagram");

% plot(SFD, 'LineWidth', 3);
% plot(BMD, 'LineWidth', 3);


max_V = max(SFD);
max_M = max(BMD);

diaphram_max = 400;

design0param = [
    0, 100, 1, 80, 73.73, 1, 5, 1, 1;
    L, 100, 1, 80, 73.73, 1, 5, 1, 1;
];

param = design0param;

% param = [
%     0, 100, 3, 65, 90, 1, 5, 1, 1;
%     400, 80, 3, 65, 70, 1, 5, 1, 1;
%     800, 80, 3, 65, 70, 1, 5, 1, 1;
%     L, 100, 3, 65, 90, 1, 5, 1, 1;
% ];

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

total_volume = deck_area + walls_area + tab_area + bottom_area;
total_volume = sum(total_volume);

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
minFOS = [min(FOS_tens), min(FOS_comp), min(FOS_shear), min(FOS_glue), min(FOS_buck1), min(FOS_buck2), min(FOS_buck3), min(FOS_buckV)];
minFOS = min(minFOS);
fail_P = minFOS * 400;


%% 8. Vfail and Mfail
Mf_tens  = max_M * FOS_tens;
Mf_comp  = max_M * FOS_comp;
Vf_shear = max_V * FOS_shear;
Vf_glue  = max_V * FOS_glue;
Mf_buck1 = max_M * FOS_buck1;
Mf_buck2 = max_M * FOS_buck2;
Mf_buck3 = max_M * FOS_buck3;
Vf_buckV = max_V * FOS_buckV;

% Remove line on right, it looks nicer
max_abs_SFD = max(abs(SFDi));
max_abs_SFD(L + 1) = max_abs_SFD(L);

%% 9. Output plots of Vfail and Mfail

set(gca,'fontname','Inter')

subplot(2,3,1)
hold on; grid on; grid minor;
plot(x, Vf_shear, 'r', 'LineWidth', 2)
plot(x, max_abs_SFD, 'k', 'LineWidth', 2);
legend('Matboard Shear Failure', 'FontName', 'Inter') 
xlabel('Distance along bridge (mm)', 'FontName', 'Inter') 
ylabel('Shear Force (N)', 'FontName', 'Inter') 

subplot(2,3,2)
hold on; grid on; grid minor;
plot(x, Vf_glue, 'r', 'LineWidth', 2)
plot(x, max_abs_SFD, 'k', 'LineWidth', 2);
legend('Glue Shear Failure', 'FontName', 'Inter') 
xlabel('Distance along bridge (mm)', 'FontName', 'Inter') 
ylabel('Shear Force (N)', 'FontName', 'Inter') 

subplot(2,3,3)
hold on; grid on; grid minor;
plot(x, Vf_buckV, 'r', 'LineWidth', 2)
plot(x, max_abs_SFD, 'k', 'LineWidth', 2);
legend('Matboard Shear Buckling Failure', 'FontName', 'Inter') 
xlabel('Distance along bridge (mm)', 'FontName', 'Inter') 
ylabel('Shear Force (N)', 'FontName', 'Inter') 

subplot(2,3,4)
hold on; grid on; grid minor;
plot(x, Mf_tens, 'r', 'LineWidth', 2)
plot(x, Mf_comp, 'b', 'LineWidth', 2)
plot(x, max(abs(BMDi)), 'k', 'LineWidth', 2);
legend('Matboard Tension Failure', 'Matboard Compression Faliure', 'FontName', 'Inter') 
xlabel('Distance along bridge (mm)', 'FontName', 'Inter') 
ylabel('Bending Moment (Nmm)', 'FontName', 'Inter') 


subplot(2,3,5)
hold on; grid on; grid minor;
plot(x, Mf_buck1, 'r', 'LineWidth', 2)
plot(x, Mf_buck2, 'b', 'LineWidth', 2)
plot(x, max(abs(BMDi)), 'k', 'LineWidth', 2);
legend('Matboard Buckling Failure, Top Flange - Mid', 'Matboard Buckling Failure, Top Flange - Sides', 'FontName', 'Inter') 
xlabel('Distance along bridge (mm)', 'FontName', 'Inter') 
ylabel('Bending Moment (Nmm)', 'FontName', 'Inter') 


subplot(2,3,6)
hold on; grid on; grid minor;
plot(x, Mf_buck3, 'r', 'LineWidth', 2)
plot(x, max(abs(BMDi)), 'k', 'LineWidth', 2);
legend('Matboard Buckling Failure, Webs', 'FontName', 'Inter') 
xlabel('Distance along bridge (mm)', 'FontName', 'Inter') 
ylabel('Bending Moment (Nmm)', 'FontName', 'Inter') 


volume_text_col = 'g';
vol_allow = 1.27 * 812.8 * 1016;
if total_volume > vol_allow
    volume_text_col = 'r';
elseif total_volume > vol_allow * 0.8
    volume_text_col = 'y';
end

fail_text_col = 'k';
if fail_P < P
    fail_text_col = 'r';
elseif fail_P < P * 1.5
    fail_text_col = 'y';
elseif fail_P > 1000
    fail_text_col = 'g';
end


annotation('textbox', [0.01, 0.47, 0.1, 0], 'string', "Volume: " + round(total_volume) + " mm^3", 'FontSize', 18, 'Color', volume_text_col, 'FontWeight','bold', 'FontName', 'Inter');
annotation('textbox', [0.01, 0.57, 0.1, 0], 'string', "Failure Load: " + round(fail_P) + " N", 'FontSize', 18, 'Color', fail_text_col, 'FontWeight','bold', 'FontName', 'Inter');

