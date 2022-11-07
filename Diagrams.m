clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis (0 to 1201)

%% 1. SFD, BMD under train loading

x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train=[1 1 1 1 1 1]*P/6; 
n_train = 3;

SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations

for i = 1:n_train
    location = L / (n_train+1) * i; % start location of train (300, 600, 900)
    
    % Reaction forces (Rearranged moment and Fy equations)
    B_y = location * P / L;
    A_y = P - B_y;

    % Create a vector with forces at locations
    w = zeros(n+1, 1);
    w(2) = A_y;
    w(location) = -P;
    w(L) = B_y;

    % Create SFD and BMD vectors with integration
    SFDi(i,:) = cumtrapz(w');
    BMDi(i,:) = cumtrapz(SFDi(i,:));
end

% Get max SFD, BMD and plot

SFD = max(abs(SFDi));
BMD = max(BMDi);

figure;
plot(SFD, "LineWidth", 3);
figure;
plot(BMD, "LineWidth", 3);