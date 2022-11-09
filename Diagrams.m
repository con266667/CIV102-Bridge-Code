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
    w(2) = A_y;
    w(locations(locations>0)) = w(locations(locations>0)) -(P/6);
    w(L) = B_y;

    % Create SFD and BMD vectors with integration
    SFDi(i,:) = cumtrapz(w');
    BMDi(i,:) = cumtrapz(SFDi(i,:));
end

% Get max SFD, BMD and plot
SFD = max(abs(SFDi));
BMD = max(BMDi);

figure;
plot(SFDi, "LineWidth", 3);
figure;
plot(BMDi, "LineWidth", 3);