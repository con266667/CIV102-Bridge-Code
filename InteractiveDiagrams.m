n = 1200; % Discretize into 1 mm seg.
L = 1200; % Length of bridge
x = linspace(0, L, n+1); % x-axis (0 to 1201)
tiledlayout(2, 1);

nexttile;
sfdplot = plot(x,0*x);
nexttile;
bmdplot = plot(x,0*x);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
addlistener(h,'ContinuousValueChange', @(hObject, event) makeplot(hObject,sfdplot, bmdplot));

function makeplot(hObject,sfdplot,bmdplot)
    n = fix(get(hObject,'Value') * 2400);
    
    %% 0. Initialize Parameters
    L = 1200; % Length of bridge
    P = 400; % Total weight of train [N]
    
    %% 1. SFD, BMD under train loading
    
    x_train = [52 228 392 568 732 908]; % Train Load Locations
    P_train=[1 1 1 1 1 1]*P/6; 
    
    % Solve for SFD and BMD with the train at different locations
    
    locations = x_train - 1200 + fix(n);
    locations(0 > locations | locations > 1201) = 0;
    loads = P_train;
    loads(locations == 0) = 0;
    
    % Reaction forces (Rearranged moment and Fy equations)
    B_y = sum(locations .* loads) / L;
    A_y = sum(loads) - B_y;
    
    % Create a vector with forces at locations
    w = zeros(L+1, 1);
    w(2) = A_y;
    w(locations(locations>0)) = w(locations(locations>0)) - (P/6);
    w(L) = B_y;
    
    % Get SFD, BMD and plot
    SFD = cumtrapz(w');
    BMD = cumtrapz(SFD);
    set(sfdplot,'ydata',SFD);
    set(bmdplot,'ydata',BMD);
    drawnow;
end