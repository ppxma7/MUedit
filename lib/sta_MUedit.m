function [stas, stas_norm, fig] = sta_MUedit(emg_grd, grid_order, bad_sig, firing, window_size)

% Input:
% emg_grd: emg signal of grid
% grid_order: electrode order of grid
% bad_sig: masked channels == 1
% firing: firing rate of a specific motor unit
% window_size: size of the time window that will use for calculating the spike triggered averaging (in frames)

% Output:
% sta: result of spike triggered averaging
% fig: figure handle
%
% Code Part:


emg_grd     = emg_grd(bad_sig==0,:);
grid_order  = grid_order(bad_sig==0,:);

electrodes = 1:size(emg_grd,1);

spans   = arrayfun(@(fir) fir-round(window_size/2)+1:1:fir+round(window_size/2),firing,'UniformOutput',false);
stas    = arrayfun(@(elc) stas_per_electrode(elc,spans,emg_grd), electrodes, 'UniformOutput', false); 

nroma = max(cellfun(@(sta) max([max(sta),abs(min(sta))]), stas));
stas_norm = cellfun(@(sta) sta.*(1/nroma), stas, 'UniformOutput', false);

fig = figure;
set(gcf, 'Color', 'w', 'WindowState', 'maximized');
t = tiledlayout(max(grid_order(:,1)), max(grid_order(:,2)), 'TileSpacing', 'tight', 'Padding', 'none'); % Create grid

% Use arrayfun to iterate over indices and plot each signal
arrayfun(@(elc) plotTile(elc, grid_order, stas_norm), electrodes);

% Define the function for plotting
function plotTile(elc, grid_order, stas_norm)
    ls = [0.05:1/max(grid_order(:,2)):1.05-1/max(grid_order(:,2))];
    us = sort([0:1/max(grid_order(:,1)):1-1/max(grid_order(:,1))],'descend');
    ax = subplot("Position", [ls(grid_order(elc,2)) us(grid_order(elc,1)) 1/max(grid_order(:,2)) 1/max(grid_order(:,1))*0.5]);
    
    
    pos     = get(ax, 'Position');
    newax   = [pos(1) pos(2) pos(3)*0.8 pos(4)*2];%us(grid_order(elc,1))
    plot(stas_norm{elc}, 'LineWidth', 2);
    ylim([-1 1]);
    axis off;
    set(ax, 'Position',newax);
end%function

function stas = stas_per_electrode(elc,spans,emg_signal)        
        sigFs = cellfun(@(spa) emg_signal(elc,spa), spans, 'UniformOutput', false);
        stas  = mean(cell2mat(sigFs'),1);
end%function

end%function