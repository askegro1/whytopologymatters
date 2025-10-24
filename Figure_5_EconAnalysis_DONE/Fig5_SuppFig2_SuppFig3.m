% =========================================================================
% File: Fig5_SuppFig2_SuppFig3.m
%
% Outputs:
%   - PlotResults_Fig5_SuppFig2_SuppFig3/Figure5_pdfformat.pdf
%   - PlotResults_Fig5_SuppFig2_SuppFig3/Figure5_figformat.fig
%   - PlotResults_Fig5_SuppFig2_SuppFig3/SuppFig2_pdfformat.pdf
%   - PlotResults_Fig5_SuppFig2_SuppFig3/SuppFig2_figformat.fig
%   - PlotResults_Fig5_SuppFig2_SuppFig3/SuppFig3_pdfformat.pdf
%   - PlotResults_Fig5_SuppFig2_SuppFig3/SuppFig3_figformat.fig
% =========================================================================

%% PREPARE
clearvars; close all; clc;

plot_dir  = "PlotResults_Fig5_SuppFig2_SuppFig3";
if ~isfolder(plot_dir), mkdir(plot_dir); end

%%%%%%%%%%%%%%%%%%%%%%% USER CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%
params                      = getParams();
paramsBaseline              = params;


% ========== PARAMETER RANGES ==========
paramRanges                 = getParamRanges();



%%%%%%%%%%%%%%%%%%%%%%% TEST IF THE FUNCTIONS WORK %%%%%%%%%%%%%%%%%%%%%%%
% ============== DIAGNOSTIC COMPUTATION ============== %
%resultBaseline             = computeNPC_DiscAndNonDisc(paramsBaseline);



%% %%%%%%%%%%%%%%%%%%%%%%% SAMPLE & PLOT MULTIVARIATE SENSITIVITY %%%%%%%%%%%%%%%%%%%%%%%
n_samples                                                               = 100000;


% ---------- Run the experiment ----------
[samples, delta_NPC, results_all, p_out_all]                            = sampleDeltaNPC_LHS(paramRanges, n_samples, paramsBaseline);

% Convert to kEUR
delta_NPC_kEuros                                                        = delta_NPC / 1000;


% ---------- Identify best, baseline, worst ----------
[maxNPC, idx_max]       = max(delta_NPC_kEuros);
resultBest              = results_all(idx_max);
pBest                   = p_out_all(idx_max);
idx.best                = idx_max;

idx_baseline            = n_samples;
resultBaseline          = results_all(idx_baseline);
baselineNPC             = resultBaseline.delta_NPC;
pBaseline               = p_out_all(idx_baseline);
idx.baseline            = idx_baseline;

[minNPC, idx_min]       = min(delta_NPC_kEuros);
resultWorst             = results_all(idx_min);
pWorst                  = p_out_all(idx_min);
idx.worst               = idx_min;



%% %%%%%%%%%%%%%%%%%%%%%%% MAIN FIGURE EXPORT %%%%%%%%%%%%%%%%%%%%%%%
% --- Margins ---    
ylims_Main                   = [-3 6];
yticks_Main                  = [-3 0 3 6];
dotSize                     = 1;
dotSize2                    = 10;

figMain = figure;

    tl = tiledlayout(4,3,'TileSpacing','compact','Padding','none');

    axAnnualCost = nexttile(tl, 1, [1 3]);  
        yLim_axAnnual = [-3 9];
        yTicks_axAnnual = [-3 0 3 6 9];          
        plotAnnualCostBreakdown_DiscNonDisc_Subplot(axAnnualCost, resultBaseline, yLim_axAnnual, yTicks_axAnnual, 'nondiscounted');
        addPanelLabel(axAnnualCost, '\textbf{(a)}', -0.065, 1, 12);     
    axBestCase = nexttile(tl);
        yLim_axBest = [10 30];
        yTicks_axBest = [10 15 20 25 30];
        plotCumulativeNPC_DiscNonDisc_Subplot(axBestCase, resultBest, true, true, "northwest", "Best", resultBest.delta_NPC, yLim_axBest, yTicks_axBest, 'discounted');
        addPanelLabel(axBestCase, '\textbf{(b)}', -0.25, 1, 12); 
    axBaselineCase = nexttile(tl);
        yLim_axBaseline = [5 25];
        yTicks_axBaseline = [5 10 15 20 25]; 
        plotCumulativeNPC_DiscNonDisc_Subplot(axBaselineCase, resultBaseline, true, true, "northwest", "Baseline", resultBaseline.delta_NPC, yLim_axBaseline, yTicks_axBaseline, 'discounted');
        addPanelLabel(axBaselineCase, '\textbf{(c)}', -0.235, 1, 12); 
    axWorstCase = nexttile(tl);
        yLim_axWorst = [0 60];
        yTicks_axWorst = [0 20 40 60];      
        plotCumulativeNPC_DiscNonDisc_Subplot(axWorstCase, resultWorst, true, true, "northwest", "Worst", resultWorst.delta_NPC, yLim_axWorst, yTicks_axWorst, 'discounted');
        addPanelLabel(axWorstCase, '\textbf{(d)}', -0.235, 1, 12); 
    axSpearman = nexttile(tl);
        tornadoPlotDeltaNPC_Subplot(axSpearman, samples, delta_NPC_kEuros);  
        addPanelLabel(axSpearman, '\textbf{(e)}', -0.235, 1, 12);          
    axEpacknom = nexttile(tl);
        x_axSupp2               = samples.E_pack_nom; 
        y_axSupp2               = delta_NPC_kEuros; 
        xLabel_axSupp2          = '$E^{\mathrm{nom}}_{\mathrm{pack}}$ [kWh]';
        xLim_axSupp2            = [10 130];
        xTicks_axSupp2          = 20:20:120;
        [x_zero_E, R_squared_E, RMSE_E] = plotSensitivity1D_NoPerc_Subplot(axEpacknom, x_axSupp2, y_axSupp2, xLabel_axSupp2, xLim_axSupp2, xTicks_axSupp2, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize2);
        addPanelLabel(axEpacknom, '\textbf{(f)}', -0.235, 1, 12);          
    axL = nexttile(tl);
        xlims_L = [5000 70000];
        xticks_L = 5000:20000:70000;      
        [x_zero_L, R_squared_L, RMSE_L] = plotSensitivity1D_NoPerc_Subplot(axL, samples.L, delta_NPC_kEuros, ...
        '$L$ [km]', xlims_L, xticks_L, ...
        ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize, "northwest");    
        addPanelLabel(axL, '\textbf{(g)}', -0.235, 1, 12); 
    axnu = nexttile(tl);
        x_axnu                  = samples.nu; 
        y_axnu                  = delta_NPC_kEuros; 
        xLabel_axnu             = '$\nu\,[\%]$';    
        xlims_axnu              = [0 16];
        xticks_axnu             = [0 5 10 15];
        [x_zero_nu, R_squared_nu, RMSE_nu] = plot1DWithRegression_Subplot(axnu, x_axnu, y_axnu, xLabel_axnu, xlims_axnu, xticks_axnu, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize); 
        addPanelLabel(axnu, '\textbf{(h)}', -0.235, 1, 12);                 
    axRectPlot = nexttile(tl);
        xLim_axRectPlot     = [0 15.1];
        xTicks_axRectPlot   = [0 5 10 15];
        yLim_axRectPlot     = [4500 82000];
        yTicks_axRectPlot   = [10000 30000 50000 70000];    
    axChemPlot = nexttile(tl);
        rect_bounds         = [paramRanges.nu.LIMITS; paramRanges.L.LIMITS];
        successThreshold    = 0.997;
        dotSize             = 1;
        [nu_best, L_best] = plotChemistryBoxScatterWithRectangle_Modified_Subplot_v2(axChemPlot, axRectPlot, samples, delta_NPC, ...
            rect_bounds, ylims_Main, yticks_Main, true, successThreshold, dotSize, ...
            xLim_axRectPlot, xTicks_axRectPlot, yLim_axRectPlot, yTicks_axRectPlot);
        addPanelLabel(axRectPlot, '\textbf{(i)}', -0.235, 1, 12); 
        addPanelLabel(axChemPlot, '\textbf{(j)}', -0.235, 1, 12);  
     
 
         

    screens = get(0,'MonitorPositions');
    targetMonitor = 2; 

    % A4 dimensions (80% height)
    paperWidth = 21.0;    % cm
    paperHeight = 29.7 * 0.9; % 23.76 cm

    % Get screen DPI
    dpi = get(0,'ScreenPixelsPerInch');
    cm_to_inch = 1/2.54;

    % Convert to pixels
    width_px = round(paperWidth * cm_to_inch * dpi);
    height_px = round(paperHeight * cm_to_inch * dpi);

    % Get monitor position
    mon = screens(targetMonitor,:);
    fig_left = mon(1);
    fig_bottom = mon(2);

    % Position figure in top-left of monitor with exact size
    set(figMain, 'Units', 'pixels');
    set(figMain, 'Position', [fig_left, fig_bottom + mon(4) - height_px - 100, width_px, height_px]);

    % Set PaperSize for export
    set(figMain, 'PaperUnits', 'centimeters');
    set(figMain, 'PaperSize', [paperWidth paperHeight]);
    set(figMain, 'PaperPosition', [0 0 paperWidth paperHeight]);
    set(figMain, 'PaperPositionMode', 'manual');

    %exportgraphics(figMain, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'Figure5_epsformat.eps'), 'ContentType','vector');
    exportgraphics(figMain, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'Figure5_pdfformat.pdf'), 'ContentType','image','Resolution',300);
    savefig(figMain, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'Figure5_figformat.fig'));

    


%% ============================ SUPPLEMENTARY FIGURE S2 EXPORT (A4 layout, EPS + FIG) ============================    
figSuppFigure2 = figure;

    tlSupp = tiledlayout(2,3,'TileSpacing','compact','Padding','none');
    axdeltaLoss = nexttile(tlSupp);
        xlims_deltaLoss = [1 5];
        xticks_deltaLoss = 1:1:5;
        [x_zero_delta, R_squared_delta, RMSE_delta] = plot1DWithRegression_Subplot(axdeltaLoss, samples.delta_loss, delta_NPC_kEuros, ...
            '$\delta_{\mathrm{loss}} \,[\%]$', xlims_deltaLoss, xticks_deltaLoss, ...
            ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize); 
        addPanelLabel(axdeltaLoss, '\textbf{(a)}', -0.2, 1, 12);
    axdeltaalpha = nexttile(tlSupp);
        samples.delta_alpha     = samples.alpha_CBP - samples.alpha_RBP;
        x_axSupp1               = samples.delta_alpha; 
        y_axSupp1               = delta_NPC_kEuros; 
        xLabel_axSupp1          = '$\Delta \alpha \,[\%]$';
        xLim_axSupp1            = [0.5 1.5];
        xTicks_axSupp1          = [0.5 1 1.5];
        plot1DWithRegression_Subplot(axdeltaalpha, x_axSupp1, y_axSupp1, xLabel_axSupp1, xLim_axSupp1, xTicks_axSupp1, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize);
        addPanelLabel(axdeltaalpha, '\textbf{(b)}', -0.2, 1, 12);   
    axr = nexttile(tlSupp);
        x_axr               = samples.r_CBP; 
        y_axr               = delta_NPC_kEuros; 
        xLabel_axr          = '$r_{\mathrm{CBP}} = r_{\mathrm{RBP}} \,[\%]$';    
        xlims_axr           = [2 4];
        xticks_axr          = [1 2 3 4 5];
        plotSensitivity1D_Subplot(axr, x_axr, y_axr, ...
        xLabel_axr, xlims_axr, xticks_axr, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize);    
        addPanelLabel(axr, '\textbf{(c)}', -0.2, 1, 12);    
    axYEV = nexttile(tlSupp);
        x_axYEV               = samples.Y_EV; 
        y_axYEV               = delta_NPC_kEuros; 
        xLabel_axYEV          = '$Y_{\mathrm{EV}}$ [years]';    
        xlims_axYEV           = [15 20];
        xticks_axYEV          = 15:1:20;
        plotSensitivity1D_Subplot(axYEV, x_axYEV/100, y_axYEV, xLabel_axYEV, xlims_axYEV, xticks_axYEV, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize);  
        addPanelLabel(axYEV, '\textbf{(d)}', -0.2, 1, 12); 
    axSOHCBP = nexttile(tlSupp);
        x_axSOHCBP               = samples.SOH_EOL_SecondLife_CBP; 
        y_axSOHCBP               = delta_NPC_kEuros; 
        xLabel_axSOHCBP          = '$SOH^{\mathrm{CBP}}_{EOL2} \,[\%]$';
        xLim_axSOHCBP            = [50, 60];
        xTicks_axSOHCBP          = [50 55 60];
        plotSensitivity1D_Subplot(axSOHCBP, x_axSOHCBP, y_axSOHCBP, xLabel_axSOHCBP, xLim_axSOHCBP, xTicks_axSOHCBP, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize);
        addPanelLabel(axSOHCBP, '\textbf{(e)}', -0.2, 1, 12);   
    axSOHRBP = nexttile(tlSupp);
        x_axSOHRBP               = samples.SOH_EOL_SecondLife_RBP; 
        y_axSOHRBP               = delta_NPC_kEuros; 
        xLabel_axSOHRBP          = '$SOH^{\mathrm{RBP}}_{EOL2} \,[\%]$';    
        xlims_axSOHRBP           = [40, 50];
        xticks_axSOHRBP          = [40 45 50];
        plotSensitivity1D_Subplot(axSOHRBP, x_axSOHRBP, y_axSOHRBP, ...
        xLabel_axSOHRBP, xlims_axSOHRBP, xticks_axSOHRBP, ylims_Main, yticks_Main, idx_max, idx_baseline, idx_min, dotSize);    
        addPanelLabel(axSOHRBP, '\textbf{(f)}', -0.2, 1, 12);   
       



    screens         = get(0,'MonitorPositions');
    targetMonitor   = 2; 

    % A4 dimensions (80% height)
    paperWidth      = 21.0;    % cm
    paperHeight     = 29.7 * 0.5; % 23.76 cm

    % Get screen DPI
    dpi             = get(0,'ScreenPixelsPerInch');
    cm_to_inch      = 1/2.54;

    % Convert to pixels
    width_px        = round(paperWidth * cm_to_inch * dpi);
    height_px       = round(paperHeight * cm_to_inch * dpi);

    % Get monitor position
    mon             = screens(targetMonitor,:);
    fig_left        = mon(1);
    fig_bottom      = mon(2);

    % Position figure in top-left of monitor with exact size
    set(figSuppFigure2, 'Units', 'pixels');
    set(figSuppFigure2, 'Position', [fig_left, fig_bottom + mon(4) - height_px - 100, width_px, height_px]);

    % Set PaperSize for export
    set(figSuppFigure2, 'PaperUnits', 'centimeters');
    set(figSuppFigure2, 'PaperSize', [paperWidth paperHeight]);
    set(figSuppFigure2, 'PaperPosition', [0 0 paperWidth paperHeight]);
    set(figSuppFigure2, 'PaperPositionMode', 'manual');

    %exportgraphics(figSuppFigure2, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure2_epsformat.eps'), 'ContentType','vector');
    exportgraphics(figSuppFigure2, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure2_pdfformat.pdf'), 'ContentType','image','Resolution',300);
    savefig(figSuppFigure2, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure2_figformat.fig'));



%% ============================ SUPPLEMENTARY FIGURE S3 EXPORT (A4 layout, EPS + FIG) ============================
%% ============================ SUPPLEMENTARY FIGURE S3 EXPORT (direct draw, no temp figs) ============================
% Spec per panel
xlims = {[0 15.1], [0 15.1], [1 5]};
ylims = {[4500 82000], [1 5.85714], [4500 82000]};
xticks = {[0 5 10 15], [0 5 10 15], [1 3 5]};
yticks = {[10000 30000 50000 70000], [1 3 5], [10000 30000 50000 70000]};
xVars = {"nu", "nu", "delta_loss"};
yVars = {"L", "delta_loss", "L"};
xLabels = {"$\nu\,[\%]$", "$\nu\,[\%]$", "$\delta_{\mathrm{loss}}\,[\%]$"};
yLabels = {"$L$ [km]", "$\delta_{\mathrm{loss}}\,[\%]$", "$L$ [km]"};
xInPerc = [true, true, true];
yInPerc = [false, true, false];
boundsList = {
    [paramRanges.nu.LIMITS;          paramRanges.L.LIMITS], ...
    [paramRanges.nu.LIMITS;          paramRanges.delta_loss.LIMITS], ...
    [paramRanges.delta_loss.LIMITS;  paramRanges.L.LIMITS]
};

dotSize = 1;  % as before
successThreshold = 0.997;

% Create final figure + layout
figSuppFigure3 = figure;
tlSupp_2D = tiledlayout(figSuppFigure3, 1, 3, 'TileSpacing','compact','Padding','none');

% Loop and draw directly into the tiles
for i = 1:3
    % Visible rectangle panel axis (goes into the tile)
    axRect = nexttile(tlSupp_2D);

    % Dummy/invisible axis for the chemistry panel (not used/seen)
    axChem = axes('Parent', figSuppFigure3, 'Visible', 'off', 'Units', 'normalized', 'Position', [0 0 0.01 0.01]);

    % Draw (function expects axChem first, axRect second)
    [~, ~] = plotDominantRectangleIn2DFeatureSpace( ...
        axChem, axRect, ...
        samples, delta_NPC, ...
        boundsList{i}, ylims_Main, yticks_Main, true, ...
        successThreshold, dotSize, ...
        xlims{i}, xticks{i}, ...
        ylims{i}, yticks{i}, ...
        xVars{i}, yVars{i}, xInPerc(i), yInPerc(i), ...
        xLabels{i}, yLabels{i}, ...
        0.95, true);

    % Optional: set TickLabelInterpreter to LaTeX consistently
    set(axRect, 'TickLabelInterpreter','latex');
end

% Positioning & export (unchanged)
screens = get(0,'MonitorPositions');
dpi = get(0,'ScreenPixelsPerInch');
cm_to_inch = 1/2.54;
paperWidth = 21.0;
paperHeight = 29.7 * 0.4;
width_px = round(paperWidth * cm_to_inch * dpi);
height_px = round(paperHeight * cm_to_inch * dpi);
mon = screens(2,:);
fig_left = mon(1);
fig_bottom = mon(2);

set(figSuppFigure3, 'Units', 'pixels', ...
    'Position', [fig_left, fig_bottom + mon(4) - height_px - 100, width_px, height_px]);
set(figSuppFigure3, 'PaperUnits', 'centimeters', ...
    'PaperSize', [paperWidth paperHeight], ...
    'PaperPosition', [0 0 paperWidth paperHeight], ...
    'PaperPositionMode', 'manual');

% Export (consider ContentType='image' for speed as noted earlier)
exportgraphics(figSuppFigure3, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure3_pdfformat.pdf'), 'ContentType','image','Resolution',300);
savefig(figSuppFigure3, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure3_figformat.fig'));
% 
% % Set up input variables for each subplot
% figHandles = {@()figure, @()figure, @()figure};
% axesHandles = cell(3,2);
% xlims = {[0 15.1], [0 15.1], [1 5]};
% ylims = {[4500 82000], [1 5.85714], [4500 82000]};
% xticks = {[0 5 10 15], [0 5 10 15], [1 3 5]};
% yticks = {[10000 30000 50000 70000], [1 3 5], [10000 30000 50000 70000]};
% xVars = {"nu", "nu", "delta_loss"};
% yVars = {"L", "delta_loss", "L"};
% xLabels = {"$\nu\,[\%]$", "$\nu\,[\%]$", "$\delta_{\mathrm{loss}}\,[\%]$"};
% yLabels = {"$L$ [km]", "$\delta_{\mathrm{loss}}\,[\%]$", "$L$ [km]"};
% xInPerc = [true, true, true];
% yInPerc = [false, true, false];
% boundsList = {
%     [paramRanges.nu.LIMITS; paramRanges.L.LIMITS],...
%     [paramRanges.nu.LIMITS; paramRanges.delta_loss.LIMITS],...
%     [paramRanges.delta_loss.LIMITS; paramRanges.L.LIMITS]
% };
% 
% % Main plotting loop
% for i = 1:3
%     fig1 = figHandles{i}();
%     axesHandles{i,1} = gca(fig1);
%     fig2 = figHandles{i}();
%     axesHandles{i,2} = gca(fig2);
% 
%     [~, ~] = plotDominantRectangleIn2DFeatureSpace( ...
%         axesHandles{i,1}, axesHandles{i,2}, ...
%         samples, delta_NPC, ...
%         boundsList{i}, ylims_Main, yticks_Main, true, ...
%         0.997, 1, ...
%         xlims{i}, xticks{i}, ...
%         ylims{i}, yticks{i}, ...
%         xVars{i}, yVars{i}, xInPerc(i), yInPerc(i), ...
%         xLabels{i}, yLabels{i}, ...
%         0.95, true);
% end
% 
% % Create final tiled layout
% figSuppFigure3 = figure;
% tlSupp_2D = tiledlayout(1,3,'TileSpacing','compact','Padding','none');
% 
% for i = 1:3
%     nexttile(tlSupp_2D);
%     axSrc = axesHandles{i,2};
%     axDest = gca;
% 
%     % Copy plot objects
%     newObjs = copyobj(allchild(axSrc), axDest);
%     % Set axis properties
%     title(axSrc.Title.String, 'Interpreter', 'latex');
%     xlabel(axSrc.XLabel.String, 'Interpreter', 'latex');
%     ylabel(axSrc.YLabel.String, 'Interpreter', 'latex');
%     xlim(axSrc.XLim);
%     ylim(axSrc.YLim);
%     set(gca, 'XTick', axSrc.XTick, 'YTick', axSrc.YTick);
%     set(gca, 'TickLabelInterpreter', 'latex');
% 
%     % Recreate legend (use findobj if specific handles not saved)
%     oldLeg = legend(axSrc);
%     if isvalid(oldLeg)
%         % Optional: use findobj to get lines if needed
%         plotLines = findobj(axDest, 'Type', 'line');
%         legend(flip(plotLines), oldLeg.String, ...
%             'Interpreter','latex', ...
%             'Location', oldLeg.Location, ...
%             'NumColumns', oldLeg.NumColumns, ...
%             'FontSize', oldLeg.FontSize);
%     end
% end
% 
% % % === Align top Y-tick ===
% % Step 1: Collect top Y-tick value for each source axis
% topTicks = cellfun(@(ax) max(ax.YTick), axesHandles(:,2));
% 
% % Step 2: Define a common top tick value (e.g., max of all top ticks)
% commonTopTick = max(topTicks);
% 
% % Step 3: For each axis, add its native ticks but ensure commonTopTick is included
% tileAxes = findall(figSuppFigure3, 'Type', 'Axes');
% 
% for k = 1:numel(tileAxes)
%     ax = tileAxes(k);
%     ticks = ax.YTick;
% 
%     % Add the commonTopTick if it's not already in the list and within YLim
%     if commonTopTick <= ax.YLim(2) && ~ismember(commonTopTick, ticks)
%         ticks(end+1) = commonTopTick; %#ok<AGROW>
%     end
% 
%     % Sort ticks and apply
%     ticks = unique(ticks);
%     ax.YTick = ticks;
% end
% 
% % Positioning (A4 top-left on monitor 2)
% screens = get(0,'MonitorPositions');
% dpi = get(0,'ScreenPixelsPerInch');
% cm_to_inch = 1/2.54;
% paperWidth = 21.0;
% paperHeight = 29.7 * 0.4;
% width_px = round(paperWidth * cm_to_inch * dpi);
% height_px = round(paperHeight * cm_to_inch * dpi);
% mon = screens(2,:);
% fig_left = mon(1);
% fig_bottom = mon(2);
% set(figSuppFigure3, 'Units', 'pixels');
% set(figSuppFigure3, 'Position', [fig_left, fig_bottom + mon(4) - height_px - 100, width_px, height_px]);
% set(figSuppFigure3, 'PaperUnits', 'centimeters', ...
%     'PaperSize', [paperWidth paperHeight], ...
%     'PaperPosition', [0 0 paperWidth paperHeight], ...
%     'PaperPositionMode', 'manual');
% 
% % Save figure
% %exportgraphics(figSuppFigure3, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure3_epsformat.eps'), 'ContentType','vector');
% exportgraphics(figSuppFigure3, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure3_pdfformat.pdf'), 'ContentType','image','Resolution',300);
% savefig(figSuppFigure3, fullfile("PlotResults_Fig5_SuppFig2_SuppFig3",'SuppFigure3_figformat.fig'));
% 


figs = findall(0, 'Type', 'figure');
close(setdiff(figs, [figMain, figSuppFigure2, figSuppFigure3]));


%% AUXILIARY FUNCTIONS
% -------------------------------------------------------------------------
% GETPARAMS
% -------------------------------------------------------------------------
function paramsOut = getParams(varargin)
%GETPARAMS Return economic parameters for battery pack modeling (with overrides).
%   PARAMS = GETPARAMS() returns a struct of default economic parameters for
%   battery systems (pack size, voltages, costs, lifetimes, etc.).
%
%   PARAMS = GETPARAMS(Name,Value,...) applies user-specified overrides for
%   any subset of supported fields (see INPUTS below). Unknown names error.
%
%   The function also sets a chemistry-specific lifetime multiplier for the
%   RBP (reconfigurable battery pack) via a percentage "chi" term.
%
%   INPUTS (Name-Value pairs)
%   All fields below can be overridden unless noted otherwise:
%     'E_pack_nom'             (1,1) double   Nominal pack energy [kWh]. Default 80.
%     'V_pack_nom'             (1,1) double   Nominal pack voltage [V].  Default 800.  
%     'V_nom_module'           (1,1) double   Nominal module voltage [V]. Default 50.
%     'Chemistry'              string|char    "LFP" (default) or "NMC".
%     'UserSelection'          (1,1) double   Placeholder for UI logic. Default 1.
%     'Y_EV'                   (1,1) double   Vehicle lifetime [years]. Default 18.8.
%     'Y_CBP'                  (1,1) double   Lifetime of conventional battery pack (CBP) [years]. Default 10.
%     'r_CBP'                  (1,1) double   Discount rate CBP (fraction). Default 0.03.
%     'r_RBP'                  (1,1) double   Discount rate RBP (fraction). Default 0.03.
%     'alpha_CBP'              (1,1) double   Annual O&M cost rate CBP. Default 0.02.
%     'alpha_RBP'              (1,1) double   Annual O&M cost rate RBP. Default 0.01.
%     'delta_loss'             (1,1) double   Relative increase in ohmic losses (RBP). Default 0.03.
%     'SOH_EOL_FirstLife_CBP'  (1,1) double   EOL SOH threshold CBP first life. Default 0.8.
%     'SOH_EOL_FirstLife_RBP'  (1,1) double   EOL SOH threshold RBP first life. Default 0.8.
%     'SOH_EOL_SecondLife_CBP' (1,1) double   EOL SOH threshold CBP second life. Default 0.55.
%     'SOH_EOL_SecondLife_RBP' (1,1) double   EOL SOH threshold RBP second life. Default 0.45.
%     'c_pack_USDperkWh'       (1,1) double   Initial pack cost [USD/kWh]. Default 115.
%     'currExch_USDtoEUR'      (1,1) double   USD→EUR exchange rate. Default 0.8554.
%     'c_en_EURperkWh'         (1,1) double   Energy price [EUR/kWh]. Default 0.24.
%     'eta_en_kWhkm'           (1,1) double   Vehicle energy use [kWh/km]. Default 0.2.
%     'L'                      (1,1) double   Annual distance [km]. Default 12000.
%     'nu'                     (1,1) double   Upfront cost difference (RBP vs CBP). Default 0.075.
%     'c_res_CBP_USDperkWh'    (1,1) double   Residual value CBP [USD/kWh]. Default 22.99.
%     'c_res_RBP_USDperkWh'    (1,1) double   Residual value RBP [USD/kWh]. Default 22.99.
%     'c_res_depreccoeff'      (1,1) double   Residual depreciation coefficient. Default 0.5.
%     'DefaultIndicator'       (1,1) logical  If TRUE use mean chi; if FALSE sample. Default true.
%
%
%   OUTPUTS
%     paramsOut  struct  All parameters after overrides, plus:
%       .c_res_CBP_EURperkWh  Converted residual value [EUR/kWh].
%       .c_res_RBP_EURperkWh  Converted residual value [EUR/kWh].
%       .Y_RBP                Chemistry/voltage-adjusted RBP lifetime [years].
%
%   NOTES
%   - If V_pack_nom is not 800, chemistry bounds are undefined and an error is thrown.
%   - Chi % mean (mu) and bounds follow the original log-linear forms by chemistry.
%

    % ====== DEFAULT PARAMETERS ======
    % Core defaults for economic inputs
    defaultParams = struct( ...
        'E_pack_nom',               80, ...            % [kWh]
        'V_pack_nom',               800, ...           % [V] 
        'V_nom_module',             50, ...            % [V]
        'Chemistry',                "LFP", ...         % "LFP" or "NMC"
        'UserSelection',            1, ...
        'Y_EV',                     18.8, ...          % [years]
        'Y_CBP',                    10, ...            % [years]
        'r_CBP',                    0.03, ...          % [fraction]
        'r_RBP',                    0.03, ...          % [fraction]
        'alpha_CBP',                0.02, ...          % [fraction/year]
        'alpha_RBP',                0.02 * 0.5, ...    % [fraction/year]
        'delta_loss',               0.03, ...          % [fraction]
        'SOH_EOL_FirstLife_CBP',    0.8, ...
        'SOH_EOL_FirstLife_RBP',    0.8, ...
        'SOH_EOL_SecondLife_CBP',   0.55, ...
        'SOH_EOL_SecondLife_RBP',   0.45, ...
        'c_pack_USDperkWh',         115, ...           % [USD/kWh]
        'currExch_USDtoEUR',        0.8554, ...
        'c_en_EURperkWh',           0.24, ...
        'eta_en_kWhkm',             0.2, ...
        'L',                        12000, ...         % [km/year]
        'nu',                       0.075, ...
        'c_res_CBP_USDperkWh',      22.99, ...         % [USD/kWh]
        'c_res_RBP_USDperkWh',      22.99, ...         % [USD/kWh]
        'c_res_depreccoeff',        0.5, ...
        'DefaultIndicator',         true...        
    );

    % Derived (currency conversion for residual values)
    defaultParams.c_res_CBP_EURperkWh = defaultParams.c_res_CBP_USDperkWh * defaultParams.currExch_USDtoEUR;
    defaultParams.c_res_RBP_EURperkWh = defaultParams.c_res_RBP_USDperkWh * defaultParams.currExch_USDtoEUR;


    % ====== HANDLE USER INPUTS ======
    % If no arguments are provided, use defaults only 
    if nargin == 0 

        overrideParams = struct(); 

    else 

        % Convert name-value pairs into a struct 
        overrideParams = struct(varargin{:}); 

    end


    % Merge defaults with overrides (overrides take precedence)
    paramsOut = mergeStructs(defaultParams, overrideParams);


    % ====== CHEMISTRY-SPECIFIC ADJUSTMENTS ======
    % Adjust RBP lifetime based on selected chemistry
    chem        = string(paramsOut.Chemistry);
    V           = paramsOut.V_pack_nom;

    % Define coefficients and bounds per chemistry
    switch chem

        case "LFP"
            chi_mu_perc = 1.54 * log(V) + 1.06;
            if V==800
                bounds = [9.4439, 12.8275];
            else 
                error("Unsupported V_pack_nom for LFP");
            end

        case "NMC"
            chi_mu_perc = 4.06 * log(V) - 1.91;
            if V==800
                bounds = [20.7396, 28.9010];
            else 
                error("Unsupported V_pack_nom for NMC");
            end

        otherwise
            error("Unknown Chemistry: %s", chem);

    end

    % Select chi%: mean (deterministic) or uniform sample (optionally seeded)
    if paramsOut.DefaultIndicator
        chi_percent = chi_mu_perc;                       % use chemistry mean
    else
        % if hadSeed 
        %     rng(randomSeed, 'twister'); 
        % end      % reproducible sampling if requested
        chi_value_perc = bounds(1) + (bounds(2)-bounds(1)) * rand; % U[a,b]
        chi_percent    = chi_value_perc;
    end

    % Apply lifetime multiplier to obtain Y_RBP
    paramsOut.Y_RBP = paramsOut.Y_CBP * (1 + chi_percent/100);

end



% -------------------------------------------------------------------------
% MERGESTRUCTS
% -------------------------------------------------------------------------
function paramsOut = mergeStructs(defaultParams, overrideParams)
%MERGESTRUCTS Combine two structs with precedence to overrides.
%   OUT = MERGESTRUCTS(DEFAULTS, OVERRIDES) copies DEFAULTS and replaces any
%   fields present in OVERRIDES.

    % Initialize the output with all default parameters
    paramsOut = defaultParams;                           

    % Get the names of all fields in the overrideParams struct
    f = fieldnames(overrideParams);

    % Loop through each field in overrideParams
    for i = 1:numel(f)
        paramsOut.(f{i}) = overrideParams.(f{i});        
    end

end



% -------------------------------------------------------------------------
% GETPARAMRANGES
% -------------------------------------------------------------------------
function ranges = getParamRanges()
%GETPARAMRANGES Define allowable ranges for key economic/model parameters.
%   R = GETPARAMRANGES() returns a struct R whose fields encode either
%   numeric LIMITS = [min max] for scalar parameters or enumerated LIST
%   values for categorical parameters. 
%
%   FIELDS (scalar parameters use .LIMITS = [min max])
%     R.Y_EV.LIMITS                      Vehicle lifetime [years].
%     R.r_CBP.LIMITS                     Discount rate for CBP [fraction].
%     R.alpha_CBP.LIMITS                 Annual O&M rate for CBP [fraction/yr].
%     R.delta_loss.LIMITS                Relative increase in ohmic losses (RBP) [fraction].
%     R.L.LIMITS                         Annual driving distance [km/year].
%     R.nu.LIMITS                        Upfront cost factor (RBP vs CBP) [fraction].
%     R.E_pack_nom.LIMITS                Nominal pack energy [kWh].
%     R.E_pack_nom_detailed              Vector grid [kWh] for detailed sweeps.
%     R.n_E                              Number of grid points in E_pack_nom_detailed.
%     R.SOH_EOL_SecondLife_CBP.LIMITS    EOL SOH threshold for CBP (second life) [fraction].
%     R.SOH_EOL_SecondLife_RBP.LIMITS    EOL SOH threshold for RBP (second life) [fraction].
%
%   CATEGORICAL
%     R.Chemistry.LIST                   Allowed chemistries: "LFP", "NMC".
%     R.n_chem                           Number of chemistries in the list.
%
%   NOTES
%   - All LIMITS are inclusive and assume min <= max.
%
%   EXAMPLE
%     R = getParamRanges();
%     E_min = R.E_pack_nom.LIMITS(1);
%     allowed_chems = R.Chemistry.LIST;
%
%   SEE ALSO: GETPARAMS

    % ====== Numeric ranges (inclusive [min max]) ======
    % Lifetime [years]
    ranges.Y_EV.LIMITS                  = [15, 20];

    % Discount rates [fraction]
    ranges.r_CBP.LIMITS                 = [0.02, 0.04];

    % O&M rates [fraction/year]
    ranges.alpha_CBP.LIMITS             = [0.01, 0.03];

    % Relative energy/ohmic loss [fraction]
    ranges.delta_loss.LIMITS            = [0.01, 0.05];

    % Driving distance [km/year]
    ranges.L.LIMITS                     = [5000, 70000];

    % Upfront cost factor (RBP vs CBP) [fraction]
    ranges.nu.LIMITS                    = [0.01, 0.15];

    % Nominal pack energy [kWh]
    ranges.E_pack_nom.LIMITS            = [20, 120];
    % Grid for sweeps
    ranges.E_pack_nom_detailed          = 20:20:120;        
    ranges.n_E                          = numel(ranges.E_pack_nom_detailed);

    % Second-life SOH thresholds [fraction]
    ranges.SOH_EOL_SecondLife_CBP.LIMITS  = [0.50, 0.60];
    ranges.SOH_EOL_SecondLife_RBP.LIMITS  = [0.40, 0.50];


    % ====== Categorical ranges ======
    % Allowed chemistries
    ranges.Chemistry.LIST               = ["LFP", "NMC"];
    ranges.n_chem                       = numel(ranges.Chemistry.LIST);


    % ====== Lightweight sanity checks (guard against inverted limits) ======
    assertLimits(ranges.Y_EV.LIMITS, 'Y_EV');
    assertLimits(ranges.r_CBP.LIMITS, 'r_CBP');
    assertLimits(ranges.alpha_CBP.LIMITS, 'alpha_CBP');
    assertLimits(ranges.delta_loss.LIMITS, 'delta_loss');
    assertLimits(ranges.L.LIMITS, 'L');
    assertLimits(ranges.nu.LIMITS, 'nu');
    assertLimits(ranges.E_pack_nom.LIMITS, 'E_pack_nom');
    assertLimits(ranges.SOH_EOL_SecondLife_CBP.LIMITS, 'SOH_EOL_SecondLife_CBP');
    assertLimits(ranges.SOH_EOL_SecondLife_RBP.LIMITS, 'SOH_EOL_SecondLife_RBP');

end



% -------------------------------------------------------------------------
% ASSERTLIMITS
% -------------------------------------------------------------------------
function assertLimits(lim, name)
%ASSERTLIMITS Ensure LIMITS is a 1x2 finite, ordered numeric vector.
    if ~(isnumeric(lim) && isvector(lim) && numel(lim)==2 && all(isfinite(lim)))
        error('getParamRanges:InvalidLimits', ...
              'Field "%s.LIMITS" must be a 1x2 finite numeric vector.', name);
    end
    if lim(1) > lim(2)
        error('getParamRanges:LimitsOrder', ...
              'Field "%s.LIMITS" must satisfy min <= max.', name);
    end
end



% -------------------------------------------------------------------------
% COMPUTENPC_DISCANDNONDSC
% -------------------------------------------------------------------------
function results = computeNPC_DiscAndNonDisc(params)
%COMPUTENPC_DISCANDNONDSC Net present costs (discounted & non-discounted) for CBP vs RBP.
%   RESULTS = COMPUTENPC_DISCANDNONDSC(PARAMS) computes annual and total costs for a
%   Conventional Battery Pack (CBP) and a Reconfigurable Battery Pack (RBP),
%   returning both discounted and non-discounted aggregates, the residual value streams,
%   and the net present cost (NPC).
%
%   REQUIRED FIELDS IN PARAMS (units)
%     Y_EV [yr], Y_CBP [yr], Y_RBP [yr]
%     E_pack_nom [kWh], c_pack_USDperkWh [USD/kWh], currExch_USDtoEUR [-]
%     c_en_EURperkWh [EUR/kWh], eta_en_kWhkm [kWh/km], L [km/yr]
%     r_CBP [-], r_RBP [-], alpha_CBP [-/yr], alpha_RBP [-/yr]
%     delta_loss [-], nu [-]
%     SOH_EOL_FirstLife_CBP, SOH_EOL_FirstLife_RBP, SOH_EOL_SecondLife_CBP, SOH_EOL_SecondLife_RBP [-]
%     c_res_depreccoeff [-]  (used inside computeResidual_DiscAndNonDisc)
%
%   OUTPUT (struct)
%     .years                : 1:ceil(Y_EV)
%     .CBP / .RBP:
%         .costs.InitCost
%         .costs.AnnualCosts.discounted / .nondiscounted (struct with fields: energy,oandm,replacement,residual,totalPerYear)
%         .costs.TotalAnnualCosts.discounted / .nondiscounted
%         .costs.RV_discounted / .RV_nondiscounted (totals)
%         .costs.NPC
%         .SOH_EOL_Pack, .SOH_residual, .E_residual, .residualStruct (details per replacement)
%     .delta_NPC            : NPC_CBP - NPC_RBP
%
%   NOTES
%   - Discounting uses per-year factors for annual costs and exact (fractional) timing
%     for replacement costs and residual values (see computeAnnualCosts_DiscAndNonDisc).
%   - Replacement costs rely on GETREPLCOST(t, E_pack_nom), which must return the capital
%     cost (in EUR) of a new pack at time t [yr]. This function is external to this file.

    % ---- Minimal input checks (helpful errors if something is missing) ----
    mustHave(params, {'Y_EV','Y_CBP','Y_RBP','E_pack_nom','c_pack_USDperkWh','currExch_USDtoEUR', ...
                      'c_en_EURperkWh','eta_en_kWhkm','L','r_CBP','r_RBP','alpha_CBP','alpha_RBP', ...
                      'delta_loss','nu','SOH_EOL_FirstLife_CBP','SOH_EOL_FirstLife_RBP', ...
                      'SOH_EOL_SecondLife_CBP','SOH_EOL_SecondLife_RBP','c_res_depreccoeff'});

    % Extract and precompute
    Y_EV        = params.Y_EV;
    Y_EV_floor  = floor(Y_EV);
    Y_EV_ceil   = ceil(Y_EV);
    d_Y_EV      = Y_EV - Y_EV_floor;
    Y_EV_array  = 1:1:Y_EV_ceil;


    % ========== UPFRONT COST ==========
    Cost_Init_CBP = params.E_pack_nom * params.c_pack_USDperkWh * params.currExch_USDtoEUR; 
    Cost_Init_RBP = (1 + params.nu) * Cost_Init_CBP;


    % ========== ENERGY COST ==========
    Cost_energy_CBP = params.c_en_EURperkWh * params.eta_en_kWhkm * params.L;             
    Cost_energy_RBP = (1 + params.delta_loss) * Cost_energy_CBP;


    % ========== RESIDUAL VALUE STREAMS ==========
    % Compute residual value at each retirement (discounted & non-discounted).
    Cost_ResidualValue_CBP = computeResidual_DiscAndNonDisc( ...
        params, ...
        params.Y_CBP, ...
        params.r_CBP, ...
        Cost_Init_CBP, ...
        params.SOH_EOL_FirstLife_CBP, ...
        params.SOH_EOL_SecondLife_CBP);

    Cost_ResidualValue_RBP = computeResidual_DiscAndNonDisc( ...
        params, ...
        params.Y_RBP, ...
        params.r_RBP, ...
        Cost_Init_RBP, ...
        params.SOH_EOL_FirstLife_RBP, ...
        params.SOH_EOL_SecondLife_RBP);


    % === CBP annual costs ===
    cbpInputs = struct( ...
        'Y_EV', Y_EV, ...
        'Y_EV_floor', Y_EV_floor, ...
        'd_Y_EV', d_Y_EV, ...
        'Y_battPack', params.Y_CBP, ...
        'initCost', Cost_Init_CBP, ...
        'alphaOM', params.alpha_CBP, ...
        'energyCost', Cost_energy_CBP, ...
        'discountRate', params.r_CBP, ...
        'E_pack_nom', params.E_pack_nom, ...
        'RV_discounted', Cost_ResidualValue_CBP.RV_discounted, ...
        'RV_nondiscounted', Cost_ResidualValue_CBP.RV_nondiscounted, ...
        't_residuals', Cost_ResidualValue_CBP.t_residuals, ...
        'nu', 0); 
    [Cost_Yearly_PerYear_CBP, Cost_Yearly_Total_CBP_d, Cost_Yearly_Total_CBP_nd] = computeAnnualCosts_DiscAndNonDisc(cbpInputs);


    % === RBP annual costs ===
    rbpInputs = struct( ...
        'Y_EV', Y_EV, ...
        'Y_EV_floor', Y_EV_floor, ...
        'd_Y_EV', d_Y_EV, ...
        'Y_battPack', params.Y_RBP, ...
        'initCost', Cost_Init_RBP, ...
        'alphaOM', params.alpha_RBP, ...
        'energyCost', Cost_energy_RBP, ...
        'discountRate', params.r_RBP, ...
        'E_pack_nom', params.E_pack_nom, ...
        'RV_discounted', Cost_ResidualValue_RBP.RV_discounted, ...
        'RV_nondiscounted', Cost_ResidualValue_RBP.RV_nondiscounted, ...
        't_residuals', Cost_ResidualValue_RBP.t_residuals, ...
        'nu', params.nu);
    [Cost_Yearly_PerYear_RBP, Cost_Yearly_Total_RBP_d, Cost_Yearly_Total_RBP_nd] = computeAnnualCosts_DiscAndNonDisc(rbpInputs);


    % === Net Present Cost (discounted totals with residuals subtracted) ===
    NPC_CBP   = Cost_Init_CBP + Cost_Yearly_Total_CBP_d - Cost_ResidualValue_CBP.RV_total_discounted;
    NPC_RBP   = Cost_Init_RBP + Cost_Yearly_Total_RBP_d - Cost_ResidualValue_RBP.RV_total_discounted;
    delta_NPC = NPC_CBP - NPC_RBP;


    % === Store results ===
    results.years = Y_EV_array;

    % CBP
    results.CBP.costs.InitCost                         = Cost_Init_CBP;
    results.CBP.costs.AnnualCosts.discounted           = Cost_Yearly_PerYear_CBP.discounted;
    results.CBP.costs.AnnualCosts.nondiscounted        = Cost_Yearly_PerYear_CBP.nondiscounted;
    results.CBP.costs.TotalAnnualCosts.discounted      = Cost_Yearly_Total_CBP_d;
    results.CBP.costs.TotalAnnualCosts.nondiscounted   = Cost_Yearly_Total_CBP_nd;
    results.CBP.costs.RV_discounted                    = Cost_ResidualValue_CBP.RV_total_discounted;
    results.CBP.costs.RV_nondiscounted                 = Cost_ResidualValue_CBP.RV_total_nondiscounted;
    results.CBP.costs.NPC                              = NPC_CBP;
    results.CBP.SOH_EOL_Pack                           = Cost_ResidualValue_CBP.SOH_EOL_Pack;
    results.CBP.SOH_residual                           = Cost_ResidualValue_CBP.SOH_residual;
    results.CBP.residualStruct                         = Cost_ResidualValue_CBP;

    % RBP
    results.RBP.costs.InitCost                         = Cost_Init_RBP;
    results.RBP.costs.AnnualCosts.discounted           = Cost_Yearly_PerYear_RBP.discounted;
    results.RBP.costs.AnnualCosts.nondiscounted        = Cost_Yearly_PerYear_RBP.nondiscounted;
    results.RBP.costs.TotalAnnualCosts.discounted      = Cost_Yearly_Total_RBP_d;
    results.RBP.costs.TotalAnnualCosts.nondiscounted   = Cost_Yearly_Total_RBP_nd;
    results.RBP.costs.RV_discounted                    = Cost_ResidualValue_RBP.RV_total_discounted;
    results.RBP.costs.RV_nondiscounted                 = Cost_ResidualValue_RBP.RV_total_nondiscounted;
    results.RBP.costs.NPC                              = NPC_RBP;
    results.RBP.SOH_EOL_Pack                           = Cost_ResidualValue_RBP.SOH_EOL_Pack;
    results.RBP.SOH_residual                           = Cost_ResidualValue_RBP.SOH_residual;
    results.RBP.residualStruct                         = Cost_ResidualValue_RBP;

    % Delta
    results.delta_NPC                                  = delta_NPC;

end



% -------------------------------------------------------------------------
% COMPUTERESIDUAL_DISCANDNONDSC
% -------------------------------------------------------------------------
function output = computeResidual_DiscAndNonDisc(params, Y_battPack, r_battPack, Cost_Init_battPack, SOH_EOL_FirstLife, SOH_EOL_SecondLife)
%COMPUTERESIDUAL_DISCANDNONDSC Residual value streams at pack retirements.
%   OUT = COMPUTERESIDUAL_DISCANDNONDSC(PARAMS, Y_BP, r_BP, C0_BP, SOH1, SOH2)
%   computes residual values (discounted and non-discounted) at each battery
%   pack retirement, including second-life potential. Residual value at a
%   retirement time t accounts for (i) remaining SOH beyond SOH2 and (ii) a
%   depreciation coefficient (PARAMS.c_res_depreccoeff) applied to the pack’s
%   replacement capital cost at that time.
%
%   INPUTS
%     params               struct containing at least: Y_EV [yr], E_pack_nom [kWh], c_res_depreccoeff [-].
%     Y_battPack           scalar [yr] pack life.
%     r_battPack           scalar discount rate (fraction).
%     Cost_Init_battPack   scalar [EUR] initial capital cost at t=0 for this pack type.
%     SOH_EOL_FirstLife    scalar [-] EOL threshold for first-life (e.g., 0.8).
%     SOH_EOL_SecondLife   scalar [-] EOL threshold for second-life (e.g., 0.5).
%
%   OUTPUT (struct)
%     .RV_total_discounted, .RV_total_nondiscounted
%     .RV_discounted, .RV_nondiscounted        (1×N arrays per retirement)
%     .t_residuals                             (1×N retirement times [yr])
%     .SOH_EOL_Pack, .SOH_residual, 
%
%   EXTERNAL DEPENDENCY
%     GETREPLCOST(t, E_nom) must return the capital cost [EUR] of a new pack
%     at time t [yr] for nominal energy E_nom [kWh]. 

    % Basic guards
    mustHave(params, {'Y_EV','E_pack_nom','c_res_depreccoeff'});
    Y_EV  = params.Y_EV;
    E_nom = params.E_pack_nom;

    % Replacement schedule (retirement times). Include final partial pack at Y_EV.
    max_n     = floor(Y_EV / Y_battPack);
    t_replace = Y_battPack * (1:max_n);
    if mod(Y_EV, Y_battPack) > 0 || isempty(t_replace)
        t_replace = [t_replace, Y_EV];
    end

    n_packs      = numel(t_replace);
    t_residuals  = t_replace;

    SOH_EOL_Pack     = zeros(1, n_packs);
    SOH_residual     = zeros(1, n_packs);
    RV_nondiscounted = zeros(1, n_packs);
    RV_discounted    = zeros(1, n_packs);

    % Iterate through each retired pack, compute residual value
    for i = 1:n_packs

        % Years elapsed for the pack that retires at t_replace(i)
        if i == 1
            y_elapsed = t_replace(i);
        else
            y_elapsed = t_replace(i) - t_replace(i-1);
        end

        % Linear SOH trajectory from 1 to SOH_EOL_FirstLife over Y_battPack
        SOH_EOL_Pack(i) = 1 - (y_elapsed / Y_battPack) * (1 - SOH_EOL_FirstLife);

        % Fraction of second-life potential realized (normalized to [0,1])
        SOH_residual(i) = (SOH_EOL_Pack(i) - SOH_EOL_SecondLife) / (1 - SOH_EOL_SecondLife);
        %SOH_residual(i) = max(0, min(1, SOH_residual(i))); % clip to [0,1] for safety

        % Residual value is proportional to contemporary capital cost times a depreciation coeff
        if (i == 1)
            RV_nondiscounted(i) = SOH_residual(i) * Cost_Init_battPack * params.c_res_depreccoeff;
            Cost_Init_battPack_new  = getReplCost(t_replace(i), E_nom);
        else
            RV_nondiscounted(i)  = SOH_residual(i) * Cost_Init_battPack_new * params.c_res_depreccoeff;            
        end
        RV_discounted(i)    = RV_nondiscounted(i) / ((1 + r_battPack)^t_replace(i));
    end

    % Outputs
    output.RV_total_discounted     = sum(RV_discounted);
    output.RV_total_nondiscounted  = sum(RV_nondiscounted);
    output.RV_discounted           = RV_discounted;
    output.RV_nondiscounted        = RV_nondiscounted;
    output.t_residuals             = t_residuals;
    output.SOH_EOL_Pack            = SOH_EOL_Pack;
    output.SOH_residual            = SOH_residual;
end



% -------------------------------------------------------------------------
% GETREPLCOST
% -------------------------------------------------------------------------
function replCost_EUR = getReplCost(t, E_pack_nom)
%GETREPLCOST Projected battery pack replacement cost in EUR at time t.
%   C = GETREPLCOST(t, E_pack_nom) returns the projected replacement cost C
%   (in EUR) for a battery pack with nominal energy E_pack_nom [kWh], at a
%   time t [years] after the simulation start year 2025.
%
%   The model assumes a USD/kWh price trajectory following an exponential
%   decay toward an asymptote:
%       price_USDperkWh(year) = p1 * exp(-p2 * (year - t0)) + p3
%   where year = 2025 + t, with parameters:
%       p1 = 172.0648,  p2 = 0.1420,  p3 = 63.0950,  t0 = 2017.
%   The total pack cost in USD is price_USDperkWh * E_pack_nom, then
%   converted to EUR with a fixed FX rate (USD→EUR) of 0.8554.
%
%   INPUTS
%     t           numeric scalar/array
%                 Time after 2025 in years (t >= 0). Supports vector input.
%     E_pack_nom  (1x1) double, kWh
%                 Nominal pack energy; must be positive.
%
%   OUTPUT
%     replCost_EUR  numeric scalar/array, EUR
%                   Replacement cost(s) corresponding to t (same size as t).
%
%   NOTES
%   - This function is deterministic given the constants below.
%   - The exchange rate is fixed for simplicity.
%   - Vectorized over t for efficient evaluation of multiple replacement times.
%
%   EXAMPLES
%     % Single replacement 6.5 years after 2025 (i.e., mid-2031)
%     c = getReplCost(6.5, 80);
%
%     % Vector of replacement times
%     t = [5, 10, 12.5];
%     c = getReplCost(t, 60);

    % ---- Input validation (lightweight, compatible with older MATLAB) ----
    validateattributes(t, {'numeric'}, {'nonempty','real','finite','>=',0}, mfilename, 't', 1);
    validateattributes(E_pack_nom, {'double','single'}, {'scalar','real','finite','>',0}, mfilename, 'E_pack_nom', 2);

    % ---- Price trajectory parameters (USD/kWh) ----
    p1   = 172.0648;
    p2   = 0.1420;
    p3   = 63.0950;
    t0   = 2017;
    yr0  = 2025;         % simulation start year
    FX   = 0.8554;       % USD -> EUR

    % ---- Compute USD/kWh at calendar year (vectorized over t) ----
    year                    = yr0 + t;
    replCost_USDperkWh      = p1 .* exp(-p2 .* (year - t0)) + p3;

    % ---- Convert to total USD, then to EUR ----
    replCost_USD            = replCost_USDperkWh .* E_pack_nom;
    replCost_EUR            = FX .* replCost_USD;
end



% -------------------------------------------------------------------------
% COMPUTEANNUALCOSTS_DISCANDNONDSC
% -------------------------------------------------------------------------
function [costs, totalCost, totalCost_nondiscounted, replacementSchedule] = computeAnnualCosts_DiscAndNonDisc(inputs)
%COMPUTEANNUALCOSTS_DISCANDNONDSC Annual cost streams (discounted & non-discounted).
%   [COSTS, TOTAL_D, TOTAL_ND, REPL] = COMPUTEANNUALCOSTS_DISCANDNONDSC(INPUTS)
%   builds per-year cost arrays for energy, O&M, replacement, and residual value,
%   then returns discounted and non-discounted totals.
%
%   REQUIRED INPUT FIELDS (units)
%     Y_EV [yr], Y_EV_floor [yr], d_Y_EV [- in (0,1)]
%     Y_battPack [yr], energyCost [EUR/yr], alphaOM [-/yr], discountRate [-]
%     initCost [EUR], E_pack_nom [kWh], nu [-]
%   OPTIONAL (if residuals should be allocated to year bins)
%     RV_discounted (1×N), RV_nondiscounted (1×N), t_residuals (1×N [yr])
%
%   OUTPUTS
%     costs.discounted / .nondiscounted : struct with fields
%         .energy, .oandm, .replacement, .residual, .totalPerYear (1×Y_max)
%     totalCost                         : sum(costs.discounted.totalPerYear)
%     totalCost_nondiscounted           : sum(costs.nondiscounted.totalPerYear)
%     replacementSchedule               : struct with t_replace, replacementYears, frac_old, frac_new
%
%   NOTES
%   - Replacement costs use fractional-year timing for discounting:
%       t_replace = Y_battPack * (1:max_n), discounted by (1+r)^t_replace.
%   - The capital cost of a replacement is GETREPLCOST(t, E_pack_nom) * (1+nu).
%   - Year Y_EV_floor+1 is included only when d_Y_EV>0, weighted by d_Y_EV.

    % Guards
    mustHave(inputs, {'Y_EV','Y_EV_floor','d_Y_EV','Y_battPack','energyCost','alphaOM', ...
                      'discountRate','initCost','E_pack_nom','nu'});

    % Unpack
    Y_EV         = inputs.Y_EV;
    Y_EV_floor   = inputs.Y_EV_floor;
    d_Y_EV       = inputs.d_Y_EV;
    Y_battPack   = inputs.Y_battPack;
    energyCost   = inputs.energyCost;
    alphaOM      = inputs.alphaOM;
    discountRate = inputs.discountRate;
    nu           = inputs.nu;
    initCost     = inputs.initCost;
    E_pack_nom   = inputs.E_pack_nom;

    % Residuals (optional)
    if ((isfield(inputs, 'RV_discounted')) && (isfield(inputs,'RV_nondiscounted')) && (isfield(inputs,'t_residuals')))
        
        RV_discounted    = inputs.RV_discounted;
        RV_nondiscounted = inputs.RV_nondiscounted;
        t_residuals      = inputs.t_residuals;
        
        % Assign to calendar years
        RV_years         = ceil(t_residuals);  

    else

        RV_discounted    = [];
        RV_nondiscounted = [];
        RV_years         = [];

    end

    % Replacement schedule and fractions for O&M blending
    max_n         = floor(Y_EV / Y_battPack);
    t_replace     = Y_battPack * (1:max_n);
    replacementYears = ceil(t_replace);
    frac_old      = mod(t_replace, 1);
    frac_new      = 1 - frac_old;

    replacementSchedule.t_replace        = t_replace;
    replacementSchedule.replacementYears  = replacementYears;
    replacementSchedule.frac_old          = frac_old;
    replacementSchedule.frac_new          = frac_new;

    % Number of calendar years included
    Y_max = Y_EV_floor + (d_Y_EV > 0);

    % Preallocate discounted (D) and non-discounted (ND) arrays
    D.energy = zeros(1, Y_max);
    D.oandm = zeros(1, Y_max);
    D.replacement = zeros(1, Y_max);
    D.residual = zeros(1, Y_max);
    ND = D;

    % Track capital cost level for O&M scaling within each year
    currentCapitalCost = initCost;

    % Year loop
    for year = 1:Y_max

        discountFactor = (1 + discountRate)^year;

        % Energy cost
        if year <= Y_EV_floor

            yearWeight = 1.0;

        else

            % Partial last year
            yearWeight = d_Y_EV; 

        end
        D.energy(year)  = yearWeight * energyCost / discountFactor;
        ND.energy(year) = yearWeight * energyCost;

        % Replacement event happening in this calendar year?
        idx_repl = find(replacementYears == year, 1, 'first');
        if (~isempty(idx_repl))

            % Fractional time of replacement
            actualTime = t_replace(idx_repl);                   
            discountFactor_actual = (1 + discountRate)^actualTime;

            newCapitalCost = getReplCost(actualTime, E_pack_nom) * (1 + nu); 

            % O&M blends the old and new pack costs within the same calendar year
            omCost = yearWeight * ( ...
                frac_old(idx_repl) * alphaOM * currentCapitalCost + ...
                frac_new(idx_repl) * alphaOM * newCapitalCost ...
            );

            D.oandm(year)        = omCost / discountFactor;
            ND.oandm(year)       = omCost;
            D.replacement(year)  = newCapitalCost / discountFactor_actual;
            ND.replacement(year) = newCapitalCost;

            % Update capital cost baseline for subsequent years
            currentCapitalCost = newCapitalCost;

        else

            % No replacement in this year: O&M based on current pack
            omCost = yearWeight * alphaOM * currentCapitalCost;
            D.oandm(year)  = omCost / discountFactor;
            ND.oandm(year) = omCost;

        end

        % Residual value(s) realized in this calendar year: subtract from costs
        if ~isempty(RV_years)

            idx_rv = find(RV_years == year);
            if (~isempty(idx_rv))
                D.residual(year)  = -sum(RV_discounted(idx_rv));
                ND.residual(year) = -sum(RV_nondiscounted(idx_rv));
            end

        end

    end

    % Aggregate
    D.totalPerYear  = D.energy + D.oandm + D.replacement + D.residual;
    ND.totalPerYear = ND.energy + ND.oandm + ND.replacement + ND.residual;

    costs.discounted    = D;
    costs.nondiscounted = ND;

    totalCost                = sum(D.totalPerYear);
    totalCost_nondiscounted  = sum(ND.totalPerYear);

end



% -------------------------------------------------------------------------
% SAMPLEDELTANPC_LHS
% -------------------------------------------------------------------------
function [samples_out, delta_NPC_out, results_out, p_out] = sampleDeltaNPC_LHS(ranges, n_samples, params_baseline)
%SAMPLEDELTANPC_LHS Latin Hypercube Sampling for ΔNPC (CBP vs RBP) experiments.
%   [SAMPLES, DELTA_NPC, RESULTS, P] = SAMPLEDeltANPC_LHS(RANGES, N, BASELINE)
%   generates N−1 Latin Hypercube samples across predefined parameter ranges,
%   appends the BASELINE parameter set as the Nth sample, then evaluates the
%   discounted/non-discounted NPC model for each sample.
%
%   INPUTS
%     ranges          struct
%                     Produced by GETPARAMRANGES. Must contain fields:
%                       .Y_EV.LIMITS, .delta_loss.LIMITS, .L.LIMITS, .nu.LIMITS,
%                       .r_CBP.LIMITS, .alpha_CBP.LIMITS,
%                       .E_pack_nom_detailed (vector) and .n_E,
%                       .Chemistry.LIST (string array) and .n_chem,
%                       .SOH_EOL_SecondLife_CBP.LIMITS, .SOH_EOL_SecondLife_RBP.LIMITS.
%     n_samples       scalar integer (>=1)
%                     Total number of samples returned (includes baseline as last row).
%     params_baseline struct
%                     Baseline parameter set (fields matching those sampled here).
%
%   OUTPUTS
%     samples_out     struct of sampled parameters (column vectors, length = n_samples).
%                     The final row corresponds to BASELINE, and includes:
%                         Y_EV, delta_loss, L, nu, r_CBP, r_RBP,
%                         alpha_CBP, alpha_RBP, E_pack_nom, Chemistry,
%                         SOH_EOL_SecondLife_CBP, SOH_EOL_SecondLife_RBP, DefaultIndicator
%                     (DefaultIndicator=false for LHS draws, true for baseline).
%     delta_NPC_out   (n_samples×1) double
%                     ΔNPC = NPC_CBP − NPC_RBP for each sample.
%     results_out     (n_samples×1) struct
%                     Full outputs of COMPUTENPC_DISCANDNONDSC per sample.
%     p_out           (n_samples×1) struct
%                     The resolved parameter structs returned by GETPARAMS.
%
%   NOTES
%   - Randomness: rng(42) is used for reproducibility.
%   - Sampling space: 10 dimensions (Y_EV, δ_loss, L, ν, r_CBP, α_CBP,
%     E_pack_nom (discrete), Chemistry (discrete), SOH_CBP^2nd, SOH_RBP^2nd).
%   - r_RBP is set equal to r_CBP, and α_RBP = 0.5·α_CBP.
%   - Requires Statistics and Machine Learning Toolbox for LHSDESIGN.
%
%   EXAMPLE
%     R = getParamRanges();
%     base = getParams();
%     [S, dNPC, Rz, P] = sampleDeltaNPC_LHS(R, 500, base);

    % ------------------------ Validation & setup -------------------------
    validateattributes(n_samples, {'numeric'}, {'scalar','integer','>=',1}, mfilename, 'n_samples', 2);
    mustHave(ranges, {'Y_EV','delta_loss','L','nu','r_CBP','alpha_CBP', ...
                      'E_pack_nom','E_pack_nom_detailed','n_E', ...
                      'Chemistry','n_chem', ...
                      'SOH_EOL_SecondLife_CBP','SOH_EOL_SecondLife_RBP'});

    % Reproducible design
    rng(42);

    % We sample N-1 points from the hypercube; the last one is the baseline.
    n_batch = n_samples - 1;

    % 10D LHS (see NOTES for dimensions), values in (0, 1)
    X = lhsdesign(n_batch, 10); 


    % ------------------------ Draw samples -------------------------------
    samples_all = struct();

    % Continuous directly from LIMITS
    samples_all.Y_EV       = interp1([0 1], ranges.Y_EV.LIMITS,       X(:,1));
    samples_all.delta_loss = interp1([0 1], ranges.delta_loss.LIMITS, X(:,2));
    samples_all.L          = interp1([0 1], ranges.L.LIMITS,          X(:,3));
    samples_all.nu         = interp1([0 1], ranges.nu.LIMITS,         X(:,4));

    % Rates: r_CBP sampled; r_RBP set equal
    samples_all.r_CBP      = interp1([0 1], ranges.r_CBP.LIMITS,      X(:,5));
    samples_all.r_RBP      = samples_all.r_CBP;

    % O&M: alpha_CBP sampled; alpha_RBP = 0.5*alpha_CBP
    samples_all.alpha_CBP  = interp1([0 1], ranges.alpha_CBP.LIMITS,  X(:,6));
    samples_all.alpha_RBP  = samples_all.alpha_CBP * 0.5;

    % Discrete: E_pack_nom via indexed grid
    n_E        = ranges.n_E;
    E_indices  = min(ceil(X(:,7) * n_E), n_E);
    samples_all.E_pack_nom = transpose(ranges.E_pack_nom_detailed(E_indices));

    % Discrete: Chemistry from allowed list
    n_chem       = ranges.n_chem;
    chem_indices = min(ceil(X(:,8) * n_chem), n_chem);
    samples_all.Chemistry = transpose(ranges.Chemistry.LIST(chem_indices));

    % SOH thresholds (second life) for CBP and RBP
    samples_all.SOH_EOL_SecondLife_CBP = interp1([0 1], ranges.SOH_EOL_SecondLife_CBP.LIMITS, X(:,9));
    samples_all.SOH_EOL_SecondLife_RBP = interp1([0 1], ranges.SOH_EOL_SecondLife_RBP.LIMITS, X(:,10));

    % Output struct (start with LHS draws)
    samples_out = samples_all;

    % DefaultIndicator flags: false for sampled points, true for baseline
    samples_out.DefaultIndicator = false(n_batch, 1);


    % ------------------------ Append baseline ---------------------------
    % Capture list of fields we created above (so we only fill those from baseline)
    fieldnames_list = fieldnames(samples_all);

    % Extend each field to include the baseline as the last entry    
    for k = 1:numel(fieldnames_list)

        field = fieldnames_list{k};

        % Ensure the vector exists for the batch; extend by one row for the baseline
        samples_out.(field)(n_samples, 1) = params_baseline.(field);

    end
    samples_out.DefaultIndicator(n_samples, 1) = true;
    

    % ------------------------ Evaluate the model ------------------------
    n_total        = n_samples;
    delta_NPC_out  = zeros(n_total, 1);

    % First sample initializes struct arrays for results/p
    p_first     = getParams(getSampleStruct(samples_out, 1));
    result_first = computeNPC_DiscAndNonDisc(p_first);
    delta_NPC_out(1) = result_first.delta_NPC;

    results_out = repmat(result_first, n_total, 1);
    p_out       = repmat(p_first, n_total, 1);

    % Remaining samples
    for i = 2:n_total

        p_i        = getParams(getSampleStruct(samples_out, i));
        result_i   = computeNPC_DiscAndNonDisc(p_i);

        p_out(i)         = p_i;
        results_out(i)   = result_i;
        delta_NPC_out(i) = result_i.delta_NPC;

    end

    % Run diagnostics:
    % paramNames = {'Y_{EV}','\delta_{loss}','L','\nu','r_{CBP}','r_{RBP}', ...
    %               '\alpha_{CBP}','\alpha_{RBP}','E_{pack,nom}', ...
    %               'SOH^{2nd}_{CBP}','SOH^{2nd}_{RBP}'};
    % runLHSdiagnostics(samples_out, delta_NPC_out, paramNames);

end



% -------------------------------------------------------------------------
% GETSAMPLESTRUCT
% -------------------------------------------------------------------------
function s = getSampleStruct(samples_out, i)
%GETSAMPLESTRUCT Extract a single parameter set from a sampled structure.
%   S = GETSAMPLESTRUCT(SAMPLES_OUT, I) returns a scalar struct S with
%   the parameter values corresponding to sample index I. This struct
%   can be passed directly to GETPARAMS to construct a complete parameter
%   input for computeNPC or computeNPC_DiscAndNonDisc.
%
%   INPUTS
%     samples_out   struct
%                   Output from sampleDeltaNPC_LHS. Each field is a
%                   column vector (size n_samples × 1) of sampled values.
%     i             scalar integer
%                   Index of the sample to extract (1 ≤ i ≤ n_samples).
%
%   OUTPUT
%     s             scalar struct with fields:
%       .Y_EV, .alpha_CBP, .alpha_RBP, .delta_loss
%       .r_CBP, .r_RBP, .L, .nu
%       .E_pack_nom, .Chemistry
%
%   EXAMPLE
%     s1 = getSampleStruct(samples_out, 1);
%     params = getParams(s1);

    % --- Basic validation ---
    validateattributes(i, {'numeric'}, ...
        {'scalar','integer','>=',1,'<=',numel(samples_out.Y_EV)}, ...
        mfilename, 'i', 2);

    % --- Extract fields (names must match those expected by getParams) ---
    s.Y_EV       = samples_out.Y_EV(i);
    s.alpha_CBP  = samples_out.alpha_CBP(i);
    s.alpha_RBP  = samples_out.alpha_RBP(i);
    s.delta_loss = samples_out.delta_loss(i);
    s.r_CBP      = samples_out.r_CBP(i);
    s.r_RBP      = samples_out.r_RBP(i);
    s.L          = samples_out.L(i);
    s.nu         = samples_out.nu(i);
    s.E_pack_nom = samples_out.E_pack_nom(i);
    s.Chemistry  = samples_out.Chemistry(i);

end



% -------------------------------------------------------------------------
% PLOTMULTIVARIATESENSITIVITYPANELS
% -------------------------------------------------------------------------
function plotAnnualCostBreakdown_DiscNonDisc_Subplot(ax, results, yLim, yTicks, costType)
%PLOTANNUALCOSTBREAKDOWN_DISCNONDISC_SUBPLOT Stacked per-year cost bars (CBP vs RBP).
%   PLOTANNUALCOSTBREAKDOWN_DISCNONDISC_SUBPLOT(AX, RESULTS, YLIM, YTICKS, COSTTYPE)
%   draws a stacked bar chart comparing five cost categories for CBP and RBP
%   in each calendar year: Energy, O&M, Upfront (at t=0), Replacement, Residual.
%   Values are plotted in kilo-euros (kEUR).
%
%   INPUTS
%     ax        matlab.graphics.axis.Axes
%               Target axes (subplot) to draw into.
%     results   struct
%               Output from COMPUTENPC_DISCANDNONDSC with fields:
%                 .years
%                 .CBP.costs.InitCost
%                 .CBP.costs.AnnualCosts.discounted / .nondiscounted.{energy,oandm,replacement,residual}
%                 .RBP.costs.InitCost
%                 .RBP.costs.AnnualCosts.discounted / .nondiscounted.{energy,oandm,replacement,residual}
%     yLim      1×2 double | []   (optional) y-axis limits in kEUR.
%     yTicks    1×N double | []   (optional) y-tick positions in kEUR.
%     costType  char|string        'discounted' (default) or 'nondiscounted'.
%
%   NOTES
%   - The “Upfront” bar appears at t=0 only (first pair of bars).
%   - CBP and RBP bars at each tick are offset for visual comparison.
%   - Residual values are included (negative bars when subtracting costs).

    % -------- Defaults & validation --------
    if nargin < 5 || isempty(costType), costType = 'discounted'; end
    costType = validatestring(string(costType), ["discounted","nondiscounted"], mfilename, 'costType', 5);

    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('plotAnnualCost:InvalidAxes','First argument must be a valid axes handle.');
    end

    % Required fields quick check 
    if ~isfield(results,'CBP') || ~isfield(results,'RBP') || ~isfield(results,'years')
        error('plotAnnualCost:BadResults','RESULTS must contain .CBP, .RBP, and .years.');
    end


    % -------- Prepare data (convert to kEUR) --------
    years     = results.years;
    Y         = numel(years);      % number of calendar years shown
    Y_full    = Y + 1;             % include t=0 for upfront cost
    tickStep  = 1.5;               % horizontal spacing between years
    x_base    = tickStep * (0:Y);  % include t=0
    x_cbp     = x_base - 0.3;
    x_rbp     = x_base + 0.3;

    % Extract annual cost streams for selected type
    cbpCosts = results.CBP.costs.AnnualCosts.(costType);
    rbpCosts = results.RBP.costs.AnnualCosts.(costType);

    % Compose the 5 stacks (Energy, O&M, Upfront@t0, Replacement, Residual)
    % Prepend a zero for t=0 to align year bars; Upfront placed at index 1.
    CBP_e       = [0, cbpCosts.energy]        / 1000;
    CBP_m       = [0, cbpCosts.oandm]         / 1000;
    CBP_upfront = zeros(1, Y_full);  CBP_upfront(1) = results.CBP.costs.InitCost / 1000;
    CBP_r       = [0, cbpCosts.replacement]   / 1000;
    CBP_res     = [0, cbpCosts.residual]      / 1000;

    RBP_e       = [0, rbpCosts.energy]        / 1000;
    RBP_m       = [0, rbpCosts.oandm]         / 1000;
    RBP_upfront = zeros(1, Y_full);  RBP_upfront(1) = results.RBP.costs.InitCost / 1000;
    RBP_r       = [0, rbpCosts.replacement]   / 1000;
    RBP_res     = [0, rbpCosts.residual]      / 1000;

    % Matrix form: rows = bars, cols = stacks (Energy, O&M, Upfront, Replacement, Residual)
    data_cbp = [CBP_e', CBP_m', CBP_upfront', CBP_r', CBP_res'];
    data_rbp = [RBP_e', RBP_m', RBP_upfront', RBP_r', RBP_res'];


    % -------- Plot --------
    axes(ax); 
    hold(ax,'on');

    % Bars
    b_cbp = bar(ax, x_cbp, data_cbp, 0.4, 'stacked');
    b_rbp = bar(ax, x_rbp, data_rbp, 0.4, 'stacked');

    % Colors (consistent palette; RBP slightly transparent)
    cbp_colors = [ 0.6 0.8 1.0; 0.4 0.7 1.0; 0.2 0.5 1.0; 1.0 0.6 0.6; 0.6 0.4 0.2 ];
    rbp_colors = [ 0.6 0.9 0.6; 0.4 0.8 0.4; 0.2 0.6 0.2; 1.0 0.6 0.6; 0.6 0.4 0.2 ];

    for j = 1:5
        b_cbp(j).FaceColor = cbp_colors(j,:);  b_cbp(j).EdgeColor = 'none';
        b_rbp(j).FaceColor = rbp_colors(j,:);  b_rbp(j).EdgeColor = 'none';
        b_rbp(j).FaceAlpha = 0.7;
    end


    % -------- Axes styling --------
    if exist('yLim','var') && ~isempty(yLim), ylim(ax, yLim); end
    if exist('yTicks','var') && ~isempty(yTicks), yticks(ax, yTicks); end

    ylabel(ax, 'Cost [kEUR]', 'Interpreter','latex', 'FontName','Latin Modern Roman');
    xticks(ax, tickStep * (0:5:20));
    xticklabels(ax, string(0:5:20));
    xlabel(ax, 'Year [-]', 'Interpreter','latex', 'FontName','Latin Modern Roman');
    grid(ax,'on'); box(ax,'on');
    set(ax, 'FontSize',12, 'TickLabelInterpreter','latex', 'FontName','Latin Modern Roman');


    % -------- Legend (grouped CBP first, then RBP) --------
    handles = [ ...
        b_cbp(3), b_cbp(1), b_cbp(2), b_cbp(4), b_cbp(5), ...
        b_rbp(3), b_rbp(1), b_rbp(2), b_rbp(4), b_rbp(5) ];

    labels = { ...
        'CBP Upfront','CBP Energy','CBP O\&M','CBP Replacement','CBP Residual', ...
        'RBP Upfront','RBP Energy','RBP O\&M','RBP Replacement','RBP Residual' };

    lgd = legend(ax, handles, labels, ...
        'Orientation','horizontal', 'Location','northeast', ...
        'FontSize',9, 'Interpreter','latex');
    lgd.NumColumns    = 5;
    lgd.ItemTokenSize = [10 6];

    hold(ax,'off');

end



% -------------------------------------------------------------------------
% ADDPANELLABEL
% -------------------------------------------------------------------------
function addPanelLabel(ax, label, offsetX, offsetY, fontSize)
%ADDPANELLABEL Add subplot/panel label outside the top-left corner.
%
%   addPanelLabel(AX, LABEL) places LABEL (e.g. '\textbf{(a)}') near the
%   top-left corner of the axes AX, using LaTeX interpreter for publication-
%   quality figure panel annotations.
%
%   addPanelLabel(AX, LABEL, offsetX, offsetY, fontSize) allows specifying:
%
%     offsetX   Horizontal offset (normalized units, default = -0.12)
%     offsetY   Vertical offset (normalized units, default = 1.02)
%     fontSize  Label font size (default = 12)
%
%   INPUTS
%     ax        - Valid axes handle (subplot or standalone axes)
%     label     - Text to display (string or char, LaTeX compatible)
%
%   EXAMPLE
%     figure; ax = subplot(2,2,1);
%     plot(ax, 1:10, rand(1,10));
%     addPanelLabel(ax, '\textbf{(a)}');
%
%   NOTES
%   - Coordinates are in normalized units relative to the axes.
%   - LaTeX interpreter and Latin Modern Roman.
%   - Does not alter axis limits or layout spacing.
%
%   SEE ALSO: text, subplot, tiledlayout

    % -------- Input validation --------
    if nargin < 2
        error('addPanelLabel:MissingInputs', ...
            'At least two inputs (ax and label) are required.');
    end
    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('addPanelLabel:InvalidAxes', 'First argument must be an axes handle.');
    end
    if ~(ischar(label) || isstring(label))
        error('addPanelLabel:InvalidLabelType', 'Label must be a string or character array.');
    end

    % -------- Default values for offsets and size --------
    if nargin < 3 || isempty(offsetX), offsetX = -0.12; end
    if nargin < 4 || isempty(offsetY), offsetY = 1.02;  end
    if nargin < 5 || isempty(fontSize), fontSize = 12;  end

    % -------- Place the label --------
    text(ax, offsetX, offsetY, label, ...
        'Units', 'normalized', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'latex', ...
        'FontName', 'Latin Modern Roman', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');
end



% -------------------------------------------------------------------------
% PLOTCUMULATIVENPC_DISCNONDISC_SUBPLOT
% -------------------------------------------------------------------------
function plotCumulativeNPC_DiscNonDisc_Subplot(ax, results, yLabelOn, legendOn, legendLocation, label1, deltaNPC, yLim, yTicks, costType)
%PLOTCUMULATIVENPC_DISCNONDISC_SUBPLOT Cumulative NPC trajectories (CBP vs RBP).
%   PLOTCUMULATIVENPC_DISCNONDISC_SUBPLOT(AX, RESULTS, YLABELON, LEGENDON,
%   LEGENDLOCATION, LABEL1, DELTANPC, YLIM, YTICKS, COSTTYPE) plots the
%   cumulative sum of annual costs for CBP and RBP as kEUR time series
%   (Energy + O&M + Upfront + Replacement + Residual), optionally adding a
%   legend, axis labels, and a label box with ΔNPC and break-even year.
%
%   INPUTS
%     ax              axes handle (target subplot)
%     results         struct from COMPUTENPC_DISCANDNONDSC (uses .years and .CBP/.RBP.costs)
%     yLabelOn        logical  show y-label if true
%     legendOn        logical  show legend if true
%     legendLocation  char/string  e.g., 'northwest'
%     label1          char/string  scenario label printed in the corner
%     deltaNPC        double       ΔNPC to print (interpreted as EUR in the label)
%     yLim            1x2 double | []  optional y limits (kEUR)
%     yTicks          1xN double | []  optional y ticks (kEUR)
%     costType        'discounted' (default) or 'nondiscounted'
%
%   NOTES
%   - Curves are in kEUR.
%   - Break-even is linearly interpolated between adjacent samples where the
%     ΔNPC curve crosses zero.
%   - Year tick spacing uses 1.5 units per year to align with the panel grid.

    % ---- defaults & validation ----
    if nargin < 10 || isempty(costType), costType = 'discounted'; end
    costType = validatestring(string(costType), ["discounted","nondiscounted"], mfilename, 'costType', 10);

    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('plotCumulativeNPC:InvalidAxes','First argument must be a valid axes handle.');
    end
    if ~isfield(results,'years') || ~isfield(results,'CBP') || ~isfield(results,'RBP')
        error('plotCumulativeNPC:BadResults','RESULTS must contain fields: years, CBP, RBP.');
    end

    % ---- time base ----
    years      = results.years;
    Y          = numel(years);
    Y_full     = Y + 1;
    years_full = [0, years];
    tick_spacing = 1.5;
    x_base     = tick_spacing * (0:Y);  % includes t=0


    % ---- select annual cost streams ----
    cbp = results.CBP.costs.AnnualCosts.(costType);
    rbp = results.RBP.costs.AnnualCosts.(costType);


    % ---- build per-year totals (EUR), include t=0 upfront at index 1 ----
    CBP_e   = [0, cbp.energy];
    CBP_m   = [0, cbp.oandm];
    CBP_r   = [results.CBP.costs.InitCost, cbp.replacement];
    CBP_res = [0, cbp.residual];  % residual already signed

    RBP_e   = [0, rbp.energy];
    RBP_m   = [0, rbp.oandm];
    RBP_r   = [results.RBP.costs.InitCost, rbp.replacement];
    RBP_res = [0, rbp.residual];

    CBP_total = CBP_e + CBP_m + CBP_r + CBP_res;
    RBP_total = RBP_e + RBP_m + RBP_r + RBP_res;


    % ---- cumulative (kEUR) ----
    cum_CBP = cumsum(CBP_total) / 1000;
    cum_RBP = cumsum(RBP_total) / 1000;


    % ΔNPC curve in kEUR (for finding zero crossing only)
    delta_curve = cum_CBP - cum_RBP;


    % detect last sign change from + to − (closest to end)
    sgn = sign(delta_curve);
    idx_cross = find(diff(sgn) ~= 0, 1, 'last');


    % ---- plot ----
    axes(ax); 
    hold(ax,'on');

    hCBP = plot(ax, x_base, cum_CBP, '-', 'LineWidth',2, 'Color', [0.6, 0.8, 1]);
    hRBP = plot(ax, x_base, cum_RBP, '-', 'LineWidth',2, 'Color', [0.6, 0.9, 0.6]);

    % final markers
    hFinalCBP = plot(ax, x_base(end), cum_CBP(end), '*', 'MarkerSize',10, 'LineWidth',2, 'Color', [0 0.3 0.8]);
    hFinalRBP = plot(ax, x_base(end), cum_RBP(end), '*', 'MarkerSize',10, 'LineWidth',2, 'Color', [0 0.5 0]);

    % break-even (linear interpolation between years_full(idx_cross:idx_cross+1))
    haveBreakEven = ~isempty(idx_cross) && idx_cross < Y_full;
    if haveBreakEven
        x1 = years_full(idx_cross);   x2 = years_full(idx_cross+1);
        y1 = delta_curve(idx_cross);  y2 = delta_curve(idx_cross+1);
        if y2 ~= y1
            year_be = x1 - y1*(x2 - x1)/(y2 - y1); % in years
            % interpolate CBP/RBP values at that fractional year (kEUR)
            t_query = year_be * tick_spacing;
            CBP_be  = interp1(years_full*tick_spacing, cum_CBP, t_query, 'linear', 'extrap');
            %RBP_be  = interp1(years_full*tick_spacing, cum_RBP, t_query, 'linear', 'extrap');
            hBE1 = plot(ax, t_query, CBP_be, 'rx', 'MarkerSize',10, 'LineWidth',2);
            %hBE2 = plot(ax, t_query, RBP_be, 'rx', 'MarkerSize',10, 'LineWidth',2);
        else
            haveBreakEven = false; % avoid zero division edge case
        end
    end


    % ---- axes formatting ----
    if yLabelOn
        ylabel(ax, 'NPC [kEUR]', 'Interpreter','latex', 'FontName','Latin Modern Roman');
    end
    xlabel(ax, 'Year [-]', 'Interpreter','latex', 'FontName','Latin Modern Roman');
    xlim(ax, tick_spacing*[0 20]);
    xticks(ax, tick_spacing*(0:5:20));
    xticklabels(ax, string(0:5:20));
    if exist('yLim','var') && ~isempty(yLim),   ylim(ax, yLim);   end
    if exist('yTicks','var') && ~isempty(yTicks), yticks(ax, yTicks); end
    grid(ax,'on'); box(ax,'on');
    set(ax, 'FontSize',12, 'TickLabelInterpreter','latex', 'FontName','Latin Modern Roman');


    % ---- legend ----
    if legendOn
        legHandles = [hCBP, hRBP, hFinalCBP, hFinalRBP];
        legLabels  = {'CBP','RBP','Final CBP','Final RBP'};
        if haveBreakEven
            legHandles(end+1) = hBE1; 
            legLabels{end+1}  = 'Break-even'; 
        end
        lgd = legend(ax, legHandles, legLabels, ...
            'Location', legendLocation, ...
            'Interpreter','latex', 'FontSize',9);
        if contains(string(legendLocation), 'outside')
            lgd.NumColumns = numel(legLabels);
        end
        lgd.ItemTokenSize = [10 4];
    end


    % ---- label box (bottom-right) ----
    if nargin >= 6 && ~isempty(label1) && ~isempty(deltaNPC)
        if exist('year_be','var') == 1 && haveBreakEven
            breakEvenStr = sprintf('%.1f', year_be);
        else
            breakEvenStr = '--';
        end
        deltaNPCStr = sprintf('%.1f', deltaNPC); % printed in EUR (by design)
        label_text = sprintf('\\textbf{%s} scenario\n\\(\\Delta\\)NPC: %s EUR\nBreak-even: %s years', ...
                             label1, deltaNPCStr, breakEvenStr);
        text(ax, 0.98, 0.02, label_text, ...
             'Units','normalized', 'Interpreter','latex', ...
             'FontSize',9, 'HorizontalAlignment','right', ...
             'VerticalAlignment','bottom', 'FontName','Latin Modern Roman');
    end

    hold(ax,'off');

end



% -------------------------------------------------------------------------
% TORNADOPLOTDELTANPC_SUBPLOT
% -------------------------------------------------------------------------
function sorted_corrs = tornadoPlotDeltaNPC_Subplot(ax, samples, delta_NPC_kEuros)
%TORNADOPLOTDELTANPC_SUBPLOT Tornado (Spearman) sensitivities for ΔNPC.
%   SORTED_CORRS = TORNADOPLOTDELTANPC_SUBPLOT(AX, SAMPLES, DELTA_NPC_KEUROS)
%   computes Spearman rank correlations between ΔNPC (kEUR) and selected
%   input factors, sorts them by absolute magnitude (descending), and
%   renders a “tornado” bar plot in AX. The bars are labeled with ρ values.
%
%   INPUTS
%     ax                 Axes handle (subplot or tiledlayout axis)
%     samples            Struct with fields:
%                        Y_EV, alpha_CBP, alpha_RBP, delta_loss, nu, L, r_CBP, E_pack_nom
%     delta_NPC_kEuros   (n×1) or (1×n) numeric vector of ΔNPC in kEUR
%
%   OUTPUT
%     sorted_corrs       (k×1) double of |ρ| sorted descending (k = #factors)
%
%   NOTES
%   - We use Spearman (rank) correlation with pairwise handling of NaNs.
%   - Factors shown (LaTeX):
%       { Y_EV, Δα, δ_loss, ν, L, r, E_pack_nom }
%     where Δα = alpha_CBP − alpha_RBP and r = r_CBP.


    % ---- Validate inputs ----
    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('tornadoPlot:InvalidAxes','First input must be a valid axes handle.');
    end
    req = {'Y_EV','alpha_CBP','alpha_RBP','delta_loss','nu','L','r_CBP','E_pack_nom'};
    for k = 1:numel(req)
        if ~isfield(samples, req{k})
            error('tornadoPlot:MissingField','SAMPLES missing field "%s".', req{k});
        end
    end
    y = delta_NPC_kEuros(:);
    n = numel(y);

    % Ensure all fields are column vectors of length n
    for k = 1:numel(req)
        v = samples.(req{k});
        if ~isvector(v) || numel(v) ~= n
            error('tornadoPlot:SizeMismatch','Field "%s" must be a vector of length %d.', req{k}, n);
        end
    end

    % ---- Define labels (LaTeX) & construct data matrix ----
    factors = {'$Y_{\mathrm{EV}}$', ...
               '$\Delta\alpha$', ...
               '$\delta_{\mathrm{loss}}$', ...
               '$\nu$', ...
               '$L$', ...
               '$r$', ...
               '$E^{\mathrm{nom}}_{\mathrm{pack}}$'};

    alpha_diff = samples.alpha_CBP(:) - samples.alpha_RBP(:);

    X = [ samples.Y_EV(:), ...
          alpha_diff, ...
          samples.delta_loss(:), ...
          samples.nu(:), ...
          samples.L(:), ...
          samples.r_CBP(:), ...
          samples.E_pack_nom(:) ];


    % ---- Spearman correlation (NaN-safe) ----
    % corr returns a vector when X is (n×k) and y is (n×1)
    corrs = corr(X, y, 'Type','Spearman', 'Rows','pairwise');


    % ---- Sort by absolute value (descending) ----
    [sorted_corrs, idx] = sort(abs(corrs), 'descend'); 
    corrs_sorted = corrs(idx);
    factors_sorted = factors(idx);


    % ---- Plot ----
    axes(ax); 
    cla(ax); hold(ax, 'on');

    x_pos = 1:numel(factors_sorted);
    b = bar(ax, x_pos, corrs_sorted, 0.6, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor','none'); %#ok<NASGU>

    % Zero line for reference
    plot(ax, [0.5, numel(factors_sorted)+0.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75);

    % Axes labels/formatting
    ylabel(ax, '$\rho$ [-]', 'Interpreter','latex', 'FontName','Latin Modern Roman');
    xticks(ax, x_pos);
    xticklabels(ax, factors_sorted);
    yticks(ax, -1.0:0.5:1.0);
    ylim(ax, [-1.0 1.0]);
    xlim(ax, [0.5, numel(factors_sorted)+0.5]);

    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12, ...
            'FontName', 'Latin Modern Roman', ...
            'TickLabelInterpreter', 'latex', ...
            'LabelFontSizeMultiplier', 1);

    % Value labels on bars
    for i = 1:numel(corrs_sorted)
        val = corrs_sorted(i);
        if val >= 0
            y_offset = 0.02;  va = 'bottom';
        else
            y_offset = -0.02; va = 'top';
        end
        text(ax, x_pos(i), val + y_offset, sprintf('%.2f', val), ...
            'HorizontalAlignment','center', 'VerticalAlignment',va, ...
            'FontSize', 9, 'Interpreter','latex');
    end

    hold(ax, 'off');

end



% -------------------------------------------------------------------------
% PLOTSENSITIVITY1D_NOPERC_SUBPLOT
% -------------------------------------------------------------------------
function [x_zero, R_squared, RMSE] = plotSensitivity1D_NoPerc_Subplot(ax, x, y, xLabel, xlims, xticks, ylims, yticks, idxBest, idxBaseline, idxWorst, dotSize, legendLocation)
%PLOTSENSITIVITY1D_NOPERC_SUBPLOT 1D sensitivity scatter with best-fit line (kEUR).
%   [X_ZERO, R2, RMSE] = PLOTSENSITIVITY1D_NOPERC_SUBPLOT(AX, X, Y, XLABEL, XLIMS,
%   XTICKS, YLIMS, YTICKS, IDXBEST, IDXBASELINE, IDXWORST, DOTSIZE, LEGENDLOCATION)
%   plots Y vs X with:
%     • scatter of all samples (size=DOTSIZE),
%     • OLS best-fit line (red),
%     • break-even marker (Y=0) from the linear fit,
%     • highlighted points: Best (green *), Baseline (cyan *), Worst (magenta *).
%
%   INPUTS
%     ax        axes handle (subplot/tiledlayout target)
%     x, y      numeric vectors (same length), Y is ΔNPC in kEUR
%     xLabel    LaTeX string for x-axis label
%     xlims     [xmin xmax] (empty => auto from data)
%     xticks    vector of x tick positions
%     ylims     [ymin ymax] (empty => auto from data)
%     yticks    vector of y tick positions
%     idxBest   index of best sample
%     idxBaseline index of baseline sample
%     idxWorst  index of worst sample
%     dotSize   scalar marker size for scatter
%     legendLocation (optional) legend location string (e.g., "northwest")
%
%   OUTPUTS
%     x_zero    x at which fitted line crosses Y=0  (NaN if slope ~ 0)
%     R_squared coefficient of determination of the linear fit
%     RMSE      root-mean-squared error of the linear fit
%
%   NOTES
%   - NaN values in X or Y are removed pairwise.
%   - If X is (near) constant or variance is ~0, the fit is skipped and
%     X_ZERO, R_SQUARED, RMSE are set to NaN.

    % ---- Validation & cleaning ----
    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('plot1D:InvalidAxes','First input must be a valid axes handle.');
    end
    x = x(:); y = y(:);
    if numel(x) ~= numel(y)
        error('plot1D:SizeMismatch','x and y must have the same number of elements.');
    end
    mask = isfinite(x) & isfinite(y);
    x = x(mask); y = y(mask);
    if numel(x) < 3
        error('plot1D:TooFewPoints','Need at least 3 finite points for regression.');
    end

    % Colors
    MATLAB_blue = [0 0.4470 0.7410];
    MATLAB_red  = [0.8500 0.3250 0.0980];


    % ---- Axes & scatter ----
    axes(ax); 
    hold(ax, 'on');
    scatter(ax, x, y, dotSize, MATLAB_blue, 'filled');


    % ---- Linear regression (guard against near-constant X) ----
    if range(x) < max(1e-9, 1e-6*max(abs(x)))  % near-constant x

        %p = [NaN NaN];
        x_fit = [min(x) max(x)];
        y_fit = [mean(y) mean(y)];
        %y_pred = y*NaN;
        x_zero = NaN;
        R_squared = NaN;
        RMSE = NaN;
        plot(ax, x_fit, y_fit, '-', 'LineWidth', 2, 'Color', MATLAB_red); % flat line

    else

        p = polyfit(x, y, 1);                 % y = p1*x + p2
        x_fit = linspace(min(x), max(x), 200);
        y_fit = polyval(p, x_fit);
        plot(ax, x_fit, y_fit, '-', 'LineWidth', 2, 'Color', MATLAB_red);

        % Metrics
        y_pred = polyval(p, x);
        SS_res = sum((y - y_pred).^2);
        SS_tot = sum((y - mean(y)).^2);
        if SS_tot > 0
            R_squared = 1 - SS_res/SS_tot;
        else
            R_squared = NaN;
        end
        RMSE = sqrt(mean((y - y_pred).^2));

        % Break-even x at y=0 (if slope not ~0)
        if abs(p(1)) > eps
            x_zero = -p(2)/p(1);
        else
            x_zero = NaN;
        end

        % Break-even marker if inside plotting range
        if isfinite(x_zero) && (isempty(xlims) || (x_zero >= min(xlims) && x_zero <= max(xlims)))
            plot(ax, x_zero, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
        end

    end


    % ---- Highlight samples (ensure indices are valid after NaN filtering) ----
    % Map original indices to filtered mask.
    idxAll = find(mask);
    ib = idxAll(idxBest);    % original -> filtered
    ibl = idxAll(idxBaseline);
    iw = idxAll(idxWorst);

    if ib >= 1 && ib <= numel(x),   plot(ax, x(ib),  y(ib),  'g*', 'MarkerSize',10, 'LineWidth',2); end
    if ibl >= 1 && ibl <= numel(x), plot(ax, x(ibl), y(ibl), 'c*', 'MarkerSize',10, 'LineWidth',2); end
    if iw >= 1 && iw <= numel(x),   plot(ax, x(iw),  y(iw),  'm*', 'MarkerSize',10, 'LineWidth',2); end


    % ---- Labels & limits ----
    xlabel(ax, xLabel, 'Interpreter','latex', 'FontName','Latin Modern Roman');
    ylabel(ax, '$\Delta \mathrm{NPC}\ [\mathrm{kEUR}]$', 'Interpreter','latex', 'FontName','Latin Modern Roman');

    if isempty(xlims), xlims = [min(x) max(x)]; end
    xlim(ax, xlims);
    if ~isempty(xticks), set(ax, 'XTick', xticks); end

    if isempty(ylims), ylims = [min(y) max(y)]; end
    ylim(ax, ylims);
    if ~isempty(yticks), set(ax, 'YTick', yticks); end

    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12, 'FontName', 'Latin Modern Roman', ...
        'TickLabelInterpreter', 'latex', 'LabelFontSizeMultiplier', 1);


    % ---- Legend ----
    if ~exist('legendLocation','var') || isempty(legendLocation)
        legendLocation = "best";
    end

    % Assemble handles in the same order as labels
    hSamples   = findobj(ax, 'Type','scatter', '-and', 'Marker','o', '-depth',1);
    hFit       = findobj(ax, 'Type','line', '-and', 'Color', MATLAB_red, '-depth',1);
    hBreakEven = findobj(ax, 'Type','line', '-and', 'Color', [1 0 0], '-and', 'Marker','x', '-depth',1);
    hBest      = findobj(ax, 'Type','line', '-and', 'Color', [0 1 0], '-and', 'Marker','*', '-depth',1);
    hBase      = findobj(ax, 'Type','line', '-and', 'Color', [0 1 1], '-and', 'Marker','*', '-depth',1);
    hWorst     = findobj(ax, 'Type','line', '-and', 'Color', [1 0 1], '-and', 'Marker','*', '-depth',1);

    handles = [hSamples(1), hFit(1)];
    labels  = {'Samples','Best-fit line'};
    if ~isempty(hBreakEven), handles(end+1) = hBreakEven(1); labels{end+1} = 'Break-even'; end 
    if ~isempty(hBest),     handles(end+1) = hBest(1);      labels{end+1} = 'Best';       end 
    if ~isempty(hBase),     handles(end+1) = hBase(1);      labels{end+1} = 'Baseline';   end 
    if ~isempty(hWorst),    handles(end+1) = hWorst(1);     labels{end+1} = 'Worst';      end 

    lgd = legend(ax, handles, labels, 'Location', legendLocation, ...
        'FontSize', 9, 'Interpreter','latex', 'NumColumns', 2);
    lgd.ItemTokenSize = [10, 4];

    hold(ax, 'off');

end



% -------------------------------------------------------------------------
% PLOT1DWITHREGRESSION_SUBPLOT
% -------------------------------------------------------------------------
function [x_zero, R_squared, RMSE] = plot1DWithRegression_Subplot(ax, x, y, xLabel, xlims, xticks, ylims, yticks, idxBest, idxBaseline, idxWorst, dotSize)
%PLOT1DWITHREGRESSION_SUBPLOT 1D sensitivity plot with OLS fit (kEUR).
%   [X_ZERO, R2, RMSE] = PLOT1DWITHREGRESSION_SUBPLOT(AX, X, Y, XLABEL, XLIMS,
%   XTICKS, YLIMS, YTICKS, IDXBEST, IDXBASELINE, IDXWORST, DOTSIZE) plots Y vs X
%   using:
%       • scatter of all samples (size=DOTSIZE),
%       • OLS best-fit line (red),
%       • break-even marker at Y=0 from the fitted line,
%       • highlighted points: Best (green *), Baseline (cyan *), Worst (magenta *).
%
%   INPUTS
%     ax           axes handle (subplot/tiledlayout target)
%     x, y         numeric vectors (same length). Y is ΔNPC in kEUR.
%     xLabel       LaTeX string for x-axis label
%     xlims        [xmin xmax] for x (already in *percent* units if provided)
%     xticks       tick positions for x (percent units)
%     ylims        [ymin ymax] in kEUR
%     yticks       tick positions for y (kEUR)
%     idxBest      index of best sample
%     idxBaseline  index of baseline sample
%     idxWorst     index of worst sample
%     dotSize      scatter size (points^2)
%
%   OUTPUTS
%     x_zero       x (in *fraction* units) where fitted line crosses Y=0 (NaN if slope≈0)
%     R_squared    coefficient of determination of the linear fit
%     RMSE         root-mean-squared error of the fit
%
%   NOTES
%   - The x-axis is displayed in percent (×100), but X_ZERO is returned in fractional units.
%   - NaNs in X or Y are removed pairwise before fitting.

    % ---- Validate & clean ----
    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('plot1DReg:InvalidAxes','First input must be a valid axes handle.');
    end
    x = x(:); y = y(:);
    if numel(x) ~= numel(y)
        error('plot1DReg:SizeMismatch','x and y must have the same number of elements.');
    end
    mask = isfinite(x) & isfinite(y);
    x = x(mask); y = y(mask);
    if numel(x) < 3
        error('plot1DReg:TooFewPoints','Need at least 3 finite points for regression.');
    end

    % Colors
    MATLAB_blue = [0 0.4470 0.7410];
    MATLAB_red  = [0.8500 0.3250 0.0980];


    % ---- Fit (guard against near-constant X) ----
    nearConstX = range(x) < max(1e-9, 1e-6*max(abs(x)));
    axes(ax); 
    cla(ax); hold(ax,'on');

    scatter(ax, 100*x, y, dotSize, MATLAB_blue, 'filled'); % plot in percent

    if nearConstX

        % Flat fit: y = mean(y)
        x_fit = linspace(min(x), max(x), 2);
        y_fit = [mean(y) mean(y)];
        plot(ax, 100*x_fit, y_fit, '-', 'LineWidth', 2, 'Color', MATLAB_red);
        x_zero   = NaN;
        R_squared = NaN;
        RMSE      = sqrt(mean((y - mean(y)).^2));

    else

        p    = polyfit(x, y, 1);               % y = p1*x + p2
        x_fit = linspace(min(x), max(x), 200);
        y_fit = polyval(p, x_fit);
        plot(ax, 100*x_fit, y_fit, '-', 'LineWidth', 2, 'Color', MATLAB_red);

        y_pred  = polyval(p, x);
        SS_res  = sum((y - y_pred).^2);
        SS_tot  = sum((y - mean(y)).^2);
        R_squared = 1 - SS_res / SS_tot;
        RMSE      = sqrt(mean((y - y_pred).^2));

        if abs(p(1)) > eps

            x_zero = -p(2)/p(1);   % in fractional units

            % mark break-even only if it sits within shown x-limits (if provided)
            if isempty(xlims) || (100*x_zero >= min(xlims) && 100*x_zero <= max(xlims))
                plot(ax, 100*x_zero, 0, 'rx', 'MarkerSize',10, 'LineWidth',2);
            end

        else

            x_zero = NaN;

        end

    end

    % ---- Highlights (map original indices to filtered) ----
    idxAll = find(mask);
    iBest = idxAll(idxBest);
    iBase = idxAll(idxBaseline);
    iWorst = idxAll(idxWorst);

    if iBest >= 1 && iBest <= numel(x),   plot(ax, 100*x(iBest),  y(iBest),  'g*', 'MarkerSize',10, 'LineWidth',2); end
    if iBase >= 1 && iBase <= numel(x),   plot(ax, 100*x(iBase),  y(iBase),  'c*', 'MarkerSize',10, 'LineWidth',2); end
    if iWorst >= 1 && iWorst <= numel(x), plot(ax, 100*x(iWorst), y(iWorst), 'm*', 'MarkerSize',10, 'LineWidth',2); end


    % ---- Labels & axes ----
    xlabel(ax, xLabel, 'Interpreter','latex', 'FontName','Latin Modern Roman');
    ylabel(ax, '$\Delta \mathrm{NPC}\ [\mathrm{kEUR}]$', 'Interpreter','latex', 'FontName','Latin Modern Roman');

    if isempty(xlims), xlims = [100*min(x) 100*max(x)]; end
    xlim(ax, xlims);
    if ~isempty(xticks), set(ax, 'XTick', xticks); end

    if isempty(ylims), ylims = [min(y) max(y)]; end
    ylim(ax, ylims);
    if ~isempty(yticks), set(ax, 'YTick', yticks); end

    grid(ax,'on'); box(ax,'on');
    set(ax, 'FontSize', 12, 'FontName', 'Latin Modern Roman', ...
        'TickLabelInterpreter','latex', 'LabelFontSizeMultiplier',1);

    % Legend
    lgd = legend(ax, {'Samples','Best-fit line','Break-even','Best','Baseline','Worst'}, ...
        'Location','best', 'FontSize',9, 'Interpreter','latex', 'NumColumns',2);
    lgd.ItemTokenSize = [10, 4];

    hold(ax,'off');

end



% -------------------------------------------------------------------------
% PLOTCHEMISTRYBOXSCATTERWITHRECTANGLE_MODIFIED_SUBPLOT_V2
% -------------------------------------------------------------------------
function [nu_best, L_best] = ...
    plotChemistryBoxScatterWithRectangle_Modified_Subplot_v2 ...
        (axChemPlot, axRectPlot, samples, delta_NPC, ...
         rect_bounds, ylims_Main, yticks_Main, feasible, successThreshold, dotSize, ...
         xLim_axRectPlot, xTicks_axRectPlot, yLim_axRectPlot, yTicks_axRectPlot)
%PLOTCHEMISTRYBOXSCATTERWITHRECTANGLE_MODIFIED_SUBPLOT_V2
% Identify and visualize a dominant rectangle in (nu, L) with high ΔNPC success.
%
%   [NU_BEST, L_BEST] = PLOTCHEMISTRYBOXSCATTERWITHRECTANGLE_MODIFIED_SUBPLOT_V2( ... )
%   searches over axis-aligned rectangles anchored at (nu_min, L_min) up to
%   (nu_max, L_max) within provided bounds, and finds the largest (by area in
%   normalized space) whose share of successful outcomes (ΔNPC>0) meets or
%   exceeds SUCCESS THRESHOLD. It then:
%     • plots success/failure scatter in (nu,L) with the dominant rectangle;
%     • plots chemistry-wise box+scatter of ΔNPC (kEUR) inside/outside the rectangle.
%
%   INPUTS
%     axChemPlot      axes  – target axes for chemistry box+scatter panel
%     axRectPlot      axes  – target axes for (nu,L) rectangle/scatter panel
%     samples         struct with fields: nu, L, Chemistry (string array: "LFP"/"NMC")
%     delta_NPC       (n×1) double, ΔNPC in EUR
%     rect_bounds     2×2 double: [nu_min nu_max; L_min L_max]
%     ylims_Main      1×2 double (kEUR) – y-limits for chemistry panel
%     yticks_Main     vector    – y-ticks for chemistry panel
%     feasible        (n×1) logical – feasibility mask for (nu,L) search
%     successThreshold scalar in (0,1) – required share of ΔNPC>0
%     dotSize         scalar – scatter size
%     xLim_axRectPlot, xTicks_axRectPlot, yLim_axRectPlot, yTicks_axRectPlot
%                      axis formatting for rectangle panel (x in %, y in km)
%
%   OUTPUTS
%     nu_best, L_best  scalars – dominant rectangle upper-right corner;
%                              NaN if none meets the threshold.
%
%   NOTES
%   • Search is on a 101×101 grid in normalized coordinates for unbiased area.
%   • ΔNPC is converted to kEUR for the chemistry panel; rectangle panel uses EUR sign only.

    % ---------- Basic validation ----------
    if ~ishandle(axChemPlot) || ~strcmp(get(axChemPlot,'Type'),'axes')
        error('chemRect:InvalidAxes','axChemPlot must be an axes handle.');
    end
    if ~ishandle(axRectPlot) || ~strcmp(get(axRectPlot,'Type'),'axes')
        error('chemRect:InvalidAxes','axRectPlot must be an axes handle.');
    end
    req = {'nu','L','Chemistry'};
    for k = 1:numel(req)
        if ~isfield(samples, req{k})
            error('chemRect:MissingField','samples.%s is required.', req{k});
        end
    end
    if size(rect_bounds,1)~=2 || size(rect_bounds,2)~=2
        error('chemRect:BoundsShape','rect_bounds must be 2x2: [nu_min nu_max; L_min L_max].');
    end
    if ~(isscalar(successThreshold) && successThreshold>0 && successThreshold<1)
        error('chemRect:BadThreshold','successThreshold must be in (0,1).');
    end


    % ---------- Prep data ----------
    delta_NPC_k = delta_NPC(:) / 1000;      % kEUR for chemistry panel
    nu = samples.nu(:);
    L  = samples.L(:);
    feasible = logical(feasible(:));
    n = numel(nu);
    if ~isvector(delta_NPC_k) || numel(delta_NPC_k)~=n
        error('chemRect:SizeMismatch','delta_NPC length must match samples.nu/L.');
    end

    is_LFP = (samples.Chemistry(:) == "LFP");
    is_NMC = (samples.Chemistry(:) == "NMC");

    % Feasible region within provided bounds
    nu_bounds = rect_bounds(1,:);
    L_bounds  = rect_bounds(2,:);
    nu_min = nu_bounds(1); nu_max = nu_bounds(2); nu_range = nu_max - nu_min;
    L_min  = L_bounds(1);  L_max  = L_bounds(2);  L_range  = L_max  - L_min;

    idx_mu_safe = feasible & ...
                  (nu >= nu_min) & (nu <= nu_max) & ...
                  (L  >= L_min ) & (L  <= L_max );

    nu_safe   = nu(idx_mu_safe);
    L_safe    = L(idx_mu_safe);
    dNPC_safe = delta_NPC(idx_mu_safe);     % EUR for success computation


    % ---------- Normalized grid search over rectangles ----------
    nu_grid_norm = linspace(0, 1, 101);
    L_grid_norm  = linspace(0, 1, 101);

    R = numel(nu_grid_norm);
    A = numel(L_grid_norm);
    AreaMat        = -inf(R, A);   % initialize to -Inf so invalid stay out
    SuccessRateMat = NaN(R, A);

    if ~isempty(nu_safe)
        % Precompute normalized safe coords
        nu_norm_safe = (nu_safe - nu_min) / max(eps, nu_range);
        L_norm_safe  = (L_safe  - L_min ) / max(eps, L_range);

        for i = 1:R
            nu_max_norm = nu_grid_norm(i);
            for j = 1:A
                L_max_norm = L_grid_norm(j);

                in_idx = (nu_norm_safe <= nu_max_norm) & (L_norm_safe <= L_max_norm);
                count  = sum(in_idx);
                if count == 0
                    continue;
                end
                success_rate = mean(dNPC_safe(in_idx) > 0);
                SuccessRateMat(i,j) = success_rate;

                if success_rate >= successThreshold
                    AreaMat(i,j) = nu_max_norm * L_max_norm; % unbiased area
                end
            end
        end
    end


    % ---------- Pick dominant rectangle ----------
    [maxArea, linIdx] = max(AreaMat(:)); 
    if isfinite(maxArea)
        [i_best, j_best] = ind2sub(size(AreaMat), linIdx);
        nu_best = nu_min + nu_grid_norm(i_best) * nu_range;
        L_best  = L_min  + L_grid_norm(j_best)  * L_range;
        % Actual success inside dominant rect (recompute on safe set)
        in_best = (nu_safe <= nu_best) & (L_safe <= L_best);
        success_rate_best = mean(dNPC_safe(in_best) > 0);
    else
        nu_best = NaN; L_best = NaN;
        success_rate_best = NaN;
    end


    % ---------- Rectangle panel: scatter + rectangle overlay ----------
    axes(axRectPlot); cla(axRectPlot); hold(axRectPlot,'on');

    if ~isempty(nu_safe)
        success_mask = dNPC_safe > 0;
        hSucc = scatter(axRectPlot, nu_safe(success_mask)*100, L_safe(success_mask), ...
                        dotSize, [0 0.6 0.8], 'filled', 'DisplayName', '$\Delta \mathrm{NPC} > 0$');
        hFail = scatter(axRectPlot, nu_safe(~success_mask)*100, L_safe(~success_mask), ...
                        dotSize, [0.8 0 0], 'filled', 'DisplayName', '$\Delta \mathrm{NPC} < 0$'); %#ok<NASGU>

        % Rectangle (only if valid)
        if isfinite(nu_best) && isfinite(L_best)
            x_rect = [nu_min, nu_best, nu_best, nu_min, nu_min] * 100;
            y_rect = [L_min,  L_min,  L_best,  L_best,  L_min];
            rectLabel = sprintf('$%.1f\\,\\%%\\ \\mathrm{of\\ safe\\ samples}$', 100*success_rate_best);
            hpatch = patch(axRectPlot, 'XData', x_rect, 'YData', y_rect, ...
                'FaceColor','none','EdgeColor','k','LineWidth',2, ...
                'LineStyle','-');%,'DisplayName', rectLabel);
            hpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    else
        % No safe points in bounds
        text(axRectPlot, 0.5, 0.5, '\emph{No feasible points within bounds}', ...
             'Units','normalized','HorizontalAlignment','center', ...
             'Interpreter','latex','FontSize',10);
    end

    xlabel(axRectPlot, '$\nu$ [\%]', 'Interpreter','latex');
    ylabel(axRectPlot, '$L$ [km]', 'Interpreter','latex');
    xlim(axRectPlot, xLim_axRectPlot);
    xticks(axRectPlot, xTicks_axRectPlot);
    ylim(axRectPlot, yLim_axRectPlot);
    yticks(axRectPlot, yTicks_axRectPlot);
    grid(axRectPlot,'on'); box(axRectPlot,'on');
    set(axRectPlot, 'FontSize',12, 'FontName','Latin Modern Roman', 'TickLabelInterpreter','latex');

    % Legend (succ/fail; rectangle label appears in legend if drawn)
    if exist('hSucc','var')
        lgd2 = legend(axRectPlot, 'Interpreter','latex', 'Location','north', ...
                      'Orientation','horizontal','FontSize',9,'NumColumns',2);
        lgd2.ItemTokenSize = [10, 4];
    end

    hold(axRectPlot,'off');


    % ---------- Chemistry panel: box+scatter ----------
    axes(axChemPlot); cla(axChemPlot); hold(axChemPlot,'on');
    rng(42); jitter = 0.05;

    % Feasible-by-chemistry ΔNPC (kEUR) + inside-dominant-rect (if any)
    dNPC_LFP_feas = delta_NPC_k(is_LFP & feasible);
    dNPC_NMC_feas = delta_NPC_k(is_NMC & feasible);

    if isfinite(nu_best) && isfinite(L_best)
        in_dom_rect = feasible & (nu <= nu_best) & (L <= L_best);
        dNPC_LFP_dom = delta_NPC_k(is_LFP & in_dom_rect);
        dNPC_NMC_dom = delta_NPC_k(is_NMC & in_dom_rect);
    else
        dNPC_LFP_dom = []; dNPC_NMC_dom = [];
    end

    % Helper (provided externally)
    scatterPoints_subplot(axChemPlot, 0.6, dNPC_LFP_feas, [0 0.4470 0.7410], jitter);
    scatterPoints_subplot(axChemPlot, 1.4, dNPC_LFP_dom,  [0.4660 0.6740 0.1880], jitter);
    h1 = boxchart(axChemPlot, 0.6*ones(size(dNPC_LFP_feas)), dNPC_LFP_feas, ...
                  'BoxFaceAlpha',0.5,'BoxFaceColor',[0 0.4470 0.7410],'HandleVisibility','off');
    h2 = boxchart(axChemPlot, 1.4*ones(size(dNPC_LFP_dom)),  dNPC_LFP_dom,  ...
                  'BoxFaceAlpha',0.5,'BoxFaceColor',[0.4660 0.6740 0.1880],'HandleVisibility','off');
    h1.JitterOutliers='on'; h1.MarkerStyle='none';
    h2.JitterOutliers='on'; h2.MarkerStyle='none';

    scatterPoints_subplot(axChemPlot, 2.1, dNPC_NMC_feas, [0.8500 0.3250 0.0980], jitter);
    scatterPoints_subplot(axChemPlot, 2.9, dNPC_NMC_dom,  [0.9290 0.6940 0.1250], jitter);
    h3 = boxchart(axChemPlot, 2.1*ones(size(dNPC_NMC_feas)), dNPC_NMC_feas, ...
                  'BoxFaceAlpha',0.5,'BoxFaceColor',[0.8500 0.3250 0.0980],'HandleVisibility','off');
    h4 = boxchart(axChemPlot, 2.9*ones(size(dNPC_NMC_dom)),  dNPC_NMC_dom,  ...
                  'BoxFaceAlpha',0.5,'BoxFaceColor',[0.9290 0.6940 0.1250],'HandleVisibility','off');
    h3.JitterOutliers='on'; h3.MarkerStyle='none';
    h4.JitterOutliers='on'; h4.MarkerStyle='none';

    xlabel(axChemPlot, 'Chemistry', 'Interpreter','latex');
    ylabel(axChemPlot, '$\Delta \mathrm{NPC}\ [\mathrm{kEUR}]$', 'Interpreter','latex');
    xlim(axChemPlot, [0 3.5]); ylim(axChemPlot, ylims_Main);
    set(axChemPlot, 'YTick', yticks_Main);
    grid(axChemPlot,'on'); box(axChemPlot,'on');
    lgd1 = legend(axChemPlot, {'LFP','LFP (\textbf{D})','NMC','NMC (\textbf{D})'}, ...
                  'FontSize',9,'Interpreter','latex','Location','north','NumColumns',2);
    lgd1.ItemTokenSize = [10, 4];
    set(axChemPlot, 'FontSize',12, 'FontName','Latin Modern Roman', 'TickLabelInterpreter','latex');

    % Median annotations for dominant groups (if exist)
    if ~isempty(dNPC_LFP_dom)

        med_LFP_dom = median(dNPC_LFP_dom);
        label_LFP_D = sprintf('Median:\n%.1f EUR', med_LFP_dom*1000);
        edgeColor = [0.4660 0.6740 0.1880];
        bgColor   = edgeColor + (1 - edgeColor)*0.5;
        text(axChemPlot, 1.4-0.02, ylims_Main(1) + 0.18*range(ylims_Main), label_LFP_D, ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontSize',9,'Color','black','Interpreter','latex', ...
            'EdgeColor',edgeColor,'BackgroundColor',bgColor,'Margin',1);

    end

    if ~isempty(dNPC_NMC_dom)

        med_NMC_dom = median(dNPC_NMC_dom);
        label_NMC_D = sprintf('Median:\n%.1f EUR', med_NMC_dom*1000);
        edgeColor = [0.9290 0.6940 0.1250];
        bgColor   = edgeColor + (1 - edgeColor)*0.5;
        text(axChemPlot, 2.9-0.02, ylims_Main(1) + 0.18*range(ylims_Main), label_NMC_D, ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontSize',9,'Color','black','Interpreter','latex', ...
            'EdgeColor',edgeColor,'BackgroundColor',bgColor,'Margin',1);

    end

    % Tidy x ticks
    set(axChemPlot,'XTick',[]);

    hold(axChemPlot,'off');

end



% -------------------------------------------------------------------------
% SCATTERPOINTS_SUBPLOT
% -------------------------------------------------------------------------
function scatterPoints_subplot(ax, xpos, ydata, color, jitter)
%SCATTERPOINTS_SUBPLOT Jittered scatter plot for box+scatter subplots
%
%   SCATTERPOINTS_SUBPLOT(AX, XPOS, YDATA, COLOR, JITTER) creates a
%   horizontally jittered scatter plot of data points at position XPOS
%   on the x-axis. Useful for overlaying on boxplots in tiled layouts.
%
%   INPUTS:
%     ax     - Axes handle to plot into (e.g. from subplot or nexttile)
%     xpos   - Scalar x-position of the box/cluster center (e.g. 0.6, 1.4)
%     ydata  - Vector of y-values (one per sample)
%     color  - 1×3 RGB vector (e.g. [0 0.447 0.741])
%     jitter - Scalar specifying horizontal jitter width (e.g. 0.05)
%
%   EXAMPLE:
%     scatterPoints_subplot(gca, 1.4, randn(100,1), [0 0.4470 0.7410], 0.05);
%
%   NOTES:
%   - Markers are semi-transparent ('MarkerFaceAlpha' = 0.7).
%   - X-values are randomly perturbed in the range xpos ± jitter/2.


    % ---- Validate inputs ----
    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('scatterPoints_subplot:InvalidAxes', ...
              'First argument must be a valid axes handle.');
    end
    if ~isnumeric(xpos) || ~isscalar(xpos)
        error('scatterPoints_subplot:InvalidXpos', ...
              'xpos must be a numeric scalar.');
    end
    if ~isnumeric(ydata) || ~isvector(ydata)
        error('scatterPoints_subplot:InvalidYdata', ...
              'ydata must be a numeric vector.');
    end
    if ~isnumeric(color) || numel(color)~=3
        error('scatterPoints_subplot:InvalidColor', ...
              'color must be a 1x3 RGB numeric vector.');
    end
    if ~isnumeric(jitter) || ~isscalar(jitter) || jitter < 0
        error('scatterPoints_subplot:InvalidJitter', ...
              'jitter must be a nonnegative scalar.');
    end


    % ---- Plot data ----
    % Ensure we're plotting in the right axes
    axes(ax);  
    xj = xpos + (rand(size(ydata)) - 0.5) * jitter;  % random horizontal jitter
    scatter(ax, xj, ydata, 25, color, 'filled', 'MarkerFaceAlpha', 0.7);

end



% -------------------------------------------------------------------------
% POSITIONFIGUREONMONITOR
% -------------------------------------------------------------------------
function positionFigureOnMonitor(fig, monitorIdx, paperW_cm, paperH_cm)
%POSITIONFIGUREONMONITOR Place figure top-left of monitor with paper sizing.
    screens = get(0,'MonitorPositions');
    dpi     = get(0,'ScreenPixelsPerInch');
    cm2in   = 1/2.54;

    width_px  = round(paperW_cm * cm2in * dpi);
    height_px = round(paperH_cm * cm2in * dpi);

    if size(screens,1) >= monitorIdx
        mon = screens(monitorIdx,:);
    else
        mon = screens(1,:);
    end

    fig_left   = mon(1);
    fig_bottom = mon(2) + mon(4) - height_px - 50; % 50 px top margin
    set(fig, 'Units','pixels', 'Position',[fig_left, fig_bottom, width_px, height_px]);

    set(fig, 'PaperUnits','centimeters', ...
        'PaperSize',[paperW_cm paperH_cm], ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[0 0 paperW_cm paperH_cm]);
end



% -------------------------------------------------------------------------
% MUSTHAVE
% -------------------------------------------------------------------------
function mustHave(s, names)
%MUSTHAVE Assert that struct S has all fields in cellstr NAMES.
    for k = 1:numel(names)
        if ~isfield(s, names{k})
            error('computeNPC:MissingField','Required field "%s" is missing.', names{k});
        end
    end
end


function plotSensitivity1D_Subplot(ax, x, y, xLabel, xlims, xticks, ylims, yticks, idxBest, idxBaseline, idxWorst, dotSize)
% plotSensitivity1D_Subplot - 1D sensitivity scatter plot for tiledlayout/subplot use.
%
% Inputs:
%   ax      - Axes handle (subplot or tile)
%   x, y    - Data vectors
%   xLabel  - LaTeX string for x-axis label
%   ylims   - Y-axis limits [min, max]
%   xlims   - X-axis limits [min, max] (optional)

    MATLAB_blue = [0 0.4470 0.7410];
    MATLAB_red  = [0.8500 0.3250 0.0980];
    highlightColor = [0.9290, 0.6940, 0.1250]; % Orange







    
    
    
    
    
    
    
    
    
    axes(ax);
    hold(ax, 'on');

    scatter(ax, 100*x, y, dotSize, MATLAB_blue, 'filled');

    % Highlight specific samples
    plot(ax, 100*x(idxBest), y(idxBest), 'g*', 'MarkerSize',10, 'LineWidth',2);
    plot(ax, 100*x(idxBaseline), y(idxBaseline), 'c*', 'MarkerSize',10, 'LineWidth',2);
    plot(ax, 100*x(idxWorst), y(idxWorst), 'm*', 'MarkerSize',10, 'LineWidth',2);    

    % Axis labels
    xlabel(ax, xLabel, 'Interpreter','latex', 'FontName', 'Latin Modern Roman');
    ylabel(ax, '$\Delta \mathrm{NPC}\ [\mathrm{kEUR}]$', 'Interpreter','latex', 'FontName', 'Latin Modern Roman');

    % Axis limits
    if isempty(xlims)
        xlims = [min(x), max(x)];
    end
    xlim(ax, xlims);
    set(ax, 'XTick', xticks);  

    if isempty(ylims)
        ylims = [min(y), max(y)];
    end    
    ylim(ax, ylims);
    set(ax, 'YTick', yticks);       

    grid(ax, 'on');
    box(ax, 'on');

    set(ax, 'FontSize', 12, 'FontName', 'Latin Modern Roman', ...
        'TickLabelInterpreter', 'latex', ...
        'LabelFontSizeMultiplier', 1);
    lgd1 = legend(ax, {'Samples', 'Best', 'Baseline', 'Worst'}, ...
        'Location', 'best', 'FontSize', 9, 'Interpreter','latex', 'NumColumns', 2);
    lgd1.ItemTokenSize = [10, 4]; % [marker length, marker padding] in points        

    hold(ax, 'off');
end




function [x_best, y_best] = ...
    plotDominantRectangleIn2DFeatureSpace( ...
        axChemPlot, axRectPlot, samples, delta_NPC, ...
        rect_bounds, ylims_Main, yticks_Main, feasible, successThreshold, dotSize, ...
        xLim_axRectPlot, xTicks_axRectPlot, ...
        yLim_axRectPlot, yTicks_axRectPlot, ...
        xVarName, yVarName, xAxisInPerc, yAxisInPerc, xVarNameLabel, yVarNameLabel, ...
        altSuccessThreshold, plotAltCriterion)

% Extract x and y variables dynamically from samples
xVar = samples.(xVarName);
yVar = samples.(yVarName);

x_bounds = rect_bounds(1,:);
y_bounds = rect_bounds(2,:);

% Convert delta NPC to kEUR
delta_NPC_k = delta_NPC / 1000;

is_LFP = samples.Chemistry == "LFP";
is_NMC = samples.Chemistry == "NMC";

dNPC_LFP_feas = delta_NPC_k(is_LFP);
dNPC_NMC_feas = delta_NPC_k(is_NMC);

% Feasible indices within bounds
idx_safe = feasible & ...
            (xVar >= x_bounds(1)) & (xVar <= x_bounds(2)) & ...
            (yVar >= y_bounds(1)) & (yVar <= y_bounds(2));

x_safe = xVar(idx_safe);
y_safe = yVar(idx_safe);
dNPC_safe = delta_NPC_k(idx_safe);

% Normalize to [0, 1] for unbiased rectangle search
x_min = x_bounds(1);
x_range = x_bounds(2) - x_bounds(1);
y_min = y_bounds(1);
y_range = y_bounds(2) - y_bounds(1);

x_norm_safe = (x_safe - x_min) / x_range;
y_norm_safe = (y_safe - y_min) / y_range;

x_grid_norm = linspace(0, 1, 101);
y_grid_norm = linspace(0, 1, 101);

R = length(x_grid_norm);
A = length(y_grid_norm);

AreaMat = NaN(R, A);
SuccessRateMat = NaN(R, A);

for i = 1:R
    for j = 1:A
        x_max_norm = x_grid_norm(i);
        y_max_norm = y_grid_norm(j);

        x_max_real = x_min + x_max_norm * x_range;
        y_max_real = y_min + y_max_norm * y_range;

        in_idx = (x_safe <= x_max_real) & (y_safe <= y_max_real);
        if sum(in_idx) == 0
            continue;
        end

        success = dNPC_safe(in_idx) > 0;
        success_rate = sum(success) / sum(in_idx);

        SuccessRateMat(i,j) = success_rate;
        AreaMat(i,j) = x_max_norm * y_max_norm;
    end
end

valid = SuccessRateMat >= successThreshold;
AreaValid = AreaMat;
AreaValid(~valid) = -Inf;

if ~any(valid, 'all')
    fprintf('No rectangle satisfies the %.2f%% success rate constraint.\n', successThreshold*100);
end

[maxArea, idx] = max(AreaValid(:));
[i_best, j_best] = ind2sub(size(AreaValid), idx);

x_best = x_min + x_grid_norm(i_best) * x_range;
y_best = y_min + y_grid_norm(j_best) * y_range;



x_alt = NaN;
y_alt = NaN;

if plotAltCriterion
    validAlt = SuccessRateMat >= altSuccessThreshold;
    AreaValidAlt = AreaMat;
    AreaValidAlt(~validAlt) = -Inf;

    if any(validAlt, 'all')
        [~, idx_alt] = max(AreaValidAlt(:));
        [i_alt, j_alt] = ind2sub(size(AreaValidAlt), idx_alt);

        x_alt = x_min + x_grid_norm(i_alt) * x_range;
        y_alt = y_min + y_grid_norm(j_alt) * y_range;
    else
        fprintf('No rectangle satisfies the %.2f%% alt success rate constraint.\n', altSuccessThreshold*100);
    end
end



% Apply dominant rectangle
in_dominant_rectangle = feasible & ...
                        (xVar <= x_best) & (yVar <= y_best);

dNPC_LFP_dom = delta_NPC_k(in_dominant_rectangle & is_LFP);
dNPC_NMC_dom = delta_NPC_k(in_dominant_rectangle & is_NMC);

%% === Plot rectangle on feature space ===
axes(axRectPlot); hold(axRectPlot, 'on');

success = dNPC_safe > 0;
if (xAxisInPerc)
    xGood = 100*x_safe(success);
    xNotGood = 100*x_safe(~success);
else
    xGood = x_safe(success);
    xNotGood = x_safe(~success);
end

if (yAxisInPerc)
    yGood = 100*y_safe(success);
    yNotGood = 100*y_safe(~success);
else
    yGood = y_safe(success);
    yNotGood = y_safe(~success);
end
h1_ax2 = scatter(axRectPlot, xGood, yGood, dotSize, [0 0.6 0.8], 'filled', 'DisplayName', '$\Delta \mathrm{NPC} > 0$');
h2_ax2 = scatter(axRectPlot, xNotGood, yNotGood, dotSize, [0.8 0 0], 'filled', 'DisplayName', '$\Delta \mathrm{NPC} < 0\,$');

if (xAxisInPerc)
    x_rect = 100*[x_bounds(1), x_best, x_best, x_bounds(1), x_bounds(1)];
else
    x_rect = [x_bounds(1), x_best, x_best, x_bounds(1), x_bounds(1)];
end

if (yAxisInPerc)
    y_rect = 100*[y_min, y_min, y_best, y_best, y_min];
else
    y_rect = [y_min, y_min, y_best, y_best, y_min];
end
patch(axRectPlot, 'XData', x_rect, 'YData', y_rect, ...
    'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2, ...
    'LineStyle', '-', 'DisplayName', '99.7\% Y (\textbf{D})');


if plotAltCriterion && ~isnan(x_alt) && ~isnan(y_alt)
    if xAxisInPerc
        x_rect_alt = 100 * [x_bounds(1), x_alt, x_alt, x_bounds(1), x_bounds(1)];
    else
        x_rect_alt = [x_bounds(1), x_alt, x_alt, x_bounds(1), x_bounds(1)];
    end
    if yAxisInPerc
        y_rect_alt = 100 * [y_min, y_min, y_alt, y_alt, y_min];
    else
        y_rect_alt = [y_min, y_min, y_alt, y_alt, y_min];
    end
    patch(axRectPlot, 'XData', x_rect_alt, 'YData', y_rect_alt, ...
        'FaceColor', 'none', 'EdgeColor', [0.3 0.3 0.3], ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.1f%% Threshold', altSuccessThreshold*100));
end



xlabel(axRectPlot, xVarNameLabel, 'Interpreter', 'latex');
ylabel(axRectPlot, yVarNameLabel, 'Interpreter', 'latex');
    lgd2 = legend(axRectPlot, [h1_ax2, h2_ax2], 'Interpreter','latex', ...    
        'Location','north', 'Orientation','horizontal', ...
        'FontSize', 9, 'NumColumns', 2);
    %lgd2 = legend(axRectPlot, [h1_ax2, h2_ax2], 'Interpreter','latex', ...        
    lgd2.ItemTokenSize = [10, 2]; % [marker length, marker padding] in points
xlim(axRectPlot, xLim_axRectPlot);        
xticks(axRectPlot, xTicks_axRectPlot);  
ylim(axRectPlot, yLim_axRectPlot);      
yticks(axRectPlot, yTicks_axRectPlot);    
grid(axRectPlot, 'on'); box on;
set(axRectPlot, 'FontSize', 12, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter','latex');

hold(axRectPlot, 'off');

%% === Box + Scatter Plot for Chemistries ===
axes(axChemPlot); hold(axChemPlot, 'on'); rng(42); jitter = 0.05;

scatterPoints_subplot(axChemPlot, 0.6, dNPC_LFP_feas, [0 0.4470 0.7410], jitter);
scatterPoints_subplot(axChemPlot, 1.4, dNPC_LFP_dom, [0.4660, 0.6740, 0.1880], jitter);
boxchart(axChemPlot, 0.6*ones(size(dNPC_LFP_feas)), dNPC_LFP_feas, 'BoxFaceAlpha',0.5, 'BoxFaceColor',[0 0.4470 0.7410],'HandleVisibility', 'off');
boxchart(axChemPlot, 1.4*ones(size(dNPC_LFP_dom)), dNPC_LFP_dom, 'BoxFaceAlpha',0.5, 'BoxFaceColor',[0.4660, 0.6740, 0.1880],'HandleVisibility', 'off');

scatterPoints_subplot(axChemPlot, 2.1, dNPC_NMC_feas, [0.8500 0.3250 0.0980], jitter);
scatterPoints_subplot(axChemPlot, 2.9, dNPC_NMC_dom, [0.9290 0.6940 0.1250], jitter);
boxchart(axChemPlot, 2.1*ones(size(dNPC_NMC_feas)), dNPC_NMC_feas, 'BoxFaceAlpha',0.5, 'BoxFaceColor',[0.8500 0.3250 0.0980],'HandleVisibility', 'off');
boxchart(axChemPlot, 2.9*ones(size(dNPC_NMC_dom)), dNPC_NMC_dom, 'BoxFaceAlpha',0.5, 'BoxFaceColor',[0.9290 0.6940 0.1250],'HandleVisibility', 'off');

xlabel(axChemPlot, 'Chemistry', 'Interpreter','latex');
ylabel(axChemPlot, '$\Delta \mathrm{NPC}\ [\mathrm{kEUR}]$', 'Interpreter','latex');
xlim(axChemPlot, [0 3.5]); ylim(axChemPlot, ylims_Main);
set(axChemPlot, 'YTick', yticks_Main); grid(axChemPlot, 'on');
legend(axChemPlot, {'LFP','LFP (\textbf{D})','NMC','NMC (\textbf{D})'}, ...
       'FontSize',9,'Interpreter','latex','Location','north', 'NumColumns', 2);
set(axChemPlot, 'FontSize', 12, 'FontName', 'Latin Modern Roman', 'TickLabelInterpreter','latex');
box on;
set(axChemPlot, 'XTick', []);

% Add median labels
med_LFP_dom  = median(dNPC_LFP_dom);
med_NMC_dom  = median(dNPC_NMC_dom);
offset = 0.18 * range(ylims_Main);
label_LFP_D  = sprintf('Median:\n%.1f EUR', med_LFP_dom*1000);
label_NMC_D  = sprintf('Median:\n%.1f EUR', med_NMC_dom*1000);
bgColor1 = [0.4660, 0.6740, 0.1880] + (1 - [0.4660, 0.6740, 0.1880]) * 0.5;
bgColor2 = [0.9290, 0.6940, 0.1250] + (1 - [0.9290, 0.6940, 0.1250]) * 0.5;

text(axChemPlot, 1.4-0.02, ylims_Main(1) + offset, label_LFP_D, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'Interpreter', 'latex', 'EdgeColor', [0.4660, 0.6740, 0.1880], ...
    'BackgroundColor', bgColor1, 'Margin', 1);

text(axChemPlot, 2.9-0.02, ylims_Main(1) + offset, label_NMC_D, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'Interpreter', 'latex', 'EdgeColor', [0.9290 0.6940 0.1250], ...
    'BackgroundColor', bgColor2, 'Margin', 1);

end