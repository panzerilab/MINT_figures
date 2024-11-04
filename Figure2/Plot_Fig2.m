clc, clear all, close all;
% -------------------------------------------------------------------------
% Plot function for Figure 2
% -------------------------------------------------------------------------

plotRows = 5;
plotColumns = 4;
width_cm = 19.05;
height_cm = 22.23;
width_inch = width_cm / 2.54;
height_inch = (width_inch * plotRows) / plotColumns;
figure_handle = figure('Units', 'inches', 'Position', [1, 1, width_inch, height_inch]);
t = tiledlayout(5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% Setzen der Standardwerte für Text
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 8);  % Für Text
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 8);  % Für Achsentitel und -beschriftungen
set(0, 'DefaultLegendFontName', 'Arial', 'DefaultLegendFontSize', 8);  % Für Legenden
set(0, 'DefaultAxesFontSize', 8);
fontsize = 8;

colorJoint = [102 102 102]./255;
colorRed = [184 209 240]./255;
colorSyn = [134 150 167]./255;

colorLin_2 = [197 197 197]./255;

colorInfoBreak = [182 186 183]./255;

colorLin = [203 191 133]./255;
colorLinSh = [234 224 182]./255;
colorRBF = [124 75 74]./255;
colorRBFSh = [176 122 121]./255;

% -------------------------------------------------------------------------
% Step 1: Load data
% -------------------------------------------------------------------------
for col= 1:plotColumns
    row = 1;
    tileIndex = (row - 1) * plotColumns + col;
    Joint_all = [];
    ILIN_all = [];
    ISS_all = [];
    ICI_all = [];
    ICD_all = [];

    Lin_all = [];
    RBF_all = [];
    LinSh_all = [];
    RBFSh_all = [];

    Syn_all = [];
    Red_all = [];
    Unq1_all = [];
    Unq2_all = [];
    PIDJoint_all = [];
    if col == 1 || col == 2
        load('Results/Simulation_Figure2.mat');

        if col == 1
            scenario_labels = {'onlyNoise'};
        else
            scenario_labels = {'limitNoise'};
        end
        scLab = scenario_labels{:};
        for bdwIdx = 1:numel(info_bdw_terms)
            bdwLab = info_bdw_terms{bdwIdx};
            eval([bdwLab '_all = reshape(MI_breakdown.(scLab).(bdwLab), 1, []);']);
        end
        Lin_all = vertcat(MI_pop.(scLab).linear);
        RBF_all = vertcat(MI_pop.(scLab).RBF);
        Lin_all = Lin_all(:)';
        RBF_all = RBF_all(:)';

        LinSh_all = vertcat(MISh_pop.(scLab).linear);
        RBFSh_all = vertcat(MISh_pop.(scLab).RBF);
        LinSh_mean = (LinSh_all(:,1:5) + LinSh_all(:,6:10)) / 2;
        RBFSh_mean = (RBFSh_all(:,1:5) + RBFSh_all(:,6:10)) / 2;
        LinSh_mean = LinSh_mean(:)';
        RBFSh_mean = RBFSh_mean(:)';
        LinSh_all = LinSh_all(:)';
        RBFSh_all = RBFSh_all(:)';


        PIDJoint_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,1,:)), 1, []));
        Syn_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,2,:)), 1, []));
        Red_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,3,:)), 1, []));
        Unq1_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,4,:)), 1, []));
        Unq2_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,5,:)), 1, []));

    elseif col == 4
        load('Results/A1_Figure2.mat');
        for sessIdx = 1:numel(field_names)
            sessLab = field_names{sessIdx};
            Joint_all = [Joint_all, (MI_breakdown.(sessLab).Joint)];
            ILIN_all = [ILIN_all, (MI_breakdown.(sessLab).ILIN)];
            ISS_all = [ISS_all, (MI_breakdown.(sessLab).ISS)];
            ICI_all = [ICI_all, (MI_breakdown.(sessLab).ICI)];
            ICD_all = [ICD_all, (MI_breakdown.(sessLab).ICD)];

            Lin_all = [Lin_all, MI_pop.(sessLab).linear];
            RBF_all = [RBF_all, MI_pop.(sessLab).RBF];

            LinSh_all = [LinSh_all; MISh_pop.(sessLab).linear];
            RBFSh_all = [RBFSh_all; MISh_pop.(sessLab).RBF];

            % PID
            PIDJoint_all = [PIDJoint_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,1,:)), 1, []))];
            Syn_all = [Syn_all, cell2mat(reshape(squeeze(PID_result.(sessLab)(:,2,:)), 1, []))];
            Red_all = [Red_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,3,:)), 1, []))];
            Unq1_all = [Unq1_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,4,:)), 1, []))];
            Unq2_all = [Unq2_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,5,:)), 1, []))];
        end
        LinSh_mean = (LinSh_all(:,1:2) + LinSh_all(:,3:4)) / 2;
        RBFSh_mean = (RBFSh_all(:,1:2) + RBFSh_all(:,3:4)) / 2;
        LinSh_mean = LinSh_mean(:)';
        RBFSh_mean = RBFSh_mean(:)';
        LinSh_all = LinSh_all(:)';
        RBFSh_all = RBFSh_all(:)';
    elseif col == 3
        load('Results/CA1_Figure2.mat');
        for sessIdx = 1:numel(field_names)
            sessLab = field_names{sessIdx};
            Joint_all = [Joint_all, (MI_breakdown.(sessLab).Joint)];
            ILIN_all = [ILIN_all, (MI_breakdown.(sessLab).ILIN)];
            ISS_all = [ISS_all, (MI_breakdown.(sessLab).ISS)];
            ICI_all = [ICI_all, (MI_breakdown.(sessLab).ICI)];
            ICD_all = [ICD_all, (MI_breakdown.(sessLab).ICD)];

            Lin_all = [Lin_all, MI_pop.(sessLab).linear];
            RBF_all = [RBF_all, MI_pop.(sessLab).RBF];

            LinSh_all = [LinSh_all; MISh_pop.(sessLab).linear];
            RBFSh_all = [RBFSh_all; MISh_pop.(sessLab).RBF];

            % PID
            PIDJoint_all = [PIDJoint_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,1,:)), 1, []))];
            Syn_all = [Syn_all, cell2mat(reshape(squeeze(PID_result.(sessLab)(:,2,:)), 1, []))];
            Red_all = [Red_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,3,:)), 1, []))];
            Unq1_all = [Unq1_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,4,:)), 1, []))];
            Unq2_all = [Unq2_all,cell2mat(reshape(squeeze(PID_result.(sessLab)(:,5,:)), 1, []))];
        end
        LinSh_mean = (LinSh_all(:,1:5) + LinSh_all(:,6:10)) / 2;
        RBFSh_mean = (RBFSh_all(:,1:5) + RBFSh_all(:,6:10)) / 2;
        LinSh_mean = LinSh_mean(:)';
        RBFSh_mean = RBFSh_mean(:)';
        LinSh_all = LinSh_all(:)';
        RBFSh_all = RBFSh_all(:)';

    end
    % ---------------------------------------------------------------------
    % Step 2: Compute mean, SEM and p-value
    % ---------------------------------------------------------------------
    Co_info_all =  Joint_all - ILIN_all;

    ratios = RBF_all ./ Lin_all;
    mean_ratio1 = mean(ratios);

    ratios = Lin_all ./ RBF_all;
    mean_ratio2 = mean(ratios);

    ratios = Lin_all ./LinSh_mean;
    mean_ratioshLin = mean(ratios);

    ratios = RBF_all./RBFSh_mean;
    mean_ratioshRBF = mean(ratios);

    [meanJoint, errJoint] = compute_stats(Joint_all);
    [meanPIDJoint, errPIDJoint] = compute_stats(PIDJoint_all);

    [meanILin, errILin] = compute_stats(ILIN_all);
    [meanISS, errISS] = compute_stats(ISS_all);
    [meanICI, errICI] = compute_stats(ICI_all);
    [meanICD, errICD] = compute_stats(ICD_all);
    meanISS = (meanISS/meanJoint)*100;
    errISS = (errISS/meanJoint)*100;
    meanICI = (meanICI/meanJoint)*100;
    errICI = (errICI/meanJoint)*100;
    meanICD = (meanICD/meanJoint)*100;
    errICD = (errICD/meanJoint)*100;

    [meanRBF, errRBF] = compute_stats(RBF_all);
    [meanLin, errLin] = compute_stats(Lin_all);
    [meanRBFSh, errRBFSh] = compute_stats(RBFSh_all);
    [meanLinSh, errLinSh] = compute_stats(LinSh_all);

    [meanSyn, errSyn] = compute_stats((Syn_all));
    [meanRed, errRed] = compute_stats((Red_all));
    meanSyn = (meanSyn/meanPIDJoint)*100;
    errSyn = (errSyn/meanPIDJoint)*100;
    meanRed = (meanRed/meanPIDJoint)*100;
    errRed = (errRed/meanPIDJoint)*100;

    [meanCo, errCo] = compute_stats(Co_info_all);


    [~, p_Lin_RBF] = ttest(Lin_all, RBF_all);
    [~, p_RBF_RBFSh] = ttest(RBF_all, RBFSh_mean);
    [~, p_Lin_LinSh]   = ttest(Lin_all,  LinSh_mean);
    [~, p_Lin_Joint] = ttest(ILIN_all, Joint_all);
    [~, p_Red_Syn] = ttest(Red_all, Syn_all);

    % ---------------------------------------------------------------------
    % Step 3: Plot column
    % ---------------------------------------------------------------------

    scIdx = 1;
    %Joint vs. Linear Information Panel 1
    xShift = 0;
    fSize = 12;
    minYText = 0.03;
    scIdx = 1;

    nexttile(tileIndex);
    hold on
    if col == 1
        title('Pair Info');
    end
    h{1} = bar(scIdx,meanJoint,'FaceColor',colorJoint);
    errorbar(scIdx, meanJoint,errJoint,'k--');
    h{2} = bar(scIdx+1,meanILin,'FaceColor',colorLin_2);
    errorbar(scIdx+1,meanILin,errILin,'k--')
    h{3} = bar(scIdx+2,meanCo,'FaceColor',colorJoint);
    errorbar(scIdx+2, meanCo,errCo,'k--');
    ylabel('Information [bits]', 'FontSize', fontsize)
    xticks([scIdx scIdx+1 scIdx+2])
    set(gca, 'XTickLabel', {}, 'FontSize', fontsize)
    if col == 2
        maxY = 0.3;
        minY = -0.05;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.005, 0.27)
    elseif col == 1
        maxY = 0.08;
        minY = -0.0133;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.0025, 0.065)
    elseif col == 4
        maxY = 0.3;
        minY = -0.05;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.005, 0.25)
    elseif col == 3
        maxY = 0.15;
        minY = -0.025;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.005, 0.12)
    end


    ylim([minY,maxY])

    xticklabels({'Joint', 'Sum', 'Co-Info'})

    hold off

    % PID
    row = 2;
    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on
    if col == 1
        title('PID');
    end
    h{1} = bar(scIdx,meanSyn,'FaceColor',colorSyn);
    errorbar(scIdx, meanSyn,errSyn,'k--')
    h{2} = bar(scIdx+1,meanRed,'FaceColor',colorRed);
    errorbar(scIdx+1,meanRed,errRed,'k--')
    ylabel('% of Joint Information', 'FontSize', fontsize)
    xticks([scIdx scIdx+1])
    set(gca, 'XTickLabel', {}, 'FontSize', fontsize)
    if col ==2
        maxY = 60;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 2, 55)
    elseif col ==1
        maxY = 110;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 2, 105)
    elseif col == 4
        maxY = 50;
        add_pvalue(p_Red_Syn, scIdx-0.25, scIdx+1.25, 1, 45)
    elseif col == 3
        maxY = 30;
        add_pvalue(p_Red_Syn, scIdx-0.25, scIdx+1.25, 1, 20)
    end
    minY = 0;
    ylim([minY,maxY])
    xticklabels({'Syn', 'Red'})
    hold off

    % Information Breakdown
    row = 3;
    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on;
    if col == 1
        title('Info Breakdown');
    end
    ylabel('% of Joint Information', 'FontSize', fontsize)
    % Sig. sim.
    h{1} = bar(scIdx,meanISS,'FaceColor',colorInfoBreak);
    errorbar(scIdx,meanISS,errISS,'k--')
    % Cor-indep
    h{2} = bar(scIdx+1,meanICI,'FaceColor',colorInfoBreak);
    errorbar(scIdx+1,meanICI,errICI,'k--')
    % Cor-dep
    h{3} = bar(scIdx+2,meanICD,'FaceColor',colorInfoBreak);
    errorbar(scIdx+2,meanICD,errICD,'k--')
    if col ==2
        maxY = 10;
    elseif col ==1
        maxY = 110;
    elseif col == 4
        maxY = 10;
    elseif col == 3
        maxY = 10;
    end
    minY = -12;
    xticks([scIdx scIdx+1 scIdx+2])
    set(gca, 'XTickLabel', {}, 'FontSize', fontsize)
    ylim([minY,maxY])
    xticklabels({'ISS', 'ICI', 'ICD'});
    hold off;

    scattersize = 10;
    row = 4;
    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on;
    ylabel('Info SVM linear', 'FontSize', fontsize)
    xlabel('Info SVM RBF', 'FontSize', fontsize)
    if col ==1 || col == 3
        scatter(RBF_all, Lin_all,  scattersize,'filled');
        if col ==1
            ymax = 0.6;
            xmax = 0.6;
            ymin = -0.05;
            xmin = -0.05;
        else
            ymax = 1.4;
            xmax = 1.4;
            ymin = 0;
            xmin = 0;
        end
        plot([ymin, ymax], [xmin, xmax], 'k--', 'LineWidth', 1);
    elseif col ==2 || col == 4
        scatter(RBF_all, Lin_all, scattersize, 'filled');
        ymax = 1;
        xmax = 1;
        ymin = 0;
        xmin = 0;
        plot([ymin, ymax], [xmin, xmax], 'k--', 'LineWidth', 1);
    end
    ylim([ymin, ymax]);
    xlim([xmin, xmax]);
    axis square;
    hold off
    row = 5;

    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on;
    ylabel('Info SVM orig. data', 'FontSize', fontsize)
    xlabel('Info SVM shuff. data', 'FontSize', fontsize)
    if col ==1 || col == 3
        scatter(RBFSh_mean, RBF_all, scattersize, 'filled');
        if col ==1
            ymax = 0.6;
            xmax = 0.6;
            ymin = -0.05;
            xmin = -0.05;
        else
            ymax = 1.4;
            xmax = 1.4;
            ymin = 0;
            xmin = 0;
        end
        plot([ymin, ymax], [xmin, xmax], 'k--', 'LineWidth', 1);
    elseif col ==2 || col == 4
        scatter(LinSh_mean, Lin_all, scattersize, 'filled');
        ymax = 1.2;
        xmax = 1.2;
        ymin = 0;
        xmin = 0;
        plot([ymin, ymax], [xmin, xmax], 'k--', 'LineWidth', 1);
    end
    ylim([ymin, ymax]);
    xlim([xmin, xmax]);
    axis square;
    hold off
end

% saveas(gcf, 'Figure2/Plots/Figure2_Draft.svg')
% saveas(gcf, 'Figure2/Plots/Figure2_Draft.jpg')
function [mean_val, sem_val] = compute_stats(data)
mean_val = mean(data);
std_val = std(data);
n = length(data);
sem_val = std_val / sqrt(n);
end

function add_pvalue(p_values, x1, x2, text_heights, line_heights)
for i = 1:numel(p_values)
    stars = get_stars(p_values(i));
    line([x1(i), x2(i)], [line_heights(i),line_heights(i)], 'Color', 'black', 'LineWidth', 1);
    text((x1(i) + x2(i)) / 2, line_heights(i) + text_heights, stars, 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
end
end

function stars = get_stars(p_value)
if p_value < 0.001
    stars = '***';
elseif p_value < 0.01
    stars = '**';
elseif p_value < 0.05
    stars = '*';
else
    stars = 'n.s';
end
end