clc, clear all, close all;



plotRows = 6;
plotColumns = 4;
width_cm = 19.05; 
height_ratio = 6 / 4;  
height_cm = width_cm * height_ratio;
width_inch = width_cm / 2.54;  
height_inch = height_cm / 2.54; 
figure_handle = figure('Units', 'inches', 'Position', [1, 1, width_inch, height_inch]);
t = tiledlayout(6, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
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
        load('Results/Simulation_Results_200_10.mat');

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

        LinSh_all2 = [];
        RBFSh_all2 = [];
        LinSh_all2_2 = [];
        RBFSh_all2_2 = [];
        for i = 1:nShuff
            LinSh_all2 = [LinSh_all2, MI_popSh.linear.(scLab)(:,i).all];
            RBFSh_all2 = [RBFSh_all2, MI_popSh.RBF.(scLab)(:,i).all];
            LinSh_all2_2 = [LinSh_all2_2; MI_popSh.linear.(scLab)(:,i).all];
            RBFSh_all2_2 = [RBFSh_all2_2; MI_popSh.RBF.(scLab)(:,i).all];
        end
        LinSh_all = mean(LinSh_all2,1);
        RBFSh_all = mean(RBFSh_all2,1);
        Lin_all =   horzcat(MI_pop.linear.(scLab)(:).all);
        RBF_all =   horzcat(MI_pop.RBF.(scLab)(:).all);
        LinSh_mean = mean(LinSh_all2_2,1);
        RBFSh_mean = mean(RBFSh_all2_2,1);

        PIDJoint_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,1,:)), 1, []));
        Syn_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,2,:)), 1, []));
        Red_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,3,:)), 1, []));
        Unq1_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,4,:)), 1, []));
        Unq2_all = cell2mat(reshape(squeeze(PID_result.(scLab)(:,5,:)), 1, []));

    elseif col == 4
        load('Results/FrancisData_Results_btsp_naive.mat');

        LinSh_all2_2 = [];
        RBFSh_all2_2 = [];
        for sessIdx = 1:numel(field_names)
            sessLab = field_names{sessIdx};
            Joint_all = [Joint_all, (MI_breakdown.(sessLab).Joint)'];
            ILIN_all = [ILIN_all, (MI_breakdown.(sessLab).ILIN)'];
            ISS_all = [ISS_all, (MI_breakdown.(sessLab).ISS)'];
            ICI_all = [ICI_all, (MI_breakdown.(sessLab).ICI)'];
            ICD_all = [ICD_all, (MI_breakdown.(sessLab).ICD)'];

            % Population Information
            for i = 1:nShuff
                LinSh_all  = [LinSh_all, MISh.linear.(sessLab)(:,i).all];
                RBFSh_all = [ RBFSh_all, MISh.RBF.(sessLab)(:,i).all];
                LinSh_all2_2 = [LinSh_all2_2; MISh.linear.(sessLab)(:,i).all];
                RBFSh_all2_2 = [RBFSh_all2_2;  MISh.RBF.(sessLab)(:,i).all];
            end
            Lin_all = [Lin_all,MI.linear.(sessLab).all];
            RBF_all = [RBF_all,MI.RBF.(sessLab).all];


            % PID
            Syn_all = [Syn_all, (PID_result.(sessLab).corrected(:,4))'];
            Red_all = [Red_all, (PID_result.(sessLab).corrected(:,1))'];
            Unq1_all = [Unq1_all, (PID_result.(sessLab).corrected(:,2))'];
            Unq2_all = [Unq2_all, (PID_result.(sessLab).corrected(:,3))'];
            PIDJoint_all = [PIDJoint_all, (PID_result.(sessLab).corrected(:,5))'];
        end
        LinSh_mean = mean(LinSh_all2_2',1);
        RBFSh_mean = mean(RBFSh_all2_2',1);

        ratios = RBF_all ./ Lin_all;
        mean_ratio = mean(ratios);

    elseif col == 3
        load('Results/HippocampusData_Results_btsp_naive_6bins.mat');
        LinSh_tmp = [];
        RBFSh_tmp = [];
        LinSh_all = [];
        RBFSh_all = [];
        for sessIdx = 1:numel(field_names)
            sessLab = field_names{sessIdx};
            Joint_all = [Joint_all, (MI_breakdown.(sessLab).Joint)'];
            ILIN_all = [ILIN_all, (MI_breakdown.(sessLab).ILIN)'];
            ISS_all = [ISS_all, (MI_breakdown.(sessLab).ISS)'];
            ICI_all = [ICI_all, (MI_breakdown.(sessLab).ICI)'];
            ICD_all = [ICD_all, (MI_breakdown.(sessLab).ICD)'];

            LinSh_all2 = [];
            RBFSh_all2 = [];
            for i = 1:nShuff
                LinSh_all2  = [LinSh_all2; MISh.linear.(sessLab)(:,i).all];
                RBFSh_all2 = [ RBFSh_all2; MISh.RBF.(sessLab)(:,i).all];
                LinSh_all = [LinSh_all, MISh.linear.(sessLab)(:,i).all];
                RBFSh_all = [RBFSh_all,  MISh.RBF.(sessLab)(:,i).all];
            end
            Lin_all = [Lin_all,MI.linear.(sessLab).all];
            RBF_all = [RBF_all,MI.RBF.(sessLab).all];
            LinSh_tmp = [LinSh_tmp, LinSh_all2];
            RBFSh_tmp = [RBFSh_tmp, RBFSh_all2];
            % PID
            Syn_all = [Syn_all, (PID_result.(sessLab).corrected(:,4))'];
            Red_all = [Red_all, (PID_result.(sessLab).corrected(:,1))'];
            Unq1_all = [Unq1_all, (PID_result.(sessLab).corrected(:,2))'];
            Unq2_all = [Unq2_all, (PID_result.(sessLab).corrected(:,3))'];
            PIDJoint_all = [PIDJoint_all, (PID_result.(sessLab).corrected(:,5))'];
        end
        LinSh_mean = mean(LinSh_tmp,1);
        RBFSh_mean = mean(RBFSh_tmp,1);

    end
    Co_info_all =  Joint_all - ILIN_all;

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
    % meanCo = (meanCo/meanJoint)*100;
    % errCo = (errCo/meanJoint)*100;
    % p_values

    [~, p_Lin_RBF] = ttest(Lin_all, RBF_all);
    [~, p_RBF_RBFSh] = ttest(RBF_all, RBFSh_mean);
    [~, p_Lin_LinSh]   = ttest(Lin_all,  LinSh_mean);
    [~, p_Lin_Joint] = ttest(ILIN_all, Joint_all);
    [~, p_Red_Syn] = ttest(Red_all, Syn_all);

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
    % ylabel('Information [bits]', 'FontSize', fontsize)
    xticks([scIdx scIdx+1 scIdx+2])
    set(gca, 'XTickLabel', {}, 'FontSize', fontsize)
    if col ==2
        maxY = 0.3;
        minY = -0.05;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.005, 0.27)
    elseif col ==1
        maxY = 0.08;
        minY = 0;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.0025, 0.065)
    elseif col == 4
        maxY = 0.3;
        minY = -0.05;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.005, 0.25)
    elseif col == 3
        maxY = 0.15;
        minY = 0;
        add_pvalue(p_Lin_Joint, scIdx-0.25, scIdx+1.25, 0.005, 0.12)
    end

    ylim([minY,maxY])
    if col == 4
        %xticklabels({'Joint', 'Sum'})
    end
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
    if col == 4
        %xticklabels({'Syn', 'Red'})
    end
    hold off

    % Information Breakdown
    row = 3;
    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on;
    if col == 1
        title('Info Breakdown');
    end
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
        maxY = 100;
    elseif col == 4
        maxY = 10;
    elseif col == 3
        maxY = 10;
    end
    minY = -12;
    xticks([scIdx scIdx+1 scIdx+2])
    set(gca, 'XTickLabel', {}, 'FontSize', fontsize)
    ylim([minY,maxY])
    if col == 4
        %xticklabels({'ISS', 'ICI', 'ICD'});
    end
    hold off;


    % %Population Information
    % nexttile;
    % hold on;
    % h{1} = bar(scIdx,meanLin,'FaceColor',colorLin);
    % errorbar(scIdx,meanLin,errLin,'k--')
    % h{2} = bar(scIdx+1,meanLinSh,'FaceColor',colorLinSh);
    % errorbar(scIdx+1,meanLinSh,errLinSh,'k--')
    % h{3} = bar(scIdx+2,meanRBF,'FaceColor',colorRBF);
    % errorbar(scIdx+2,meanRBF,errRBF,'k--')
    % h{4} = bar(scIdx+3,meanRBFSh,'FaceColor',colorRBFSh);
    % errorbar(scIdx+3,meanRBFSh,errRBFSh,'k--')
    % xticks([scIdx scIdx+1 scIdx+2 scIdx+3 scIdx+4])
    % set(gca, 'XTickLabel', {}, 'FontSize', fontsize)
    % ylabel('Information [bits]', 'FontSize', fontsize)
    %
    % if col ==2
    %     maxY = 1.2;
    %     add_pvalue([p_Lin_LinSh, p_RBF_RBFSh, p_Lin_RBF], [1, 3, 1]-.25, [1, 3, 2]+1.25, 0.04, [0.95,  0.95, 1.1]);
    % elseif col ==1
    %     maxY = 0.8;
    %     add_pvalue([p_Lin_LinSh, p_RBF_RBFSh, p_Lin_RBF], [1, 3, 1]-.25, [1, 3, 2]+1.25, 0.04, [0.5,  0.5, 0.7]);
    % elseif col == 4
    %     maxY = 0.8;
    %     add_pvalue([p_Lin_LinSh, p_RBF_RBFSh, p_Lin_RBF], [1, 3, 1]-.25, [1, 3, 2]+1.25, 0.02, [0.62,  0.62, 0.72]);
    % elseif col == 3
    %     maxY = 1.4;
    %     add_pvalue([p_Lin_LinSh, p_RBF_RBFSh, p_Lin_RBF], [1, 3, 1]-.25, [1, 3, 2]+1.25, 0.02, [1,  1.1, 1.25]);
    % end
    % if col == 4
    %     xticklabels({'Linear', 'Linear shuff.', 'RBF', 'RBF shuff.'})
    % end
    % minY = 0;
    % ylim([minY, maxY])
    % hold off;
    scattersize = 10;
    row = 4;
    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on;
    if col ==1 || col == 3
        scatter(RBFSh_mean, RBF_all, scattersize, 'filled');
        % xlabel('Info SVM RBF');
        % ylabel('Info SVM RBFSh');
              
        if col ==1
            ymax = 0.6;
            xmax = 0.6;
        else
            ymax = 1.4;
            xmax = 1.4;
        end
        plot([0, ymax], [0, xmax], 'k--', 'LineWidth', 2);
    elseif col ==2 || col == 4
        scatter(LinSh_mean, Lin_all, scattersize, 'filled');
        % xlabel('Info SVM Lin');
        % ylabel('Info SVM LinSh');
        ymax = 1.2;
        xmax = 1.2;
        plot([0, ymax], [0, xmax], 'k--', 'LineWidth', 2);
    end
    ylim([0, ymax]);
    xlim([0, xmax]);
     axis equal;
    hold off

    row = 5;
    tileIndex=(row-1)*plotColumns+col;
    nexttile(tileIndex);
    hold on;
    if col ==1 || col == 3
        scatter(Lin_all, RBF_all,  scattersize,'filled');
        % xlabel('Info SVM RBF');
        % ylabel('Info SVM Lin');
        if col ==1
            ymax = 0.6;
            xmax = 0.6;
        else
            ymax = 1.4;
            xmax = 1.4;
        end
        plot([0, ymax], [0, xmax], 'k--', 'LineWidth', 2);
    elseif col ==2 || col == 4
        scatter(RBF_all, Lin_all, scattersize, 'filled');
        % xlabel('Info SVM Lin');
        % ylabel('Info SVM RBF');
        ymax = 1;
        xmax = 1;
        plot([0, ymax], [0, xmax], 'k--', 'LineWidth', 2);
    end
    ylim([0, ymax]);
    xlim([0, xmax]);
     axis equal;
    hold off
end


if ~exist('Results', 'dir')
    mkdir('Results');
end
saveas(gcf, 'Results/Figs/Figure2_200Trials_10Reps.svg')

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