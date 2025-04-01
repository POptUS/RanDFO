%problem = 'rosenbrock';
problem = 'cube';
test_types = {'balanced','progressive','imbalanced'};
%m = 16;
m = 64;
FS = 16;
seeds_to_use = 3;

figure; hlt = tiledlayout(1,3);

for j = 1:length(test_types)
    test_type = test_types{j};
    for macro_seed = 1:seeds_to_use
        filename = strcat('results/pounders_compare_',problem,'_',test_type,'_',num2str(macro_seed),'_',num2str(m),'.mat');
        load(filename);
        if macro_seed == 1
            H = [Hf(2,:); Hf(1,:)];
        else
            H = cat(1, H, Hf(1,:));
        end    
    end
    nexttile;
    [hl, hlb] = plot_trajectory_percentiles2(H,40);
    if strcmp(test_type, 'balanced')
        legend(hl,{'POUNDERS', 'SAM-POUNDERS'},'Location','SouthWest','FontSize',FS)
    end
    if strcmp(test_type, 'balanced')
        ylabel('$f(\mathbf{x}^k) - f(\mathbf{x}^*)$','interpreter','latex','FontSize',FS);
    end
    if strcmp(test_type, 'progressive')
        xlabel('Component function evaluations','FontSize',FS);
    end

end
    hlt.Padding = 'compact'; hlt.TileSpacing = 'compact'; 
    pause() % MAKE IT CAMERA READY DURING THE PAUSE
    outputname = strcat('images/pounders_compare_',problem,'_',num2str(m),'.eps');
    saveas(gcf,outputname,'epsc');
