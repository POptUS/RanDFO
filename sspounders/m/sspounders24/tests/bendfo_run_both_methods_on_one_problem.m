function bendfo_run_both_methods_on_one_problem(prob_no)

	prob_no = str2num(prob_no);

    [objective, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver, hfun, combinemodels] = ... 
        setup_bendfo_test(prob_no);

    Prior = []; 

    Options.hfun = hfun;
    Options.combinemodels = combinemodels; 
    Options.printf = 1; 

    seeds = 1;

        % run pounders first
    filename = strcat('results/bendfo_pounders_',num2str(prob_no),'.mat');
    %if ~exist(filename)
    [~, ~, hF, ~, xkin] = pounders(objective, X0, n, nfmax, gtol, delta, m, L, U, Prior, Options);
    %save pounders result
    hF = hF(1:xkin);
    save(filename, 'hF');
    %end
    
% only ucb
    ucb_sketchsize = 1;
    random_sketchsize = 0; 
    memory = n; 
    for seed = 1:seeds    
	    savestr = strcat('results/bendfo_sspounders_',num2str(prob_no),'_',num2str(seed),'.mat');
	    rng(seed);
        %if ~exist(savestr)
            [~, ~, hF, ~, xk_in] = subspace_pounders(objective, X0, n, nfmax, gtol, delta, m, random_sketchsize, ucb_sketchsize, memory, L, U, Prior, Options);
            hF = hF(1:xk_in);     
            save(savestr, 'hF');
        %end
    end


    % no ucb
%     ucb_sketchsize = 0;
%     random_sketchsize = 2; 
%     memory = n; 
%     for seed = 1:seeds    
% 	    savestr = strcat('results/bendfo_sspounders_',num2str(prob_no),'_',num2str(seed),'.mat');
% 	    rng(seed);
%         %if ~exist(savestr)
%             [~, ~, hF, ~, xk_in] = subspace_pounders(objective, X0, n, nfmax, gtol, delta, m, random_sketchsize, ucb_sketchsize, memory, L, U, Prior, Options);
%             hF = hF(1:xk_in);     
%             save(savestr, 'hF');
%         %end
%     end
% 
%         % with ucb
%     ucb_sketchsize = 1;
%     random_sketchsize = 1; 
%     memory = n; 
%     for seed = 1:seeds    
% 	    savestr = strcat('results/bendfo_sspounders_',num2str(prob_no),'_',num2str(seed),'.mat');
% 	    rng(seed);
%         %if ~exist(savestr)
%             [~, ~, hF, ~, xk_in] = subspace_pounders(objective, X0, n, nfmax, gtol, delta, m, random_sketchsize, ucb_sketchsize, memory, L, U, Prior, Options);
%             hF = hF(1:xk_in);     
%             save(savestr, 'hF');
%         %end
%     end
end
