function run_both_methods_on_one_problem(prob_no)

    seeds = 3;

	prob_no = str2num(prob_no);

    [objective, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver, hfun, combinemodels] = ... 
        setup_yatsop_test('midscale',prob_no);

    
    % run pounders first
    filename = strcat('results/pounders_',num2str(prob_no),'.mat');
    if ~exist(filename)
    [~, ~, ~, xkin, Fs] = pounders(objective, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver, hfun, combinemodels);
    % save pounders result
    Fs = Fs(1:xkin);
    save(filename, 'Fs');
    end

    % now run sspounders with each of 30 seeds
    sketchsize = 1; % this will be ignored in this test if adaptive=true
    sketchtype = 'grad';
    adaptive = true; 
    full_iter = true; % check this. 
    for seed = 1:seeds    
	xkin = 1;
	filename = strcat('results/sspounders_',num2str(prob_no),'_',num2str(seed),'.mat');
	rng(seed);
	%if ~exist(filename)
		try
        	[~, ~, ~, xkin, Fs] = sspounders3(objective, X0, n, npmax, nfmax, gtol, delta, nfs, m, F0, xkin, L, U, printf, spsolver, hfun, combinemodels, sketchsize, sketchtype, adaptive, full_iter);
        	Fs = Fs(1:xkin);
        	save(filename, 'Fs');
		catch
			fprintf('This run did not complete correctly. \n')
		end
	%end
    end
end
