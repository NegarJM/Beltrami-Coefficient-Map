function step3_compute_bcm(batch, roi_confidence)
    if batch == '1'
        step3_compute_bcm_batch(roi_confidence)
    else
        step3_compute_bcm_select(roi_confidence)
    end
end

function step3_compute_bcm_batch(roi_confidence)
    dirinfo_smooth = dir('data_step2/smooth/*h_prf_smooth*.m');
    dirinfo_unsmooth = dir('data_step2/unsmooth/*h_prf_unsmooth*.m');

    mkdir('data_step3/');
    mkdir('data_step3/figures/');

    for i=1:length(dirinfo_smooth)
        infn_sm = ['data_step2/smooth/' dirinfo_smooth(i).name];
        infn_us = ['data_step2/unsmooth/' dirinfo_unsmooth(i).name]; 

        input_str = dirinfo_smooth(i).name;
        regexpr = '_prf_smooth_(\w*).m';
        replace_str = '_prf_smooth.m';
        resultsfn_sm_filename = regexprep(input_str, regexpr, replace_str);

        resultfn_left_bcm_stats_ecc = ['data_step3/' 'results_left_beltrami_stats_eccentricity.txt'];
        resultfn_left_bcm_stats_ang = ['data_step3/' 'results_left_beltrami_stats_polarangle.txt'];
        resultfn_right_bcm_stats_ecc = ['data_step3/' 'results_right_beltrami_stats_eccentricity.txt'];
        resultfn_right_bcm_stats_ang = ['data_step3/' 'results_right_beltrami_stats_polarangle.txt'];

        resultfn_sv_figs = 'data_step3/figures/';

        try
            run_step3(infn_sm, infn_us, resultsfn_sm_filename, resultfn_sv_figs, resultfn_left_bcm_stats_ecc, resultfn_left_bcm_stats_ang, resultfn_right_bcm_stats_ecc, resultfn_right_bcm_stats_ang, roi_confidence)
        catch
            errorfn = [dirinfo_smooth(i).name '.err'];
            fid = fopen(errorfn, 'w'); fclose(fid); 
            errorfn = [dirinfo_unsmooth(i).name '.err'];
            fid = fopen(errorfn, 'w'); fclose(fid); 
        end
    end
end

function step3_compute_bcm_select(roi_confidence)
	[file,path] = uigetfile('*.m','Select smoothed pRF file');

    if isequal(file,0)
        disp('User selected Cancel');
        return
    else
        disp(['User selected ', fullfile(path,file)]);

        mkdir('data_step3/');
        mkdir('data_step3/figures/');   

        infn_sm = ['data_step2/smooth/' file];
        infn_us = ['data_step2/unsmooth/' strrep(file, 'prf_smooth', 'prf_unsmooth')]; 

        input_str = file;
        regexpr = '_prf_smooth_(\w*).m';
        replace_str = '_prf_smooth.m';
        resultsfn_sm_filename = regexprep(input_str, regexpr, replace_str);

        resultfn_left_bcm_stats_ecc = ['data_step3/' 'results_left_beltrami_stats_eccentricity.txt'];
        resultfn_left_bcm_stats_ang = ['data_step3/' 'results_left_beltrami_stats_polarangle.txt'];
        resultfn_right_bcm_stats_ecc = ['data_step3/' 'results_right_beltrami_stats_eccentricity.txt'];
        resultfn_right_bcm_stats_ang = ['data_step3/' 'results_right_beltrami_stats_polarangle.txt'];

        resultfn_sv_figs = 'data_step3/figures/';

        try
            run_step3(infn_sm, infn_us, resultsfn_sm_filename, resultfn_sv_figs, resultfn_left_bcm_stats_ecc, resultfn_left_bcm_stats_ang, resultfn_right_bcm_stats_ecc, resultfn_right_bcm_stats_ang, roi_confidence)
        catch
            errorfn = [file '.err'];
            fid = fopen(errorfn, 'w'); fclose(fid); 
            errorfn = [strrep(file, 'prf_smooth', 'prf_unsmooth') '.err'];
            fid = fopen(errorfn, 'w'); fclose(fid); 
        end
    end
end

function run_step3(infn_s, infn_u, outfn_sm_filename, outfn_sv_figs, outfn_left_ecc, outfn_left_ang, outfn_right_ecc, outfn_right_ang, in_roi_confidence)
    if contains(infn_s,'lh')
        LoR = 0;
    else
        LoR = 1;
    end

    generate_bcm_results(infn_s, infn_u, outfn_sm_filename, outfn_sv_figs, in_roi_confidence, outfn_left_ecc, outfn_left_ang, outfn_right_ecc, outfn_right_ang)
end