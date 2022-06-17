
function step2_smooth(batch, roi_confidence)
    if batch == '1'
        step2_smooth_batch(roi_confidence)
    else
        step2_smooth_select(roi_confidence)
    end
end

function step2_smooth_batch(roi_confidence)
    dirinfo = dir('data_step1/*h_ecc_cut.m');

    mkdir('data_step2/');
    mkdir('data_step2/smooth/');
    mkdir('data_step2/unsmooth/');
    mkdir('data_step2/figures/');

    for i=1:length(dirinfo)
        infn = ['data_step1/' dirinfo(i).name];

        %
        % output filename appended with roi_confidence value
        resultfn_sm = ['data_step2/smooth/' strrep(dirinfo(i).name, '_ecc_cut.m', ['_prf_smooth_' num2str(roi_confidence) '.m'])];
        resultfn_usm = ['data_step2/unsmooth/' strrep(dirinfo(i).name, '_ecc_cut.m', ['_prf_unsmooth_' num2str(roi_confidence) '.m'])];

        resultsfn_sm_filename = strrep(dirinfo(i).name, '_ecc_cut.m', '_prf_smooth.m');

        resultdir_sv_figs = 'data_step2/figures/';

        % subject 169444 has very low pRF rsquared confidence at 5% still only
        % has a small amount of data for smoothing.
        % NOTE: results for this subject can be discarded
        if contains(infn,'169444')
            try
                run_step2(infn, resultfn_sm, resultsfn_sm_filename, resultfn_usm, resultdir_sv_figs, 5)
            catch
                errorfn = [dirinfo(i).name '.err'];
                fid = fopen(errorfn, 'w'); fclose(fid); 
            end
            disp('NOTE: Subject 169444 has very loc pRF variance explained data: consider discarding');
            %continue % skip this subject (roi_confidence > 25 == 0)
        else
            try
                run_step2(infn, resultfn_sm, resultsfn_sm_filename, resultfn_usm, resultdir_sv_figs, roi_confidence)
            catch
                errorfn = [dirinfo(i).name '.err'];
                fid = fopen(errorfn, 'w'); fclose(fid); 
            end
        end  
    end
end

function step2_smooth_select(roi_confidence)  
	[file,path] = uigetfile('*.m');
    if isequal(file,0)
        disp('User selected Cancel');
        return
    else
        disp(['User selected ', fullfile(path,file)]);

        mkdir('data_step2/');
        mkdir('data_step2/smooth/');
        mkdir('data_step2/unsmooth/');
        mkdir('data_step2/figures/');
        
        infn = fullfile(path,file);

        %
        % output filename appended with roi_confidence value
        resultfn_sm = ['data_step2/smooth/' strrep(file, '_ecc_cut.m', ['_prf_smooth_' num2str(roi_confidence) '.m'])];
        resultfn_usm = ['data_step2/unsmooth/' strrep(file, '_ecc_cut.m', ['_prf_unsmooth_' num2str(roi_confidence) '.m'])];

        resultsfn_sm_filename = strrep(file, '_ecc_cut.m', '_prf_smooth.m');

        resultdir_sv_figs = 'data_step2/figures/';
        
        % subject 169444 has very low pRF rsquared confidence at 5% still only
        % has a small amount of data for smoothing.
        % NOTE: results for this subject can be discarded
        if contains(infn,'169444')
            try
                run_step2(infn, resultfn_sm, resultsfn_sm_filename, resultfn_usm, resultdir_sv_figs, 5)
            catch
                errorfn = [file '.err'];
                fid = fopen(errorfn, 'w'); fclose(fid); 
            end
            disp('NOTE: Subject 169444 has very loc pRF variance explained data: consider discarding');
            %continue % skip this subject (roi_confidence > 25 == 0)
        else
            try
                run_step2(infn, resultfn_sm, resultsfn_sm_filename, resultfn_usm, resultdir_sv_figs, roi_confidence)
            catch
                errorfn = [file '.err'];
                fid = fopen(errorfn, 'w'); fclose(fid); 
            end
        end
    end
end

function run_step2(infn, outfn_sm, outfilename, outfn_usm, out_sv_figs, in_roi_confidence)
    compute_prf_smoothing(infn, outfn_sm, outfilename, outfn_usm, out_sv_figs, in_roi_confidence)
end