function step_all(batch, roi_confidence)
    if batch == '1'             % Input = '1' batch mode 
        step_all_batch(roi_confidence)  % Call the Cut and Flattening function for batch mode
    else                        % Input = '0' single selected subject
        step_all_select(roi_confidence) % Call the Cut and Flattening Function for one selected single subject
    end
end

function step_all_batch(roi_confidence)
    %cut_flat
    dirinfo = dir('data_step0/*_ecc.m');  
    mkdir('data_step1');
    mkdir('data_step0_subdiv');
    for i=1:length(dirinfo)
        infn = ['data_step0/' dirinfo(i).name ];
        
        cutfn = ['data_step1/' dirinfo(i).name(1:end-2) '_cut.m']; 
        subdivfn = ['data_step0_subdiv/' dirinfo(i).name(1:end-2) '_subdiv.m']; 

        try 
            run_step1(infn, cutfn, subdivfn); % Call the main function for cut and flattening 
        catch
            errorfn = [dirinfo(i).name '.err'];
            fid = fopen(errorfn, 'w'); fclose(fid); 
        end
    end
    %smooth
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
    
    %bcm_computation
    
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

function step_all_select(roi_confidence)
    [file,path] = uigetfile('*.m'); % Ask the user to select a data
    if isequal(file,0)
        disp('User selected Cancel');
        return
    else
        disp(['User selected ', fullfile(path,file)]);

        mkdir('data_step1');
        mkdir('data_step0_subdiv');

        infn = fullfile(path,file);
        
        cutfn = ['data_step1/' file(1:end-2) '_cut.m']; 
        subdivfn = ['data_step0_subdiv/' file(1:end-2) '_subdiv.m']; 

        try 
            run_step1(infn, cutfn, subdivfn); % Call the main function for cut and flattening 
        catch
            errorfn = [file '.err'];
            fid = fopen(errorfn, 'w'); fclose(fid); 
        end
    end
    
        file = [file(1:end-2) '_cut.m'];
        path = [path(1:end-11) 'data_step1\'];
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

        
        file = strrep(file, '_ecc_cut.m', ['_prf_smooth_' num2str(roi_confidence) '.m']);
        path = [path(1:end-11) 'data_step2\smooth\'];
        infn = fullfile(path,file);
        
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

function run_step1(infn, outfn, smoothfn)  
    [F,V,E]=read_mfile(infn);  % load the Faces, Vertices and Edges 
    options.interp = 'linear'; % for vertices only, use loop
    Vmid = 0.5*(E.Vertex_Vpial - E.Vertex_Vwhite) + E.Vertex_Vwhite;
    [Vs,Fs] = perform_mesh_subdivision(Vmid', F', options); 
    Vs = Vs';
    Fs = Fs';

    attributes = fieldnames(E);
    for k=1:numel(attributes)
        ak  = attributes{k};
        if( isnumeric(E.(ak)) )
            field = E.(ak);
            if size(field,1) == size(V,1)
                if all(round(field) == field)
                    options.interp = 'nearest'; % typically for interger index
                else
                    options.interp = 'linear'; % typically for interger index
                end
                if ak == "Vertex_prf"
                    prf = field;
                    xy = [prf(:,2).*cos(prf(:,1)/180*pi) prf(:,2).*sin(prf(:,1)/180*pi) ];
                    prf(:,1:2)  = xy;                
                    prf_sub= perform_mesh_subdivision(prf', F', options)'; 
                    xy_sub = prf_sub(:,1:2);

                    theta = (atan2(xy_sub(:,2), xy_sub(:,1)))/pi*180;
                    if (contains(infn, 'rh'))
                        theta(theta<0) = theta(theta<0) + 360;
                    end

                    r = sqrt(xy_sub(:,1).^2 + xy_sub(:,2).^2);
                    prf_sub(:,1:2) = [theta  r];
                    Es.(ak) = prf_sub;
                else
                    Es.(ak) = perform_mesh_subdivision(field', F', options)';
                end
            end
        end
    end

    % slightly smooth use spherical harmonic mapping
    mVs = mean(Vs);
    surf.vertices = Vs -mVs;
    surf.faces = Fs;
    sphere.vertices = Es.Vertex_Vsphere/100;
    L = 100;
    sigma = 0.0015;
    surfsmooth = SPHARMsmooth2(surf,sphere, L, sigma);

    Vss = surfsmooth.vertices + mVs;
    dv = mean(sqrt(dot(Vss - Vs ,Vss - Vs,2)));

    % fast marching
    if contains(infn, 'lh')
        start_points = 103615;
    else
        start_points = 103359; % rh
    end
    [D,~,~] = perform_fast_marching_mesh(Vss, Fs, start_points);

    %
    %
    radius = 90;
    ind2del = find(D>radius);

    [Fcut, ~, vertex_father] = gf_remove_mesh_vertices(Fs, Vss, ind2del);
    V_cut = Vss(vertex_father,:);

    attributes = fieldnames(Es);
    for k=1:numel(attributes)
        ak  = attributes{k};
        if( isnumeric(Es.(ak)) )
            field = Es.(ak);
            if size(field,1) == size(Vss,1)
                Ecut.(ak) = field(vertex_father,:);
            end
        end
    end 

    % E.Vertex_atlaswang = E0.Vertex_atlaswang(vertex_father);

    roi =  find( Ecut.Vertex_atlashcp == 1 );
    roi_f = ismember(Fcut(:,1), roi) & ismember(Fcut(:,2), roi) & ismember(Fcut(:,3), roi);

    tic
    uv = disk_conformal_map(Fcut, V_cut, roi);
    absmu = abs(compute_bc(Fcut, uv, V_cut));
    meanmu = mean(absmu(roi_f));
    toc

    fprintf('mean mu =%f\n', meanmu)

    write_mfile(outfn,'Face',Fcut,...
       'Vertex %d %f %f %f {father=(%d) uv=(%f %f) rgb=(%f %f %f) Vwhite=(%f %f %f) Vinflate=(%f %f %f) Vpial=(%f %f %f) Vsphere=(%f %f %f) prf=(%f %f %f %f %f %f) atlaswang=(%d) atlashcp=(%d)}\n', ...
                [V_cut, vertex_father, uv, Ecut.Vertex_rgb, Ecut.Vertex_Vwhite, ...
                Ecut.Vertex_Vwhite, Ecut.Vertex_Vpial, Ecut.Vertex_Vsphere, ...
                Ecut.Vertex_prf Ecut.Vertex_atlaswang Ecut.Vertex_atlashcp]);        


    write_mfile(smoothfn,'Face',Fs,...
       'Vertex %d %f %f %f {rgb=(%f %f %f) Vwhite=(%f %f %f) Vinflate=(%f %f %f) Vpial=(%f %f %f) Vsphere=(%f %f %f) prf=(%f %f %f %f %f %f) atlaswang=(%d) atlashcp=(%d)}\n', ...
                [Vs, Es.Vertex_rgb, Es.Vertex_Vwhite, ...
                Es.Vertex_Vwhite, Es.Vertex_Vpial, Es.Vertex_Vsphere, ...
                Es.Vertex_prf Es.Vertex_atlaswang Es.Vertex_atlashcp]);
end

function run_step2(infn, outfn_sm, outfilename, outfn_usm, out_sv_figs, in_roi_confidence)
    compute_prf_smoothing(infn, outfn_sm, outfilename, outfn_usm, out_sv_figs, in_roi_confidence)
end

function run_step3(infn_s, infn_u, outfn_sm_filename, outfn_sv_figs, outfn_left_ecc, outfn_left_ang, outfn_right_ecc, outfn_right_ang, in_roi_confidence)
    if contains(infn_s,'lh')
        LoR = 0;
    else
        LoR = 1;
    end

    generate_bcm_results(infn_s, infn_u, outfn_sm_filename, outfn_sv_figs, in_roi_confidence, outfn_left_ecc, outfn_left_ang, outfn_right_ecc, outfn_right_ang)
end