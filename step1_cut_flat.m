function step1_cut_flat(batch)  % Run the code on only one single subject or batch of subjects
    if batch == '1'             % Input = '1' batch mode 
        step1_cut_flat_batch()  % Call the Cut and Flattening function for batch mode
    else                        % Input = '0' single selected subject
        step1_cut_flat_select() % Call the Cut and Flattening Function for one selected single subject
    end
end

function step1_cut_flat_batch()     % Cut and Flattening Function for batch mode
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
end

function step1_cut_flat_select()  % Cut and Flattening Function for selected single subject
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
end

%% Main Function for cut and flattening 
% Inputs :
% 1. Path info for the data which this function is going to to work on it
% 2. Path info for the output of this step
% 3. Path info for the smoothed data

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
