function calculate_meanmu()
    dirinfo = dir('data_step1/*.m');
    meanmu_all = [];

    for i=1:length(dirinfo)
       fn = ['data_step1/' dirinfo(i).name]; 
       [F,V,E]=read_mfile(fn);

       roi =  find( E.Vertex_atlashcp == 1 );
       roi_f = ismember(F(:,1), roi) & ismember(F(:,2), roi) & ismember(F(:,3), roi);

       absmu = abs(compute_bc(F, E.Vertex_uv, V));
       meanmu = mean(absmu(roi_f));

       meanmu_all(i)= meanmu;
    end

    avg_mu = mean(meanmu_all);
    min_mu = min(meanmu_all);
    max_mu = max(meanmu_all);
    
    fprintf('avg mu =%f\n', avg_mu)
    fprintf('min mu =%f\n', min_mu)
    fprintf('max mu =%f\n', max_mu)
end
