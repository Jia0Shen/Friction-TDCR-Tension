function [w, x] = LCPScaleSolve(M, q)

    % %% Solve scaled LCP
    sv = svd(M);
    scale = max(sv);

    dimM = size(M,1);

    Dw = diag([ones(1, dimM-nc_f), ones(1,nc_f)*scale]);
    Dx = diag([ones(1, dimM-nc_f)/scale, ones(1,nc_f)]);

    M2 = Dw * M * Dx;
    g2 = Dw * g;

    % [w2, x2, retcode2] = LCPSolve(M2, g2, 1e-8, dimM^2);
    % test_w2 = M2*x2+g;
    % retcode2  % 1 success; 2 max termination
    solver = 'matlab';
    if strcmp(solver, 'pythonLCP')
        M2_np= py.numpy.array(M2);
        g2_np= py.numpy.array(g2);
        MaxIter = dimM*10;
        sol = lcp.lemkelcp(M2_np, g2_np, py.int(MaxIter));
        disp(sol{3})
        x2 = double(sol{1})';
    elseif strcmp(solver, 'matlab')
        [w2, x2, retcode2] = LCPSolve(M2, g2, 1e-8, dimM^2);
        test_w2 = M2*x2+g2;
        % retcode2  % 1 success; 2 max termination
    end
    % disp('relative diff of w2/w2_ is ' + string(norm(w2-test_w2)/norm(w2)));
    % [w2, x2] = LCPSolve(M2, g2, 1e-16);
    % [w2, x2] = LCP_jh(M2, g2);

    x = Dx * x2;
end