function [shape, results, ini_sol, res] = discreteCosserat_slide(tdcr, prev_results, loads, ini_guess)
    function ys = ode_1seg(s, y)
        R = reshape(y(4:12), 3, 3);
        u = y(13:15);
        n = y(16:18);
        v = [0;0;1];

        K = tdcr.K; % diag([tdcr.kb, tdcr.kb, tdcr.kt]);

        % formulate the diff equation
        le = zeros(3,1);  % DOUBLE CHECK!!!
        du = inv(K) * (- hat(u)*K*u - hat(v)*R'*n - R'*le);
        dn = - loads.f_body ;
        ys = [R * v;
            reshape(R * hat(u), 9, 1);
            du;
            dn];
    end


    function [b, results] = residual_tdcr(ini6ncc)
        % integrate iteratively

        nn = tdcr.n;
        mm = tdcr.m;

        ini6n = ini6ncc(1:6*nn);
        inicc = ini6ncc(6*nn+1:end);

        p0 = tdcr.base(1:3, 4);
        R0 = tdcr.base(1:3, 1:3);

        K = tdcr.K; % diag([tdcr.kb, tdcr.kb, tdcr.kt]);
        mu = tdcr.mu;

        pend_prev = p0;
        Rend_prev = R0;

        vz = [0; 0; -1];   % direction of actuation.
        D = [1, -1];
        e1 = [1, 1]';

        ys = {};
        pcs = zeros(3, nn+1);
        pcs(:,1) = p0;
        Rs = zeros(3,3,nn+1);
        Rs(:,:,1) = R0;

        s_list = [];
        shape_list = [];
        R_list = [];

        % get angle from vector v1, v2
        p2theta = @(v1, v2) pi - atan2(norm(cross(v1,v2)), dot(v1, v2));

        % cc 0<=x1 perp. x2 >=0  to boundary  b(x1, x2)=0
        cc2bvp = @(x1, x2) sqrt(x1.^2+x2.^2) - x1 - x2;

        for i = 1:nn

            ini6_i = ini6n(1+6*(i-1):6*i);

            ini18_i = [pend_prev; reshape(Rend_prev,9,1); reshape(ini6_i, 6, 1)];

            sRange_i = [tdcr.l*(i-1), tdcr.l*i];
            [s_sol, y_sol] = ode45(@ode_1seg, sRange_i, ini18_i);

            % update
            pend_prev = y_sol(end, 1:3)';
            Rend_prev = reshape(y_sol(end, 4:12),3,3);

            % store
            ys{end+1} = y_sol;
            pcs(:, i+1) = y_sol(end, 1:3)';
            Rs(:, :, i+1) = reshape(y_sol(end, 4:12), 3, 3);
            s_list = [s_list, s_sol'];
            shape_list = [shape_list, y_sol(:, 1:3)'];
            R_list = [R_list, y_sol(:, 4:12)'];

        end

        R_list = reshape(R_list, 3, 3, []);

        % get state variables
        P = cell(1, mm);
        theta = cell(1, mm);
        T = cell(1, mm);
        l = cell(1, mm);
        v = cell(1, mm);
        F = cell(1, mm);
        m = cell(1, mm);

        for j = 1:mm
            nj = tdcr.nj(j);

            % get P
            Pj = zeros(3, nj+1);
            for i = 1:nj+1
                rj = tdcr.r{j}(:, i);
                Pj(:, i) = pcs(:, i) + Rs(:,:,i)*rj;
            end
            P{j} = Pj;

            % get theta
            theta_j = zeros(1,nj);
            theta_j(1) = p2theta(Pj(:,2)-Pj(:,1), vz);
            for i = 2:nj
                theta_j(i) = p2theta(Pj(:,i+1)-Pj(:,i), Pj(:,i-1)-Pj(:,i));
            end
            theta{j} = theta_j;
            
            % get T
            Tj = zeros(1, nj+1);
            Tj(1) = loads.tension(j);
            Tj0_prev = prev_results.loads.tension(j);

            if Tj(1) > Tj0_prev
                % pull down
                for i = 1:nj
                    Tj(i+1) = Tj(i)*(exp(-mu(i)*theta_j(i)));
                end
                T{j} = Tj;
            else
                % release back
                for i = 1:nj
                    Tj(i+1) = Tj(i)*(exp(mu(i)*theta_j(i)));
                end
                T{j} = Tj;
            end

            % get F, m??
            Fj = zeros(3, nj+1);
            mj = zeros(3, nj+1);

            Fj(:,1) = zeros(3,1);  mj(:,1) = zeros(3,1);

            for i = 2:nj+1
                rj = tdcr.r{j}(:, i);
                if i == nj+1
                    % end of tendon.
                    Fj(:,i) = Tj(i)*Normalize(Pj(:,i-1)-Pj(:,i));
                else
                    Fj(:,i) = Tj(i)*Normalize(Pj(:,i-1)-Pj(:,i)) + Tj(i+1)*Normalize(Pj(:,i+1)-Pj(:,i));
                end
                mj(:,i) = hat(Rs(:,:,i)*rj)*Fj(:,i);
            end

            F{j} = Fj;
            m{j} = mj;

            % get l, v
            lj = zeros(1, nj+1);
            vj = zeros(1, nj);
            for i = 1:nj
                lj(i+1) = norm(Pj(:,i+1)-Pj(:,i));
            end
            lj(1) = tdcr.oriTendonLen(j) - sum(lj);   % tendon always tensioned.
            l{j} = lj;

            lj_prev = prev_results.l{j};
            for i = 1:nj
                vj(i) = sum(lj_prev(i+1:nj+1)) - sum(lj(i+1:nj+1));
            end
            v{j} = vj;

        end

        % extract the boundary
        b6n = zeros(6*nn, 1);   % boundary vector

        % extract contacts
        contacts = loads.contacts;
        nc = length(contacts);
        contact_idx = [contacts.diskID];
        contact_force = [contacts.F];
        contact_r = [contacts.r];

        % i = 1 ~ n-1
        for i = 1:nn
            % y_curr = ys{i};

            % find total force
            idx_tendonsNum = find(i <= tdcr.nj);

            F_tendon = zeros(3,1);
            m_tendon = zeros(3,1);
            for jj = idx_tendonsNum
                F_tendon = F_tendon + F{jj}(:,i+1);
                m_tendon = m_tendon + m{jj}(:,i+1);
            end

            % external loads
            F_ext = zeros(3,1);
            m_ext = zeros(3,1);
            
            for ic = find(contact_idx == i)
                F_ext = F_ext + contact_force(:,ic);
                m_ext = m_ext + hat(Rs(:,:,i+1)*contact_r(:,ic))*contact_force(:,ic);
            end

            F_total = F_tendon + F_ext;
            m_total = m_tendon + m_ext;

            % force/moment balance b(x)=0
            if i < nn
                y_after = ys{i+1};
                y_curr = ys{i};

                u_after = y_after(1, 13:15)';
                n_after = y_after(1, 16:18)';

                u_curr = y_curr(end, 13:15)';
                n_curr = y_curr(end, 16:18)';

                % b6n(1+6*(i-1):6+6*(i-1), :) = ...
                %     [u_curr - u_after  - inv(K)*Rs(:,:,i+1)'*(m_total);
                %     n_curr - n_after  - F_total];
                b6n(1+6*(i-1):6+6*(i-1), :) = ...
                    [K*(u_curr - u_after)  - Rs(:,:,i+1)'*(m_total);
                    n_curr - n_after  - F_total];

            else
                y_curr = ys{i};

                u_end = y_curr(end, 13:15)';
                n_end = y_curr(end, 16:18)';

                % b6n(1+6*(i-1):6+6*(i-1), :) = ...
                %     [u_end - inv(K)*Rs(:,:,i+1)'*(m_total);
                %     n_end  - F_total];
                b6n(1+6*(i-1):6+6*(i-1), :) = ...
                    [K*u_end - Rs(:,:,i+1)'*(m_total);
                    n_end  - F_total];
            end

        end

        b = b6n;

        % outputs
        % if nargout == 2
            results.loads = loads;
            results.pc = pcs;
            results.R_node = Rs;
            results.P = P;
            results.shape = shape_list;
            results.s = s_list;
            results.l = l;
            results.v = v;
            results.theta = theta;

            % find the bending angle
            R_bend = R0'*Rs(1:3, 1:3, end);

            results.bendAngle = pi/2 - atan2(R_bend(3,3), R_bend(1,3));

            ny = size(shape_list, 2);
            pp = reshape(shape_list, 3, 1, ny);
            base = repmat([0 0 0 1], [1,1,ny]);
            results.TT = cat(1, cat(2, R_list, pp), base);

        % end

    end

% shooting method
opt = optimoptions('fsolve', 'Display','none',...
    'Algorithm','levenberg-marquardt');
% opt = optimoptions('fsolve', 'Display','none',...
%     'Algorithm','levenberg-marquardt', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 500);
[ini_sol, res0] = fsolve(@residual_tdcr, ini_guess, opt);

% [ini_sol, res0] = fsolve(@residual_tdcr, zeros(tdcr.n*6, 1), opt);

if sum(res0.^2) >= 1e-6
    warning('big error')
end

[res, results] = residual_tdcr(ini_sol);
results.ini_sol = ini_sol;

shape = results.shape;
% g_array = results.g_array;  % to be modified.

end