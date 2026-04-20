function [traj_state] = getTraj2disk3DSlide(tdcr, T0_list)

% generate a traj given T0_list input.

ini_guess = zeros(6*tdcr.n, 1);
ini_sol = ini_guess;

loads.f_body = zeros(3,1);
loads.tension = 0;

contact1 = struct('diskID', tdcr.n, 'F', zeros(3,1), 'r', tdcr.r(:,end, 1));
loads.contacts = [contact1];

prev_results.l = {[0, tdcr.l, tdcr.l]};
prev_results.loads = loads;

loads_prev = loads;

traj_state = {};

n_steps = size(T0_list, 2);

for iter = 1:n_steps

    T0_All_i = T0_list(:,iter);

    loads.tension = T0_All_i;

    [~, results, ini_sol, res] = discreteCosserat_slide(tdcr, prev_results, loads, ini_sol);

    if norm(res) > 1e-6
        % not converge
        results = prev_results;
        results.loads = loads;
    end

    prev_results = results;
    traj_state{end+1} = results;

    if round(10*(iter+1)/n_steps) > round(10*iter/n_steps)

        fprintf('finish iteraction %d/%d; res is %.8f \n', iter, n_steps, norm(res));
    end
    
end


end