function [Jp, Jw] = computeJacobian(u,R,p,s)


    % [~, ns] = size(p);
    % 
    % p_skew = zeros(3,3,ns);
    % 
    % p_skew(1,2,:) = -p(3,:);
    % p_skew(2,1,:) = p(3,:);
    % p_skew(1,3,:) = p(2,:);
    % p_skew(3,1,:) = -p(2,:);
    % p_skew(2,3,:) = -p(1,:);
    % p_skew(3,2,:) = p(1,:);
    % 
    % p_skew = reshape(p_skew, 3, []);
    % 
    % % compute J
    % % J = zeros(3*ns, 3*ns);
    % % for i = 2:ns
    % %     J(3*i-2:3*i, 1:3*i-3) = p_skew(:, 1:3*i-3) - repmat(p_skew(:, 3*i-2:3*i), 1, i-1);
    % % end
    % % 
    % % for j = 1:ns
    % %     J(:, 3*j-2:3*j) = J(:, 3*j-2:3*j) * R(:,:,j);
    % % end
    % 
    % J = repmat(p_skew, ns, 1);
    % J = J + J';
    % 
    % for j = 1:ns
    %     J(1:3*j, 3*j-2:3*j) = 0;
    %     J(:, 3*j-2:3*j) = J(:, 3*j-2:3*j) * R(:,:,j);
    % end

    [~, ns] = size(p);

    pSkew = zeros(3*ns,3);

    pSkew(1:3:end, 2) = -p(3,:);
    pSkew(1:3:end, 3) =  p(2,:);
    pSkew(2:3:end, 1) =  p(3,:);
    pSkew(2:3:end, 3) = -p(1,:);
    pSkew(3:3:end, 1) = -p(2,:);
    pSkew(3:3:end, 2) =  p(1,:);

    Jp= zeros(3*ns);
    Jw = zeros(3*ns);

    for j = 1:ns

        if j == ns
            continue;
        end
        
        ds = s(j+1)-s(j);

        % [dexpu, C_halfv] = diffExp2(u(:,j), ds);
        % [dexpu_t, C_halfv_t] = diffExp2(u(:,j)*ds, ds);
        [dexpu, C_halfv] = diffExp2(u(:,j)*ds, ds);
        % C_halfv = 0;   % why????????????
        
        Jp(3*j+1:end, 3*(j-1)+1:3*j) = (repmat(pSkew(3*(j-1)+1:3*j,:) * R(:,:,j) * dexpu + R(:,:,j)*C_halfv, ns-j,1) - pSkew(3*j+1:end,:) * R(:,:,j) * dexpu) * ds;
        Jw(3*j+1:end, 3*(j-1)+1:3*j) = R(:,3*j+1:end)'*R(:,:,j)*dexpu * ds;
    end


end