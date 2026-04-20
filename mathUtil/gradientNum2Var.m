function [grad_x, grad_z] = gradientNum2Var(f, x, z)
    % Calculate the gradient of a vector-valued function f(x, z) with respect to x and z
    % f: function handle representing the n-output function f(x, z)
    % x: m-dimensional vector at which to evaluate the gradient with respect to x
    % z: k-dimensional vector at which to evaluate the gradient with respect to z
    % grad_x: n x m matrix of gradients, where grad_x(i,j) is the partial derivative of 
    %         the i-th output of f with respect to the j-th component of x
    % grad_z: n x k matrix of gradients, where grad_z(i,j) is the partial derivative of 
    %         the i-th output of f with respect to the j-th component of z
    
    % Small perturbation for numerical differentiation
    epsilon = 1e-6;  
    
    % Evaluate the function at the current (x, z)
    fxz = f(x, z);
    
    % Initialize the gradient matrices
    n = length(fxz);   % Number of outputs
    m = length(x);     % Number of dimensions in x
    k = length(z);     % Number of dimensions in z
    
    grad_x = zeros(n, m);
    grad_z = zeros(n, k);
    
    % Compute gradient with respect to x
    for j = 1:m
        % Create a perturbed version of x
        x_perturbed = x;
        x_perturbed(j) = x_perturbed(j) + epsilon;
        
        % Evaluate the function at the perturbed (x, z)
        fxz_perturbed = f(x_perturbed, z);
        
        % Calculate the numerical partial derivative for each output
        grad_x(:, j) = (fxz_perturbed - fxz) / epsilon;

    end
    
    % Compute gradient with respect to z
    for j = 1:k
        % Create a perturbed version of z
        z_perturbed = z;
        z_perturbed(j) = z_perturbed(j) + epsilon;
        
        % Evaluate the function at the perturbed (x, z)
        fxz_perturbed = f(x, z_perturbed);
        
        % Calculate the numerical partial derivative for each output
        grad_z(:, j) = (fxz_perturbed - fxz) / epsilon;
    end
end