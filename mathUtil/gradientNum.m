function grad = gradientNum(f, x)
    % Calculate the gradient of a vector-valued function f(x) with respect to x
    % f: function handle representing the n-output function f(x)
    % x: m-dimensional vector at which to evaluate the gradient
    % grad: n x m matrix of gradients, where grad(i,j) is the partial
    % derivative of the i-th output of f with respect to the j-th component of x
    
    % Small perturbation for numerical differentiation
    epsilon = 1e-6;  
    
    % Evaluate the function at the current x
    fx = f(x);
    
    % Initialize the gradient matrix
    n = length(fx);   % Number of outputs
    m = length(x);    % Number of inputs
    grad = zeros(n, m);
    
    % Loop over each element of x
    for j = 1:m
        % Create a perturbed version of x
        x_perturbed = x;
        x_perturbed(j) = x_perturbed(j) + epsilon;
        
        % Evaluate the function at the perturbed x
        fx_perturbed = f(x_perturbed);
        
        % Calculate the numerical partial derivative for each output
        for i = 1:n
            grad(i, j) = (fx_perturbed(i) - fx(i)) / epsilon;
        end
    end
end