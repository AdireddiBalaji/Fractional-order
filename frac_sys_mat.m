function [A, B, c] = frac_sys_mat(x,N)
char_freq = 1; decs = 4;
    
    % Comments on parameter values:
    % N is the number of elements, and the size of the approximating matrices;
    % char_freq (characteristic frequency) and decs (width in decades) roughly
    % determine the frequency range over which the approximation is sought.
    % Large values of decs (larger than, say, 8) may cause trouble within finite
    % precision arithmetic in some analyses of some problems. 
    % Note that the role of these parameters is merely to define beta1 and beta2 below.

    yc = char_freq^x;
    beta1 = log10(yc) - 0.4 * decs;
    beta2 = beta1 + 0.6 * decs;
    d = logspace(beta1, beta2, N - 1);
    p = d.^2 ./ (1 + d.^2); 
    p = [0, p, 1];

    A = zeros(N); 
    B = A; 
    c = zeros(N, 1);
    
    for i = 1:N-1
        a = p(i); b = p(i+1); a0 = -a / (b - a);
        a1 = 1 / (b - a);
        
        A(i, i+1) = x * (a0 * (1 - a0) * mb(x - 1, -x - 1, a, b) ...
            + (1 - 2 * a0) * a1 * mb(x, -x - 1, a, b) ...
            - a1^2 * mb(x + 1, -x - 1, a, b));
        
        A(i+1, i+1) = A(i+1, i+1) + x * (a0^2 * mb(x - 1, -x - 1, a, b) ...
            + 2 * a0 * a1 * mb(x, -x - 1, a, b) ...
            + a1^2 * mb(x + 1, -x - 1, a, b));
        
        A(i, i) = A(i, i) + x * ((1 - a0)^2 * mb(x - 1, -x - 1, a, b) ...
            - 2 * (1 - a0) * a1 * mb(x, -x - 1, a, b) ...
            + a1^2 * mb(x + 1, -x - 1, a, b));

        B(i, i+1) = x * (a0 * (1 - a0) * mb(x, -x - 2, a, b) ...
            + (1 - 2 * a0) * a1 * mb(x + 1, -x - 2, a, b) ...
            - a1^2 * mb(x + 2, -x - 2, a, b));

        B(i+1, i) = B(i, i+1); 
        A(i+1, i) = A(i, i+1);

        B(i+1, i+1) = B(i+1, i+1) + x * (a0^2 * mb(x, -x - 2, a, b) ...
            + 2 * a0 * a1 * mb(x + 1, -x - 2, a, b) ...
            + a1^2 * mb(x + 2, -x - 2, a, b));

        B(i, i) = B(i, i) + x * ((1 - a0)^2 * mb(x, -x - 2, a, b) ...
            - 2 * (1 - a0) * a1 * mb(x + 1, -x - 2, a, b) ...
            + a1^2 * mb(x + 2, -x - 2, a, b));

        c(i+1) = c(i+1) + x * (a0 * mb(x - 1, -x - 1, a, b) ...
            + a1 * mb(x, -x - 1, a, b));

        c(i) = c(i) + x * ((1 - a0) * mb(x - 1, -x - 1, a, b) ...
            - a1 * mb(x, -x - 1, a, b));
    end

    a = p(N); 
    b = 1; 
    a0 = 1 / (b - a);
    
    A(N, N) = A(N, N) + a0^2 * x * mb(x - 1, -x + 1, a, b);
    B(N, N) = B(N, N) + a0^2 * x * mb(x, -x, a, b);
    c(N) = c(N) + a0 * x * mb(x - 1, -x, a, b);

    c = c / sqrt(gamma(1 - x) * gamma(1 + x));
end

function s = mb(m, n, a, b)
    if (n > -1) && (m > -1)
        bb = beta(1 + m, 1 + n);
        s = (betainc(b, 1 + m, 1 + n) - betainc(a, 1 + m, 1 + n)) * bb;
    elseif (n <= -1) && (m > 0)
        s = (a^m * (1 - a)^(n + 1) - b^m * (1 - b)^(n + 1)) / (n + 1) ...
            + m / (n + 1) * mb(m - 1, n + 1, a, b);
    else
        s = mb(m, n + 1, a, b) + mb(m + 1, n, a, b);
    end
end
