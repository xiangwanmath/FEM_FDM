dimen = 100;

L = diag(ones(dimen-1, 1), -1);
U = diag(ones(dimen-1, 1), 1);
D = -4 * diag(ones(dimen, 1));
D_inv = -(1/4) * diag(ones(dimen, 1));


b = sin(linspace(0, 2*pi, dimen))';
x0 = (1:dimen)';

M =  -D_inv * (L+U);
N = D_inv;

num_iteration = 30;

M
N
x0

x = x0;
for i = 1:num_iteration
    x = M*x + N*b;
end

plot(x)