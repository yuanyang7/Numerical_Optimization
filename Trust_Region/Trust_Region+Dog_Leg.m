clear all;
close all;

MAXITER = 1000;

n = 10;
[Q,R] = qr(rand(n,n));
lambda = linspace(1,100,n);
Q = Q'*diag(lambda)*Q;

syms x1;
X = x1;
for i = 2:n
    syms(['x',num2str(i)]);
    X = [X,['x',num2str(i)]];
end

syms p1;
P = p1;
for i = 2:n
    syms(['p',num2str(i)]);
    P = [P,['p',num2str(i)]];
end

syms p;
ff = log(1+X*Q*X');
gg = jacobian(ff,X);
bb = hessian(ff,X);
k = 0;
x_k = ones(n,1);
delta_k = 10;
error = zeros(MAXITER,1);
alpha = 1;

while k < MAXITER 
    f = eval(subs(ff,X, x_k'));
    g = eval(subs(gg,X, x_k'));
    b = eval(subs(bb,X, x_k'));
    g = g';
    gbg = g' * b * g;
    
    %alpha
    if(gbg <= 0)
    	alpha = delta_k / norm(g);
    else
    	alpha = min(delta_k / norm(g), norm(g)^2/ gbg );
    end

    %p
	p_k = -alpha * g;
    %trustregion	
 	m1 = g' * zeros(n,1) + 0.5 * zeros(n,1)' * b * zeros(n,1);
    m2 = g' * p_k + 0.5 * p_k' * b * p_k;
    rho_k = eval((subs(ff,X, x_k') - subs(ff,X,(x_k + p_k)')) / (m1 - m2));
    if rho_k < 0.25
		delta_k = 0.25 * delta_k;
    else
        if rho_k > 0.75 && norm(p_k) == delta_k
		delta_k = min(2*delta_k, 10000);
        end 
    end
    k = k+1;
    if mod(k,100) == 0
        k
    end
	error(k) = f;

	if rho_k > 0
		x_k = x_k + p_k;
	end

end

%plot
figure;
hold on;
plot((1:1:k), log2(error(1:k)));
pause;
