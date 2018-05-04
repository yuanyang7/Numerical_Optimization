clear all;
close all;

%costFunc
function y = costFunc(x, n)
	y = 0;
	for i = 1:n-1
		y = y + 100 * (x(i)^2 - x(i+1))^2 + (x(i) - 1)^2 ;
	end
end
%costFunc = @(x) ((c*x(1) - 2)^4 + x(2)^2 * (c * x(1)-2)^2 + (x(2) + 1)^2);

%gradFunc
function g = gradFunc(x, n)
	g = zeros(n,1);
	g(1) = 200*(x(1)^2 - x(2))*2*x(1) + 2*(x(1) - 1);
	for i = 2:n-1
		g(i) = 200*(x(i)^2 - x(i+1))*2*x(i) + 2*(x(i) - 1) - 200*(x(i-1)^2 - x(i));
	end
end

%FR
n = 100;
x_0 = 2*rand(n,1);
f = zeros(1000,1);
theta = zeros(1000,1);
f(1) = costFunc(x_0, n);
g_0 = gradFunc(x_0, n);
p_0 = -g_0;
k = 0;

c = 10;
c1 = 1e-4;
c2 = 0.4;
a_max = 1e+6;

x = x_0;
g = g_0;
p = p_0;
while (norm(g) > 0.001) && (k < 1000)
	if mod(k,10) == 0
		norm(g)
	end
	a = StepLength(p,x,c1,c2,a_max,n);
	%a = lineSearch(x, p, c1, c2, a_max, c,n);
	x_new = x + a * p;

	g_new = gradFunc(x_new, n);
	b_new = (g_new.' * (g_new - g)) / norm(g)^2;
	p_new = -g_new + b_new * p;
	k = k + 1;

	x = x_new;
	p = p_new;
	g = g_new;
	theta(k+1) = -g'*p /(norm(g) * norm(p));
	f(k + 1) = costFunc(x, n);

end 

%plot
figure(1);
hold on;
plot([1:1:k], f(1:k));
figure(2);
hold on;
plot([1:1:k], theta(1:k));
pause;

