clear all;
close all;

%Conjugate Gradient Methods 
%Algorithm 5.2
n = 1000;
[Q,R] = qr(rand(n,n));
dd = (1000 - 10) / 999;
D = diag(10:dd:1000);
%dd_1 = 2 / 499;
%D = diag([9:dd_1:11,999:dd_1:1001]);
A = Q.' * D * Q;




bb = rand(n,1);


x = zeros(1000,1);
r_0 = A * x - bb;
p_0 = -r_0;
k = 0;

num = 1000;
err = zeros(num,1);
r = r_0;
err(1) = norm(r);
p = p_0;
while norm(r) > 0.00001
	a = (r.' * r) / (p.' * A * p);
	x_new = x + a * p;

	r_new = r + a * A * p;
	b_new = (r_new.'*r_new)/(r.'*r);
	p_new = -r_new + b_new * p;
	k = k + 1;

	x = x_new;
	r = r_new;
	b = b_new;
	p = p_new;
	err(k) = norm(r);
end 

%plot
figure;
hold on;
plot([1:1:k], err(1:k));
pause;
