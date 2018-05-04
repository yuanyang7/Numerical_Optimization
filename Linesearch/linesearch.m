clear all;
close all;

%algorithm 3.5
function a_cur = lineSearch(x, p, c1, c2, a_max, c)
	a_0 = 0;
	a_cur = a_max / 2 ;
	a_prev = a_0;
	i = 0;
	phi0 = costFunc(x,c);
	phi0_d = gradFunc(x,c)'*p;
	phi_prev = phi0;
	while i < 1000
		i = i + 1;
		phi_cur = costFunc(x + a_cur * p,c); %todo cost

		if phi_cur > phi0 + c1 * a_cur * phi0_d || (phi_cur >= phi_prev && i > 1)
			a_cur = zoom(a_prev, a_cur, phi0, phi0_d, x, p, c1, c2, c);
			break;
		end
		phi_d_cur = gradFunc(x + a_cur * p,c)'*p;
		if abs(phi_d_cur) <= -c2*phi0_d
			break;
		end
		if phi_d_cur >= 0
			a_cur = zoom(a_cur, a_prev, phi0, phi0_d, x, p, c1, c2, c);
			break;
		end
		a_prev = a_cur;
		a_cur = (a_cur + a_max) / 2;
		phi_prev = phi_cur;
	end
end

%algorithm 3.6
function a_cur = zoom(a_low, a_high, phi0, phi0_d, x, p, c1, c2, c)
	%az_sol = (a_low + a_high) / 2;
	phi_low = costFunc(x + a_low * p,c);
	j = 0;
	a_cur = (a_low + a_high) / 2;
	while j < 1000
		j = j + 1;
		a_cur = (a_low + a_high) / 2;
		phi_cur = costFunc(x + a_cur * p,c);
		if (phi_cur > phi0 + c1 * a_cur * phi0_d) || (phi_cur >= phi_low)
			a_high = a_cur;
		else 
			phi_d_cur = gradFunc(x + a_cur * p,c)'*p; %todo
			if abs(phi_d_cur) <= -c2*phi0_d	
				break;
			end
			if 	phi_d_cur * (a_high - a_low) >= 0
				a_high = a_low;
			end
			a_low = a_cur;
			phi_low = costFunc(x + a_low * p,c);
		end
	end
end


%costFunc
function y = costFunc(x,c)
	y = ((c*x(1) - 2)^4 + x(2)^2 * (c * x(1)-2)^2 + (x(2) + 1)^2);
end
%costFunc = @(x) ((c*x(1) - 2)^4 + x(2)^2 * (c * x(1)-2)^2 + (x(2) + 1)^2);

%gradFunc
function g = gradFunc(x,c)
	g = ([(c^3 * x(1)^3 - 6 * c^2 * x(1)^2 + 12 * c * x(1) + 2 * c^2 * x(1) * x(2)^2 - 4 * c * x(2)^2 - 8 );( 2 * c^2 * x(1)^2 * x(2) - 8 * c * x(1) * x(2) + 10 * x(2) + 2 )]);
end
%gradFunc = @(x) ([(c^3 * x(1)^3 - 6 * c^2 * x(1)^2 + 12 * c * x(1) + 2 * c^2 * x(1) * x(2)^2 - 4 * c * x(2)^2 );( 2 * c^2 * x(1)^2 * x(2) - 8 * c * x(1) * x(2) + 10 * x(2) + 2 )]);

%steepest descent methods
c = 10;
c1 = 0.1;
c2 = 0.8;
num = 50;
a_max = 1;
x = zeros(num, 2);
err = zeros(num,1);
x_opt = [2/c;-1]; %todo
x(1,:) = [2;2]; %todo
err(1) = norm(x(1,:) - x_opt);
p = - gradFunc(x(1,:),c);
i = 1;
while i < num && norm(p) > 0.0001%todo
	if mod(i,100) == 0
		i
	end
	a = lineSearch(x(i,:), p, c1, c2, a_max, c);
	x(i+1,:) = x(i,:) + a * p';
	err(i+1) = norm(x(i+1,:) - x_opt);
	p = - gradFunc(x(i+1,:),c);
	i = i + 1;
end

%plot
figure;
hold on;
plot([1:1:num], err, '.');
pause;
