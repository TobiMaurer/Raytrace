% octave -q raytrace.m


% object constructors


function object = new_plane(point,normal,color)
	normal /= norm(normal);
	
	if normal(3) > 0
		normal *= -1;
	end
	
	object.type = "plane";
	object.point = point;
	object.normal = normal;
	object.color = color;
end

function object = new_sphere(point,radius,color)
	object.type = "sphere";
	object.point = point;
	object.radius = radius;
	object.color = color;
end

function object = new_algebraic(point,polynomial3,color)
	object.type = "algebraic";
	object.point = point;
	object.polynomial3 = polynomial3();
	object.color = color;
end


% define surfaces


function polynomial3 = polynomial_cube()

	% x^4+y^4+z^4-1/3
	
	polynomial3(5,1,1) = 1;
	polynomial3(1,5,1) = 1;
	polynomial3(1,1,5) = 1;
	polynomial3(1,1,1) = -1/3;
end

function polynomial3 = polynomial_cube2()

	% x^4+(y*cos(phi)+z*sin(phi))^4+(z*cos(phi)-y*sin(phi))^4-1/3
	
	phi = pi/5;
	a = cos(phi);
	b = sin(phi);
	
	%  x^4
	% +(a^4+b^4)*y^4
	% +(a^4+b^4)*z^4
	% +(4*a^3*b-4*a*b^3)*y^3*z
	% +(12*a^2*b^2)*y^2*z^2
	% +(4*a*b^3-4*a^3*b)*y*z^3
	% -1/3
	
	polynomial3(5,1,1) = 1;
	polynomial3(1,5,1) = a^4+b^4;
	polynomial3(1,1,5) = a^4+b^4;
	polynomial3(1,4,2) = 4*a^3*b-4*a*b^3;
	polynomial3(1,3,3) = 12*a^2*b^2;
	polynomial3(1,2,4) = 4*a*b^3-4*a^3*b;
	polynomial3(1,1,1) = -1/3;
end

function polynomial3 = polynomial_torus(R,r)

	%  x^4
	% +y^4
	% +z^4
	% -2*(R^2+r^2)*x^2
	% -2*(R^2+r^2)*y^2
	% +2*(R^2-r^2)*z^2
	% +2*x^2*y^2
	% +2*x^2*z^2
	% +2*y^2*z^2
	% +R^4+r^4-2*R^2*r^2
	
	polynomial3(5,1,1) = +1;
	polynomial3(1,5,1) = +1;
	polynomial3(1,1,5) = +1;
	polynomial3(3,1,1) = -2*(R^2+r^2);
	polynomial3(1,3,1) = -2*(R^2+r^2);
	polynomial3(1,1,3) = +2*(R^2-r^2);
	polynomial3(3,3,1) = +2;
	polynomial3(3,1,3) = +2;
	polynomial3(1,3,3) = +2;
	polynomial3(1,1,1) = +R^4+r^4-2*R^2*r^2;
end


% define materials


function color = color_pattern(object,point,normal)
	color = 0;
	
	point -= object.point;
	base1 = cross(normal,[0 0 1]);
	base2 = cross(base1,normal);
	base1 /= norm(base1);
	base2 /= norm(base2);
	
	coord1 = dot(point,base1);
	coord2 = dot(point,base2);
	
	color += abs(mod(coord1,1) - 0.5);
	color += abs(mod(coord2,1) - 0.5);
end

function color = color_light(object,point,normal)
	color = 0;
	
	
	ambient.intensity = 0.1;
	
	parallel.intensity = 0.8;
	parallel.direction = [-1 -1 -1];
	
	specular.intensity = 0.3;
	specular.direction = [-2 -1 -3];
	
	
	% ambient light
	intensity = ambient.intensity;
	color += max(0,intensity);
	
	% parallel light
	direction = parallel.direction;
	direction /= norm(direction);
	intensity = parallel.intensity;
	intensity *= dot(direction,normal);
	color += max(0,intensity);
	
	% specular light
	direction = specular.direction;
	direction /= norm(direction);
	intensity = specular.intensity;
	intensity *= dot(direction,normal)^30;
	color += max(0,intensity);
end

function color = color_grid(object,point,normal)
	color = 0;
	
	%point *= [1 0 0;0 cos(pi/12) -sin(pi/12);0 sin(pi/12) cos(pi/12)];
	%point *= [cos(pi/12) 0 sin(pi/12);0 1 0;-sin(pi/12) 0 cos(pi/12)];
	
	for i = 1:3
		intensity = abs(mod(point(i) * 5,2) - 1) ^ 10;
		color = max(color,intensity);
	end
end


% define scene


clear objects;

%objects{end+1} = new_plane([0 1/2 0],[1 3 2],@color_pattern);
%objects{end+1} = new_sphere([0 0 0],0.9,@color_light);
%objects{end+1} = new_sphere([1/2 0 -3/5],0.3,@color_light);
%objects{end+1} = new_sphere([-1/2 -1/2 -1/2],0.5,@color_light);
objects{end+1} = new_algebraic([0 0 0],@polynomial_torus(0.6,0.2),@color_light);
%objects{end+1} = new_algebraic([0 0 0],@polynomial_cube2(),@color_light);
%objects{end+1} = new_algebraic([1/2 1/2 -1/4],@polynomial_cube2(),@color_light);


% end of scene definition


% predefined scenes


if exist("build") == 1
	clear objects;
	
	source("res.m");
	
	if build == 1
		objects{end+1} = new_algebraic([0 0 0],@polynomial_torus(0.6,0.2),@color_light);
	elseif build == 2
		objects{end+1} = new_algebraic([0 0 0],@polynomial_torus(0.6,0.2),@color_light);
	end
end


% initialize variables


global rows = [-1 : 2/40 : 1]; % projective surface
global cols = [-1 : 2/40 : 1];

global objmap = cell(numel(rows),numel(cols));
global depths = ones(numel(rows),numel(cols)) * inf;
global pixels = zeros(numel(rows),numel(cols));


% step 1: determine intersections


function depth = depth_plane(object,coords)
	depth = inf;
	
	
	x = coords - object.point;
	n = object.normal;
%	v = [0 0 1];
	
	
	% solve linear equation
	%
	% 0 = <x + v*t,n>
	%
	% 0 =  sum_i(x_i*n_i) + n_3*t
	%
	% t = -sum_i(x_i*n_i)/n_3
	
	
	if n(3) != 0
		xn = x .* n;    % x_i*n_i
		s = sum(xn);    % sum_i(x_i*n_i)
		
		                % -sum_i(x_i*n_i)/n_3
		depth = -s / n(3);
	end
end

function depth = depth_sphere(object,coords)
	depth = inf;
	
	
	x = coords - object.point;
	r = object.radius;
%	v = [0 0 1];
	
	
	% solve quadratic equation
	%
	% 0 = ||x + v*t|| - r
	%
	% 0 = t^2 + r^2 - sum_i((x_i)^2)
	%
	% t = -sqrt(r^2 - sum_i((x_i)^2))
	
	
	x2 = x([1,2]);
	x2 = x2 .* x2;      % (x_i)^2
	r2 = r * r;         % r^2
	
	c = r2 - sum(x2);   % r^2 - sum_i((x_i)^2)
	
	
	if c >= 0
		                % -sqrt(c) - x_3
		depth = -sqrt(c) - x(3);
	end
end

function polynomial = partial(polynomial3,coords,sel)
	polynomial = zeros(size(polynomial3,3),1);
	
	
	p = permute(polynomial3,sel);
	x = coords;
%	v = [0 0 1];
	
	
	% partial polynomial of x_3 at point x
	% for distance and partial derivation
	%
	% q_k = sum_ij(p_ijk*(x^1)^i*(x^2)^j)
	
	
	for i1 = 1:size(p,1)
		
		%exp((i1 - 1) * log(x(1)));
		x1 = x(1)^(i1 - 1);
		
		for i2 = 1:size(p,2)
			
			%exp((i2 - 1) * log(x(2)));
			x2 = x(2)^(i2 - 1);
			
			for i3 = 1:size(p,3)
				
				c = p(i1,i2,i3);
				
				c *= x1;
				c *= x2;
				
				polynomial(i3) += c;
				
			end
		end
	end
end

function derivative = derive(polynomial)
	derivative = zeros(numel(polynomial) - 1,1);
	
	
	p = polynomial;
	
	
	% derivative of polynomial
	%
	% q_(i - 1) = p_i*i
	
	
	for i1 = 2:numel(p)
		
		c = p(i1);
		c *= i1 - 1;
		
		derivative(i1 - 1) = c;
		
	end
end

function gradient = derive3(polynomial3,point)
	gradient = zeros(3,1);
	
	
	p = polynomial3;
	x = point;
	
	
	% gradient of polynomial in three
	% variables at point x
	%
	% g = [d_1(p,x) d_2(p,x) d_3(p,x)]
	
	
	q1 = partial(p,[x(2) x(3)],[2 3 1]);
	q2 = partial(p,[x(3) x(1)],[3 1 2]);
	q3 = partial(p,[x(1) x(2)],[1 2 3]);
	
	r1 = derive(q1);
	r2 = derive(q2);
	r3 = derive(q3);
	
	gradient(1) = evaluate(r1,x(1));
	gradient(2) = evaluate(r2,x(2));
	gradient(3) = evaluate(r3,x(3));
end

function value = evaluate(polynomial,t)
	value = 0;
	
	
	p = polynomial;
	
	
	% value of polynomial in t
	%
	% p(t) = sum_i(p_i*t^i)
	
	
	% Horner-Schema
	%
	% sum_i(p_i*x^i) = h_0
	% h_(n + 1) = 0
	% h_(i - 1) = h_i*x + p_(i - 1)
	
	
	for i1 = numel(p):-1:1
		
		value *= t;
		value += p(i1);
		
	end
end

global stat = zeros(20,1);

function depth = depth_algebraic(object,coords)
	depth = inf;
	
	
	x = coords - object.point;
	p = partial(object.polynomial3,x,[1 2 3]);
	q = derive(p);
	
	
	% solve polynomial equation
	% (newton's method)
	%
	%
	% start at t = -1, repeat:
	%
	% determine tangent at t,
	% set maximum slope to -1,
	% intersect with x-axis.
	%
	% determine maximum slope in interval,
	% intersect with x-axis.
	%
	%
	% abort on abs(p(t))<eps or after 20
	% iterations, whichever occurs first
	
	
	value = inf;
	t1 = -1;
	eps = 1e-2;
	
	for i = 1:20
		
		value = evaluate(p,t1);
		slope = evaluate(q,t1);
		slope = min(-1,slope);
		
		if abs(value) < eps
			break;
		end
		
		t2 = t1;
		t2 += value / slope;
		
		slope = evaluate(abs(q),max(abs([t1 t2])));
		
		t1 += value / slope;
		
	end
	
	
	global stat;
	stat(i) += 1;
	
	
	if abs(value) < eps
		depth = t1 - x(3);
	end
end

for row = 1:numel(rows)
	for col = 1:numel(cols)
		for object = objects
			
			object = object{:};
			coords = [cols(col) rows(row) 0];
			depth = inf;
			
			if strcmp(object.type,"plane")
				depth = depth_plane(object,coords);
			elseif strcmp(object.type,"sphere")
				depth = depth_sphere(object,coords);
			elseif strcmp(object.type,"algebraic")
				depth = depth_algebraic(object,coords);
			end
			
			if depth < depths(row,col)
				objmap{row,col} = object;
				depths(row,col) = depth;
			end
			
		end
	end
end

stat


% step 2: determine normals and colors


function normal = normal_plane(object,point)
	normal = object.normal;
end

function normal = normal_sphere(object,point)
	normal = point - object.point;
	normal /= norm(normal);
end

function normal = normal_algebraic(object,point)
	point -= object.point;
	normal = derive3(object.polynomial3,point);
	normal /= norm(normal);
end

for row = 1:numel(rows)
	for col = 1:numel(cols)
		
		object = objmap{row,col};
		depth = depths(row,col);
		
		if depth != inf
			
			point = [cols(col) rows(row) depth];
			
			if strcmp(object.type,"plane")
				normal = normal_plane(object,point);
			elseif strcmp(object.type,"sphere")
				normal = normal_sphere(object,point);
			elseif strcmp(object.type,"algebraic")
				normal = normal_algebraic(object,point);
			end
			
			pixels(row,col) = ...
				object.color(object,point,normal);
			
		end
	end
end


% display image


close all
figure = imshow(pixels);
axis equal
axis off

if exist("build") == 1
	print(figure,sprintf("scene%d.jpg",build),"-djpg")
else
	waitfor(figure)
end


% end of code

