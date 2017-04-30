% octave -q raytrace.m


% object constructors


function object = new_plane(point,normal,color)
	normal /= norm(normal);
	
	if normal(3) < 0
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


function polynomial3 = polynomial_algebraic()

	% x^4+y^4+z^4-1
	
	polynomial3(1,1,1) = -1;
	polynomial3(5,1,1) = 1;
	polynomial3(1,5,1) = 1;
	polynomial3(1,1,5) = 1;
end


% define materials


function color = color_plane(object,point,normal)
	%#
	color = abs(mod(point(3)*1,2)-1);
end

function color = color_sphere(object,point,normal)
	%#
	parallel.intensity = 0.5;
	parallel.direction = [0 0 0];

	color = 1+normal(3);
end

function color = color_algebraic(object,point,normal)
	%#
end


% define scene


clear objects;

objects{end+1} = new_plane([-1/2 -1/2 -1/2],[1/2 -1 -1/4],@color_plane);
objects{end+1} = new_sphere([0 0 0],0.9,@color_sphere);
objects{end+1} = new_sphere([1/2 0 -3/5],0.3,@color_sphere);
objects{end+1} = new_sphere([-1/2 -1/2 -1/2],0.5,@color_sphere);
%objects{end+1} = new_algebraic([0 0 0],@polynomial_algebraic,@color_algebraic);


% end of scene definition


% initialize variables


global rows = [-1 : 2/80 : 1];
global cols = [-1 : 2/80 : 1];
global appx = 10;

global objmap = cell(numel(rows),numel(cols));
global depths = ones(numel(rows),numel(cols)) * inf;
global pixels = zeros(numel(rows),numel(cols));


% step 1: determine intersections


function depth = depth_plane(object,coords)
	depth = inf;
	
	
	x = object.point - coords;
	n = object.normal;
%	v = [0 0 1];
	
	
	% solve linear equation
	%
	% 0 = <x - v*t,n>
	%
	% 0 = sum_i(x_i*n_i) - n_3*t
	%
	% t = sum_i(x_i*n_i)/n_3
	
	
	if n(3) != 0
		
		xn = x .* n;    % x_i*n_i
		s = sum(xn);    % sum_i(x_i*n_i)
		
		
		                % sum_i(x_i*n_i)/n_3
		depth = s / n(3);
	end
end

function depth = depth_sphere(object,coords)
	depth = inf;
	
	
	x = object.point - coords;
	r = object.radius;
%	v = [0 0 1];
	
	
	% solve quadratic equation
	%
	% 0 = |x - v*t| - r
	%
	% 0 =   sum_i((v_i)^2)  *t^2  } = t^2
	%     - sum_i(2*x_i*v_i)*t    } = b*t
	%     + sum_i((x_i)^2) - r^2  } = c
	%
	% t = -b/2 - sqrt(b^2/4 - c)
	
	
	x2 = x .* x;        % (x_i)^2
	r2 = r * r;         % r^2
	
	b = x(3);           % -sum_i(2*x_i*v_i) /-2
	c = sum(x2) - r2;   %  sum_i((x_i)^2) - r^2
	
	b2 = b * b;         % b^2/4
	d = b2 - c;         % b^2/4 - c
	
	
	if d >= 0
		                % -b/2 - sqrt(b^2/4 - c)
		depth = b - sqrt(d);
	end
end

function polynomial = distance(polynomial3,coords)
	polynomial(size(polynomial3,3)) = 0;
	
	
	p = polynomial3;
	x = coords;
%	v = [0 0 1];
	
	
	% distance polynomial for x + v*t
	%
	% p(x + v*t) =
	%      sum_ijk(p_ijk*(x_1 + v_1*t)^i
	%                   *(x_2 + v_2*t)^j
	%                   *(x_3 + v_3*t)^k)
	%
	% q_k = sum_ij(p_ijk*(x^1)^i*(x^2)^j)
	
	
	for i1 = 1:size(p,1)
		for i2 = 1:size(p,2)
			for i3 = 1:size(p,3)
				
				c = p(i1,i2,i3);
				c *= x(1)^(i1 - 1);
				c *= x(2)^(i2 - 1);
				
				polynomial(i3) += c;
				
			end
		end
	end
end

function derivative = derive(polynomial)
	
	
	p = polynomial;
	
	
	% derivative of polynomial in t
	%
	% q_(i - 1) = p_i*i
	
	
	for i1 = 2:numel(p)
		
		c = p(i1);
		c *= i1 - 1;
		
		derivative(i1 - 1) = c;
		
	end
end

function value = evaluate(polynomial,t)
	value = 0;
	
	
	p = polynomial;
	
	
	% value of polynomial in t
	%
	% p(t) = sum_i(p_i*t^i)
	
	
	for i1 = 1:numel(p)
		
		v = p(i1);
		v *= t^(i1 - 1);
		
		value += v;
		
	end
end

function depth = depth_algebraic(object,coords)
	global appx;
	depth = inf;
	
	
	p = distance(object.polynomial3,coords);
	q = derive(p);
	
	%object.polynomial3
	p
	q
	%exit;
	
	x = coords;
	
	for i = 1:appx
	end
	
	
	
	
	%#
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


% step 2: determine normals and colors


function normal = normal_plane(object,point)
	normal = object.normal;
end

function normal = normal_sphere(object,point)
	normal = point - object.point;
	normal /= norm(normal);
end

function normal = normal_algebraic(object,point)
	%#
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
			
			pixels(row,col) = object.color(object,point,normal);
			
		end
	end
end


% display image


close all
figure = imshow(pixels);
axis equal
axis off
waitfor(figure)


% end of code

