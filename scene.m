% do not run this file directly


% define surfaces


function polynomial3 = polynomial_cube()

	% x^4+y^4+z^4-1/16
	
	polynomial3(5,1,1) = 1;
	polynomial3(1,5,1) = 1;
	polynomial3(1,1,5) = 1;
	polynomial3(1,1,1) = -1/16;
end

function polynomial3 = polynomial_cube2()

	% x^4+(y*cos(phi)+z*sin(phi))^4+(z*cos(phi)-y*sin(phi))^4-1/16
	
	phi = pi/5;
	a = cos(phi);
	b = sin(phi);
	
	% +x^4
	% +(a^4+b^4)*y^4
	% +(a^4+b^4)*z^4
	% +(4*a^3*b-4*a*b^3)*y^3*z
	% +(12*a^2*b^2)*y^2*z^2
	% +(4*a*b^3-4*a^3*b)*y*z^3
	% -1/16
	
	polynomial3(5,1,1) = 1;
	polynomial3(1,5,1) = a^4+b^4;
	polynomial3(1,1,5) = a^4+b^4;
	polynomial3(1,4,2) = 4*a^3*b-4*a*b^3;
	polynomial3(1,3,3) = 12*a^2*b^2;
	polynomial3(1,2,4) = 4*a*b^3-4*a^3*b;
	polynomial3(1,1,1) = -1/16;
end

function polynomial3 = polynomial_torus(R,r)

	% +x^4
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

function polynomial3 = polynomial_torus2(R,r,phi)
	
	a = cos(phi);
	b = sin(phi);

	% +x^4
	% +(a^4+b^4+2*a^2*b^2)*y^4
	% +(a^4+b^4+2*a^2*b^2)*z^4
	% -2*(R^2+r^2)*x^2
	% -2*((a^2-b^2)*R^2+(a^2+b^2)*r^2)*y^2
	% +2*((a^2-b^2)*R^2-(a^2+b^2)*r^2)*z^2
	% +2*(a^2+b^2)*x^2*y^2
	% +2*(a^2+b^2)*x^2*z^2
	% +2*(a^4+2*a^2*b^2+b^4)*y^2*z^2
	% -8*a*b*R^2*y*z
	% +R^4+r^4-2*R^2*r^2
	
	polynomial3(5,1,1) = +1;
	polynomial3(1,5,1) = +(a^4+b^4+2*a^2*b^2);
	polynomial3(1,1,5) = +(a^4+b^4+2*a^2*b^2);
	polynomial3(3,1,1) = -2*(R^2+r^2);
	polynomial3(1,3,1) = -2*((a^2-b^2)*R^2+(a^2+b^2)*r^2);
	polynomial3(1,1,3) = +2*((a^2-b^2)*R^2-(a^2+b^2)*r^2);
	polynomial3(3,3,1) = +2*(a^2+b^2);
	polynomial3(3,1,3) = +2*(a^2+b^2);
	polynomial3(1,3,3) = +2*(a^4+2*a^2*b^2+b^4);
	polynomial3(1,2,2) = -8*a*b*R^2;
	polynomial3(1,1,1) = +R^4+r^4-2*R^2*r^2;
end

function polynomial3 = polynomial_flower()
	
	a = (1 + sqrt(5)) / 2;
	
	% +(2*a+1)*x^6
	% +(2*a+1)*y^6
	% +(2*a+1)*z^6
	% -(6*a+3)*x^4
	% -(6*a+3)*y^4
	% -(6*a+3)*z^4
	% -(4*a^2-6*a-3)*x^4*y^2
	% +(4*a^4+6*a+3)*x^4*z^2
	% +(4*a^4+6*a+3)*x^2*y^4
	% -(4*a^2-6*a-3)*y^4*z^2
	% -(4*a^2-6*a-3)*x^2*z^4
	% +(4*a^4+6*a+3)*y^2*z^4
	% +(6*a+3)*x^2
	% +(6*a+3)*y^2
	% +(6*a+3)*z^2
	% -(12*a+6)*x^2*y^2
	% -(12*a+6)*x^2*z^2
	% -(12*a+6)*y^2*z^2
	% -(4*a^6-12*a-10)*x^2*y^2*z^2
	% -2*a-1
	
	polynomial3(7,1,1) = +(2*a+1);
	polynomial3(1,7,1) = +(2*a+1);
	polynomial3(1,1,7) = +(2*a+1);
	polynomial3(5,1,1) = -(6*a+3);
	polynomial3(1,5,1) = -(6*a+3);
	polynomial3(1,1,5) = -(6*a+3);
	polynomial3(5,3,1) = -(4*a^2-6*a-3);
	polynomial3(5,1,3) = +(4*a^4+6*a+3);
	polynomial3(3,5,1) = +(4*a^4+6*a+3);
	polynomial3(1,5,3) = -(4*a^2-6*a-3);
	polynomial3(3,1,5) = -(4*a^2-6*a-3);
	polynomial3(1,3,5) = +(4*a^4+6*a+3);
	polynomial3(3,1,1) = +(6*a+3);
	polynomial3(1,3,1) = +(6*a+3);
	polynomial3(1,1,3) = +(6*a+3);
	polynomial3(3,3,1) = -(12*a+6);
	polynomial3(3,1,3) = -(12*a+6);
	polynomial3(1,3,3) = -(12*a+6);
	polynomial3(3,3,3) = -(4*a^6-12*a-10);
	polynomial3(1,1,1) = -2*a-1;
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

%objects{end+1} = new_plane([0.4 0.1 0],[1 3 2],@color_pattern);
%objects{end+1} = new_algebraic([0 0 0],@polynomial_torus(0.6,0.2),@color_light);
objects{end+1} = new_algebraic([0 0 0.5],@polynomial_torus2(0.6,0.2,-2*pi/5),@color_light);

%objects{end+1} = new_algebraic([0 0 0],@polynomial_cube2(),@color_light);

%objects{end+1} = new_algebraic([0 0 0],@polynomial_flower(),@color_light);
%objects{end+1} = new_sphere([0 0 1],0.8,@color_light);

%{
objects{end+1} = new_sphere([-0.05 -0.65 -0.25],0.1,@color_light);
objects{end+1} = new_sphere([-0.35 -0.65 -0.05],0.1,@color_light);
objects{end+1} = new_sphere([-0.1 -0.4 0],0.35,@color_light);
objects{end+1} = new_sphere([-0.1 0.1 0],0.3,@color_light);
objects{end+1} = new_sphere([0 0.5 0.1],0.25,@color_light);
objects{end+1} = new_sphere([0.2 0.7 0.2],0.15,@color_light);
%}


% end of scene definition


% predefined scenes


if exist("build") == 1
	clear objects;
	
	if build == 1
		objects{end+1} = new_algebraic([0 0 0],@polynomial_torus(0.6,0.2),@color_light);
	elseif build == 2
		objects{end+1} = new_algebraic([0 0 0.5],@polynomial_torus2(0.6,0.2,-pi/5),@color_light);
	elseif build == 3
		objects{end+1} = new_algebraic([0 0 0.5],@polynomial_torus2(0.6,0.2,-2*pi/5),@color_light);
	elseif build == 4
		objects{end+1} = new_algebraic([0 0 0],@polynomial_cube2(),@color_light);
	elseif build == 5
		objects{end+1} = new_plane([0.4 0.1 0],[1 3 2],@color_pattern);
		objects{end+1} = new_algebraic([0 0 0.5],@polynomial_torus2(0.6,0.2,-2*pi/5),@color_light);
	elseif build == 6
		objects{end+1} = new_algebraic([0 0 0],@polynomial_flower(),@color_light);
		objects{end+1} = new_sphere([0 0 1],0.8,@color_light);
	elseif build == 7
		objects{end+1} = new_sphere([-0.05 -0.65 -0.25],0.1,@color_light);
		objects{end+1} = new_sphere([-0.35 -0.65 -0.05],0.1,@color_light);
		objects{end+1} = new_sphere([-0.1 -0.4 0],0.35,@color_light);
		objects{end+1} = new_sphere([-0.1 0.1 0],0.3,@color_light);
		objects{end+1} = new_sphere([0 0.5 0.1],0.25,@color_light);
		objects{end+1} = new_sphere([0.2 0.7 0.2],0.15,@color_light);
	end
end


% end of code

