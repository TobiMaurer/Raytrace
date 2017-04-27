% octave -q raytrace.m


% define sets


function contained = contained1(coords) % cube/star1
end

function contained = contained2(coords) % torus/pipe
end

function contained = contained3(coords) % cone
end

function contained = contained4(coords) % star2
end

function contained = contained_test(coords)
	contained = norm(coords) < 0.5;
	%return
	
	roty = pi/6;
	rotx = pi/6;
	coords = coords*[cos(roty),0,sin(roty);0,1,0;-sin(roty),0,cos(roty)];
	coords = coords*[1,0,0;0,cos(rotx),sin(rotx);0,-sin(rotx),cos(rotx)];
	
	if norm(coords) == 0
		contained = false;
		return
	end
	
	direction = coords/norm(coords); %# div 0
	dist = min([ ...
		norm(direction*[0,0,0;0,1,0;0,0,1]), ...
		norm(direction*[1,0,0;0,0,0;0,0,1]), ...
	]);
	
	factor = ((dist^3)/3-(dist^2)/sqrt(8))*2+1;
	contained = norm(coords(3)) < 0.35 && norm(coords) < factor^5;
	
	%contained = norm(coords,inf) < 0.35;
end


% define materials


function color = color1(coords,normal,curvature) % shiny1
end

function color = color2(coords,normal,curvature) % furry
end

function color = color3(coords,normal,curvature) % checkered
end

function color = color4(coords,normal,curvature) % shiny2
end

function color = color_test(coords,normal,curvature)
	color = 1+normal(3);
end


% define scene


clear objects;
clear ambient;
clear parallel;


objects(1).type = "sphere";
objects(1).position = [0,0,0];
objects(1).radius = 1;
objects(1).transform = eye(4);
objects(1).color = @color_test;

objects(2).type = "plane";
objects(2).normal = [0,0,0];
objects(2).distance = 0.5;
objects(2).transform = eye(4);
objects(2).color = @color_test;

objects(3).type = "set";
objects(3).contained = @contained_test;
objects(3).transform = eye(4);
objects(3).color = @color_test;


ambient.intensity = 0.5;

parallel.intensity = 0.5;
parallel.direction = [0,0,0];


% end of scene definition


% helper functions


function transform = scale(factor)
end

function transform = translate()
end

function transform = rotate(pivot,angle)
end


% initialize variables


global rows = [-1:2/40:1];
global cols = [-1:2/40:1];
global time = [-1:2/20:1];
global appx = 10;

depths = ones(numel(rows),numel(cols)) * inf;
removes = zeros(numel(rows),numel(cols));
colors = cell(numel(rows),numel(cols));
normals = zeros(numel(rows),numel(cols),3);
curvatures = zeros(numel(rows),numel(cols));
pixels = zeros(numel(rows),numel(cols));


% determine intersections


function depth = depth_sphere(object,row,col)
	global rows
	global cols
	depth = inf;
	
	if object.transform != 0 %# implement
	end
end

function depth = depth_plane(object,row,col)
	global rows
	global cols
	depth = inf;
	
	if object.transform != 0 %# implement
	end
end

function depth = depth_set(object,row,col)
	global rows
	global cols
	global time
	global appx
	depth = inf;
	
	depth2 = -inf;
	for depth1 = time
		
		coords = [cols(col),rows(row),depth1];
		if object.transform != 0
			coords *= object.transform; %# left mul
		end
		
		if object.contained(coords)
			for i = 1:appx
				
				depth3 = depth1 + depth2;
				depth3 /= 2;
				
				coords = [cols(col),rows(row),depth3];
				if object.contained(coords)
					depth1 = depth3;
				else
					depth2 = depth3;
				end
			end
		
			depth = depth1;
			
			break;
		end
		
		depth2 = depth1;
	end
end

for row = 1:numel(rows)
	for col = 1:numel(cols)
		for object = objects
			
			depth = inf;	
			
			if strcmp(object.type,"sphere")
				depth = depth_sphere(object,row,col);
			elseif strcmp(object.type,"plane")
				depth = depth_plane(object,row,col);
			elseif strcmp(object.type,"set")
				depth = depth_set(object,row,col);
			end
			
			if depth != inf
				depths(row,col) = depth;
				colors{row,col} = object.color;
			
				break;
			end
			
		end
	end
end


% determine normals for intersections


for row = 1:numel(rows)
	for col = 1:numel(cols)
		
		if depths(row,col) != inf
			
			depth = [inf,inf,inf,inf];
			normal1 = [0,-1;0,-1]; % x,z;y,z
		
			if col > 1
				depth(1) = depths(row,col - 1);
				normal1(1,:) = [-1,0];
			end
		
			if col < numel(cols)
				depth(2) = depths(row,col + 1);
				normal1(1,:) = [1,0];
			end
		
			if row > 1
				depth(3) = depths(row - 1,col);
				normal1(2,:) = [-1,0];
			end
		
			if row < numel(rows)
				depth(4) = depths(row + 1,col);
				normal1(2,:) = [1,0];
			end
			
			if (depth(1) == inf && depth(2) == inf) ...
			|| (depth(3) == inf && depth(4) == inf)
				removes(row,col) = true;
			else
				
				if depth(1) != inf && depth(2) != inf
					normal1(1,:) = [ ...
						depth(2) - depth(1), ...
						cols(col - 1) - cols(col + 1) ...
					]; % x=z2-z1,z=x1-x2
					normal1(1,:) /= norm(normal1(1,:)); %# div 0
				end
				
				if depth(3) != inf && depth(4) != inf
					normal1(2,:) = [ ...
						depth(4) - depth(3), ...
						rows(row - 1) - rows(row + 1) ...
					]; % y=z2-z1,z=y1-y2
					normal1(2,:) /= norm(normal1(2,:)); %# div 0
				end
				
				normal2 = [ ...
					normal1(1,1), ...
					normal1(2,1), ...
					min(normal1(1,2),normal1(2,2)) ...
				]; %# improve
				
				
				normal2x = [ ...
					normal1(1,1), ...
					0, ...
					normal1(1,2) ...
				]; %# improve
				
				
				normal2 /= norm(normal2); %# div 0
				
				normals(row,col,:) = normal2;
				
			end
			
		end
	end
end


% remove if normal undefined


for row = 1:numel(rows)
	for col = 1:numel(cols)
		
		if removes(row,col)
			
			depths(row,col) = inf;
			
		end
		
	end
end


% determine curvature for intersections


for row = 1:numel(rows)
	for col = 1:numel(cols)
		
		if depths(row,col) != inf
			
			curvatures(row,col) = 0;
			
		end
	end
end


% determine color for intersections


for row = 1:numel(rows)
	for col = 1:numel(cols)
		
		if depths(row,col) != inf
			
			coords = [cols(col),rows(row),depths(row,col)];
			
			pixels(row,col) = colors{row,col} ...
				(coords,normals(row,col,:),curvatures(row,col));
			
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

