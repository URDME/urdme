close all;

%// center
c = [3 3];

%// number of points
n = 1000;

%// running variable
t = linspace(0,2*pi,n);

%----Proliferating----

%// radius
r = 2;

x = c(1) + r*sin(t);
y = c(2) + r*cos(t);

%// draw polygon 
patch(x,y,graphics_color('vermillion'),'LineStyle','none');

%----Quiescent----

%// radius
r = 1.5;

x = c(1) + r*sin(t);
y = c(2) + r*cos(t);

%// draw polygon 
patch(x,y,graphics_color('bluish green'),'LineStyle','none');

%----Necrotic----

%// radius
r = 0.75;

x = c(1) + r*sin(t);
y = c(2) + r*cos(t);

%// draw polygon 
patch(x,y,[0 0 0],'LineStyle','none');

axis equal
axis off