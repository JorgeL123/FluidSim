clc;
clear;
close all;
nx=200; ny=150;
% Define sizes, positions and where the new particles will appear in the
% flow
[XX,YY]=meshgrid(1:nx,1:ny);
[ux,uy,p]=deal(zeros(ny,nx));
pxo=XX(nx/2-5:nx/2+5,10:15);
pyo=YY(ny/2-5:ny/2+5,10:15);
t=[0:3500];
N=length(t);
px=pxo;
py=pyo;
%%
fig=figure(1);
pause(2)
for i=1:N-1;
    % Velocity of new particles
    ux(ny/2-5:ny/2+5,10:15)=0.9;
    % Poisson Equation 
    r=divergence(ux,uy);
    p=fMa(r);
    % Pressure Boundary Conditions
    p(1:ny,1)=p(1:ny,2);p(1,1:nx)=p(2,1:nx);p(1:ny,nx)=p(1:ny,nx-1);p(ny,1:nx)=p(ny-1,1:nx);
   % Apply Helmholtz theorem
    [dpx,dpy]=gradient(p);
    ux=ux-dpx;
    uy=uy-dpy;
    % Interpolation of velocities and positions to see where the flow is
    % moving towards
    [XXv,YYv]=RK4(XX,YY,ux,uy,-1);
    ux=interp2(XX,YY,ux,XXv,YYv);
    uy=interp2(XX,YY,uy,XXv,YYv);
    % Boundary Conditions
    ux(1,1:nx)=0;
    ux(1:ny,1)=0;
    ux(ny,1:nx)=0;
    ux(1:ny,nx)=0;
    uy(1,1:nx)=0;
    uy(1:ny,1)=0;
    uy(ny,1:nx)=0;
    uy(1:ny,nx)=0;
    % Chain new particles showing up in the flowing fluid
   [px,py]=RK4(px,py,ux,uy,1);
    px=[px,pxo];
    py=[py,pyo];

    % Draw the particles
    plot(px,py,'.b','MarkerSize',0.07)
    axis([1 nx 1 ny])

drawnow
end