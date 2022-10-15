clc;
clear;
close all;
%dimensiones del dominio
nx=200; ny=150;
[XX,YY]=meshgrid(1:nx,1:ny); 
%velocidades, presión, densidad
[ux,uy,p,rho]=deal(zeros(ny,nx));
%iteraciones
t=[0:2400];
N=length(t);
dt=0.001;
%%
fig=figure(1);
% densidad que se añadirá
rho0=0.01;
pause(2)
for i=1:N-1;
    %paso 1: divergencia de la velocidad
    r=divergence(ux,uy);
    %paso 2: resolver ecuación de Poisson
    p=fMa(r);
     %paso 3:condiciones de contorno para la presión
    p(1:ny,1)=p(1:ny,2);p(1,1:nx)=p(2,1:nx);p(1:ny,nx)=p(1:ny,nx-1);p(ny,1:nx)=p(ny-1,1:nx);
   %paso 4: teorema de Helmholtz para incompresibilidad de la velocidad
    [dpx,dpy]=gradient(p);
    ux=ux-dpx;
    uy=uy-dpy;
    %paso 5: cálculo de las celdas hacia donde apuntan las velocidades e
    %interpolación para distribuir adecuadamente las magnitudes relevantes
    [XXv,YYv]=RK4(XX,YY,ux,uy,-1);
    ux=interp2(XX,YY,ux,XXv,YYv);
    uy=interp2(XX,YY,uy,XXv,YYv);
    rho=interp2(XX,YY,rho,XXv,YYv);
     %paso 6: condiciones de contorno para la velocidad
        ux(1,1:nx)=0;
    ux(1:ny,1)=0;
    ux(ny,1:nx)=0;
    ux(1:ny,nx)=0;
    uy(1,1:nx)=0;
    uy(1:ny,1)=0;
    uy(ny,1:nx)=0;
    uy(1:ny,nx)=0;
% paso 7: añadimos la densidad del flujo incidente y su velocidad
    ux(ny/2-5:ny/2+5,10:15)=0.9;
    rho(ny/2-5:ny/2+5,10:15)=rho0;
    %paso 9: dibujar la densidad de cada iteración: paso más lento...
pcolor(XX,YY,rho)
colormap(jetvar)
shading interp 
    axis([1 nx 1 ny])
drawnow
end