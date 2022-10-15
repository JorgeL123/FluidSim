clc;
clear;
close all;
%dimensiones del dominio
nx=400; ny=90;
[XX,YY]=meshgrid(1:nx,1:ny); 
%velocidades, presión, densidad
[ux,uy,p,rho]=deal(zeros(ny,nx));
%iteraciones
t=[0:870];
N=length(t);
dt=0.001;
%%
fig=figure(1);
% densidad y dimensiones de la placa y del cuadrado
rho0=0.1; f=18;d=3;L=5;
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
%paso 7: condiciones de contorno en las paredes del tubo 
%pared inferior
ux(ny/2-f-6:ny/2-f-1,8:nx)=0; uy(ny/2-f-6:ny/2-f-1,8:nx)=0; p(ny/2-f-1,8:nx)=p(ny/2-f,8:nx);
rho(ny/2-f-6:ny/2-f-1,8:nx)=0;
%pared superior
ux(ny/2+f+1:ny/2+f+6,8:nx)=0; uy(ny/2+f+1:ny/2+f+6,8:nx)=0; p(ny/2+f+1,8:nx)=p(ny/2+f,8:nx);
rho(ny/2+f+1:ny/2+f+6,8:nx)=0;
%pared trasera
ux(ny/2-f-6:ny/2+f+6,5:8)=0;uy(ny/2-f-6:ny/2+f+6,5:8)=0; p(ny/2-f-6:ny/2+f+6,8)=p(ny/2-f-6:ny/2+f+6,9);
rho(ny/2-f-6:ny/2+f+6,5:8)=0;
%condiciones de contorno en las paredes del cuadrado
uy(ny/2-d:ny/2+d,25:25+L)=0;ux(ny/2-d:ny/2+d,25:25+L)=0; p(ny/2-d,25:25+L)=p(ny/2-d-1,25:25+L);p(ny/2+d,25:25+L)=p(ny/2+d+1,25:25+L);p(ny/2-d:ny/2+d,25)=p(ny/2-d:ny/2+d,24);p(ny/2-d:ny/2+d,25+L)=p(ny/2-d:ny/2+d,25+1+L);
rho(ny/2-d:ny/2+d,25:25+L)=0;
% paso 8: añadimos la densidad del flujo incidente y su velocidad
ux(ny/2-f:ny/2+f,10:22)=3.5;
rho(ny/2-f:ny/2+f,10:22)=rho0;
%paso 9: dibujar la densidad de cada iteración: paso más lento...
pcolor(XX,YY,rho)
colormap(jetvar)
shading interp
hold on 
plot([25 30 30 25], [ny/2-d ny/2-d ny/2+d ny/2+d],'-k')
area([25 30],[ny/2+d ny/2+d],ny/2-d,'Facecolor','black')
    axis([10 nx-150 ny/2-f ny/2+f]) 
    drawnow
end