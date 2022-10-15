function sol=fMa(lapp);
%resolver la ecuación de Poisson definiendo las diagonales relevantes y
%aplicando la inversa. Sería más sencillo probablemente aplicar
%Gauss-Seidel "a lo bruto" con la presión
%p(i,j)=(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))/4+r(i,j), donde r=div(ux,uy)
%pero perdería precisión. Aunque funciona bien y probablemente sea más
%rápido.
[nx,ny]=size(lapp);
% diagonales importantes:
boundary_index1=[1:nx];boundary_index2=[1:nx:1+(ny-1)*nx];boundary_index3=[1+(ny-1)*nx:nx*ny];boundary_index4=[nx:nx:nx*ny];
diagonals=[-4*ones(nx*ny,1),ones(nx*ny,4)];
A=spdiags(diagonals,[0 -1 1 -nx nx],nx*ny,nx*ny);
%diagonales de los bordes
diagonals1=[-3*ones(nx*ny,1),ones(nx*ny,3)];
I=spdiags(diagonals1,[0 -1 1  nx],nx*ny,nx*ny);
A(boundary_index1,:)=I(boundary_index1,:);
I=spdiags(diagonals1,[0  1  -nx nx],nx*ny,nx*ny);
A(boundary_index2,:)=I(boundary_index2,:);
I=spdiags(diagonals1,[0  -1 1  -nx ],nx*ny,nx*ny);
A(boundary_index3,:)=I(boundary_index3,:);
I=spdiags(diagonals1,[0  -1   -nx nx],nx*ny,nx*ny);
A(boundary_index4,:)=I(boundary_index4,:);
b=zeros(nx,ny);
b=lapp;
b=reshape(b,nx*ny,1);
sol=A\b;
sol=reshape(sol,nx,ny);
end