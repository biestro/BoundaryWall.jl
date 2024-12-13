function [y,yp] = Po(z,a)
% Segunda solucion para la ecuacion parabolica y''+(1/4*z^2-a)y=0
% a es un parametro real o complejo
% z es el vector/matriz de valores deseados reales o complejos 
%(si la entrada es una matriz/vector complejo todos los numeros tienen que tener la misma phase).
% y es la funcion evaluada en z
% yp es la derivada de la funcion en z
% 
% La solución es por series con continuación analítica
% % Julio C. Gutierrez

h=max(max(abs(z)))/4000; % longitud paso para la continuacion analitica
% Se puede poner otro numero mas grande en lugar de 4000 si se desea mas precision, 
% 4001 es el numero de continuaciones analiticas que se haran en el intervalo 

a=-a; % conversion por definicion de la serie

%%%%Convertimos z a un vector
[largo,ancho]=size(z);
z=transpose(z(:));

zmax=max(abs(z));
nca2=size(z,2);
zf=sum(angle((-(real(z)<0)+(real(z)>0)).*z))/nca2;
index_paridad=-(z/exp(i*zf)<0)+(z/exp(i*zf)>=0); % Vector para despues para obtener los valores negativos de la derivada
nca=floor(zmax/h); % numero de continuaciones analiticas 
[~,index]=sort(abs(z)); % se ordenan los valores
zc=z(index);
z=abs(z);

% CARACTERISTICAS DE LA FUNCION
% xo=0; % punto de arranque
%condiciones de frontera en xo=0
%yi=[1,0]; % PARA FUNCION PAR
yi=[0,1]; %PARA FUNCION IMPAR
y=[];
yp=[];

for j=1:nca
    xo=(j-1)*h;
    xf=j*h;
    zi=z(xo<=z & z<xf); % puntos a evaluar en el intervalo deceado
    zi=[zi,xf]; %se agrega punto para la proxima expancion
    [yin,ypin]=yo_parabolic(yi,xo*exp(i*zf),   a,zi*exp(i*zf));
    pint=length(yin)-1; %numero de puntos en el intervalo 
    % se actualizan valores
    y=[y,yin(1:pint)]; 
    yp=[yp,ypin(1:pint)];       
    xo=xo+h;
    yi=[yin(pint+1),ypin(pint+1)]; 
end

%Se agrega el ultimo punto
if zmax==abs( nca*h-max(abs(z)))<10^(-10)&&(nca2~=1)
    num=sum((z==zmax));
    y=[y,yin(pint+1)*ones(1,num)];
    yp=[yp,ypin(pint+1)*ones(1,num)];
end    
    

if nca2~=1
    yv(index)=y;
    yv=yv.*index_paridad;
    ypv(index)=yp;
    ypv=ypv;
else
    if z~=0  
        yv=yin(1)*index_paridad(1);
        ypv=ypin(1);
    else
        yv=1;
        ypv=0;
    end
end

%%%Regresa el vector a matriz
y=reshape(yv,[largo,ancho]);
yp=reshape(ypv,[largo,ancho]);
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,yp] = yo_parabolic(yi,xo,a,r)
% xo  punto a partir del cual se hace la expansion
% r punto apartir de xo para el cual se quiere el valor y(xo+r)
% yl es la funcion evaluada en xo+r
%yi=[ac1,ac2] condiciones iniciales, primera y segunda derivada
% Se utiliza para hace la continuacion analitica y obtener y1 y y2

N=100; %numero par de coeficientes en la expansion

% Calcula los coeficientes
ac=parabolic_coe(yi,xo,a,N);  %se calculan los coeficientes
acn=ac./gamma(1:N); %calcula coeficiente para y
acp=ac(2:N)./gamma(1:N-1); % calcula coeficientes para y'

% Evaluacion de los polinomios
y=polyval(acn(N:-1:1),(r-xo));
yp=polyval(acp(N-1:-1:1),(r-xo));

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yc = parabolic_coe(yi,xo,a,N)
% Funcion que calcula primeros N coefficientes para la expansion de las
% funciones parabolicas y1 
% yi=[y(xo),y'(xo),y''(xo)]

yc=zeros(1,N);
yi=[yi,-((xo^2)/4+a)*yi(1)]; % se complementan las condiciones iniciales el valor de las segunda derivada en el punto

yc(1:3)=yi(1:3);
yc(4)=-(4*a+xo^2)/4*yc(2)-xo/2*yc(1);

for n=4:N-1
    % recorre vector para que n=0, sea la pocision 1
    yc(n+1)=-(4*a+xo^2)/4*yc(n+1-2)-(n-2)/2*xo*yc(n+1-3)-(n-2)*(n-3)/4*yc(n+1-4);
end
end