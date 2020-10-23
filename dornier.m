clear
clc
%datos
S=33.93;
Iy=18900;
mac=2.05;
m=5457;
g=9.81;
h=3048;%m
rho=1.225*(1-22.557e-6*h)^4.256
V0=90;
alpha=0;%ang incidencia en rad
tita_e=0;%ang planeo  en rad
M=V0/340;
W=m*g
Sh=8.33;
a=5.837;%wf  aw=5.108
a1=5.73;
Xcg_mac=0.3;
Xca_mac=0.1401;
L_h=9.2;
lt=L_h-(Xcg_mac-Xca_mac)*mac
grad_downwash=0.39;
Vol_cola=L_h*Sh/(mac*S)
Kn=0.2397;
CL=(2*W)/(rho*V0^2*S)
CD0=0.029;
K1=-0.00593;
K2=0.04464;
CD=CD0+K1*CL+K2*CL^2;
dCD_dv=0;
dCm_dv=(4*Kn*W)/(rho*S*V0^3);
dT_dv=0;
dCD_dalpha=K1*a+2*K2*CL*a
dCL_dv=(M^2*CL)/(1-M^2);

%calculo de las derivadas de estabilidad
Xu=-2*CD-V0*dCD_dv+(1/(0.5*rho*V0*S))*dT_dv;
Xw=CL-dCD_dalpha;
Xq=0;
Xdotw=0;
Zu=-2*CL-V0*dCL_dv;
Zw=-CD-a;
Zq=-Vol_cola*a1;
Zdotw=Zq*grad_downwash;
Mu=V0*dCm_dv;
Mw=-a*Kn;
Mq=Zq*lt/mac;
Mdotw=Mq*grad_downwash;


%parametros
mhu=m/(0.5*rho*S*mac);
sigma=m/(0.5*rho*V0*S);
iy=Iy/(m*mac^2);
%matrices
M=[1    -(Xdotw/mhu)   0      0;
    0  (1-Zdotw/mhu)   0      0;
    0  -(Mdotw/mhu) (iy/mhu)  0;
    0    0          0      1];
Aprima=[Xu   Xw   (Xq/mhu)   -(CL*cos(tita_e));
        Zu   Zw    1+Zq/mhu    -(CL*sin(tita_e));
        Mu   Mw   (Mq/mhu)        0;
        0    0      1             0];
A=inv(M)*Aprima;
[V,D]=eig(A);%calcula valores y vectores propios
%valores propios dim
val_prop=[D(1),D(6),D(11),D(16)]/sigma
%vectores propios
V1=V(:,1);
V2=V(:,2);
V3=V(:,3);
V4=V(:,4);
%condiciones iniciales
X_inicial=[0.1;0.1;1;0];
Z_inicial=inv(V)*X_inicial;

%modo corto periodo
U_sp=@(t) V0*real(V1(1)*Z_inicial(1)*exp(val_prop(1)*t)+V2(1)*Z_inicial(2)*exp(val_prop(2)*t));
subplot(2,2,1) 
fplot(U_sp,[0,2]);title('modo corto periodo');xlabel('t(s)'); ylabel('U(m/s)');grid on

W_sp=@(t) V0*real(V1(2)*Z_inicial(1)*exp(val_prop(1)*t)+V2(2)*Z_inicial(2)*exp(val_prop(2)*t));
subplot(2,2,2)
fplot(W_sp,[0,2]);xlabel('t(s)'); ylabel('W(m/s)');grid on

Q_sp=@(t) (1/sigma)*real(V1(3)*Z_inicial(1)*exp(val_prop(1)*t)+V2(3)*Z_inicial(2)*exp(val_prop(2)*t));
subplot(2,2,3)
fplot(Q_sp,[0,2]);xlabel('t(s)'); ylabel('Q(deg/s)');grid on

tita_sp=@(t) real(V1(4)*Z_inicial(1)*exp(val_prop(1)*t)+V2(4)*Z_inicial(2)*exp(val_prop(2)*t));
subplot(2,2,4)
fplot(tita_sp,[0,2]);xlabel('t(s)'); ylabel('tita(m/s)');grid on
amortiguamiento_sp=real(val_prop(1))
frecuencia_sp=imag(val_prop(1))

%modo fugoide
figure(2)
U_ph=@(t) V0*real(V3(1)*Z_inicial(3)*exp(val_prop(3)*t)+V4(1)*Z_inicial(4)*exp(val_prop(4)*t));
subplot(2,2,1)
fplot(U_ph,[0,600]);title('modo fugoide');xlabel('t(s)'); ylabel('U(m/s)');grid on

W_ph=@(t) V0*real(V3(2)*Z_inicial(3)*exp(val_prop(3)*t)+V4(2)*Z_inicial(4)*exp(val_prop(4)*t));
subplot(2,2,2)
fplot(W_ph,[0,600]);xlabel('t(s)'); ylabel('W(m/s)');grid on

Q_ph=@(t) (1/sigma)*real(V3(3)*Z_inicial(3)*exp(val_prop(3)*t)+V4(3)*Z_inicial(4)*exp(val_prop(4)*t));
subplot(2,2,3)
fplot(Q_ph,[0,600]);xlabel('t(s)'); ylabel('Q(deg/s)');grid on

tita_ph=@(t) real(V3(4)*Z_inicial(3)*exp(val_prop(3)*t)+V4(4)*Z_inicial(4)*exp(val_prop(4)*t));
subplot(2,2,4)
fplot(tita_ph,[0,600]);xlabel('t(s)'); ylabel('tita(m/s)');grid on

amortiguamiento_ph=real(val_prop(3))
frecuencia_ph=imag(val_prop(3))







