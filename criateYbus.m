function[Ybus,de,para,barra,lin,re,med,ybus]=criateYbus(barra,lin,shunt,conf,Sbase,traf,reg,med)
 lin1=lin;


  confN=conf(1,:);
  conf=str2double(conf(2:length(conf(:,1)),:));
 %Conversão pés=> km, e milha => km
 ftkm=0.0003048;
 milekm= 1.609344;
 %Atualização das matrizes lin e conf para Km
 conf(:,3:length(conf(1,:)))= conf(:,3:length(conf(1,:)))/milekm;
 lin(:,4)=lin(:,4)*ftkm;
 
 %Definições 
 nb= length(barra(:,1));     %num. de barras
 nl= length(lin(:,3));       %num. de linhas
 nome=barra(:,1);            %nome original da barra
 if any(traf) ~=0; 
 tapt=traf(:,[4 5 6]);
 faset=traf(:,[7 8 9]);
 end
 for i=1:nl                  %Atualização dos nomes(n°) das barras
     for j=1:nb
         if lin(i,1) == nome(j);
             de(i)= barra(j,2)';
         end
         
         if lin(i,2) == nome(j);
             para(i)= barra(j,2)';
         end
     end
 end

 
%%
 %Regulator data
if reg(1,1)~=0
vmax=1e12;
vmin=1e-10;
V=reg(:,16:18);

for i=1:length(reg(:,1))
    reg_de(i)=find(barra(:,1)==reg(i,2)); %barra de do regulador
    reg_para(i)=find(barra(:,1)==reg(i,3)); %barra para do regulador
    reg_dec(i)=find(barra(:,1)==reg(i,4));  %barra de - controle
    reg_parac(i)=find(barra(:,1)==reg(i,5)); %barra para - controle (barra que está sendo controlada)
    reg_nb(i,1)=length(barra(:,1))+i;   %aumento do numero de barras
    switch reg(i,6)
       case 123
            reg_ph(i,1:3)=[1 1 1];
        case 12
            reg_ph(i,1:3)=[1 1 0];
        case 13
            reg_ph(i,1:3)=[1 0 1];
        case 23
            reg_ph(i,1:3)=[0 1 1];
        case 1
            reg_ph(i,1:3)=[1 0 0];
        case 2
            reg_ph(i,1:3)=[0 1 0];
        case 3
            reg_ph(i,1:3)=[0 0 1];
    end;
    switch reg(i,22)
        case 123
            reg_mph(i,1:3)=[1 1 1];
        case 12
            reg_mph(i,1:3)=[1 1 0];
        case 13
            reg_mph(i,1:3)=[1 0 1];
        case 23
            reg_mph(i,1:3)=[0 1 1];
        case 1
            reg_mph(i,1:3)=[1 0 0];
        case 2
            reg_mph(i,1:3)=[0 1 0];
        case 3
            reg_mph(i,1:3)=[0 0 1];
    end;
    reg_band(i,1)=reg(i,7)/120;
    z=2e-5 +3e-4i;
 
    reg_r(i,[1 2 3])=[0 0 0];
    reg_x(i,[1 2 3])=[0 0 0];
    reg_g(i,[1 2 3])=[vmax vmax vmax];
    reg_b(i,[1 2 3])=[vmax vmax vmax];
    reg_v(i,[1 2 3])=[1 1 1];    
switch reg(i,6)
       
    case 123
        reg_r(i,1)=real(z);
        reg_x(i,1)=imag(z);
        reg_g(i,1)=real(inv(z));
        reg_b(i,1)=imag(inv(z)); 
        reg_v(i,1)=V(i,1)/120;
  
        reg_r(i,2)= real(z); 
        reg_x(i,2)= imag(z);
        reg_g(i,2)= real(inv(z)); 
        reg_b(i,2)= imag(inv(z));
        reg_v(i,2)= V(i,2)/120;

        reg_r(i,3)= real(z);
        reg_x(i,3)= imag(z);
        reg_g(i,3)= real(inv(z));
        reg_b(i,3)= imag(inv(z));
        reg_v(i,3)= V(i,3)/120;
    
    case 12
        reg_r(i,1)=real(z);
        reg_x(i,1)=imag(z);
        reg_g(i,1)=real(inv(z));
        reg_b(i,1)=imag(inv(z)); 
        reg_v(i,1)=V(i,1)/120;
  
        reg_r(i,2)= real(z); 
        reg_x(i,2)= imag(z);
        reg_g(i,2)= real(inv(z)); 
        reg_b(i,2)= imag(inv(z));
        reg_v(i,2)= V(i,2)/120;
    case 13
         reg_r(i,1)=real(z);
        reg_x(i,1)=imag(z);
        reg_g(i,1)=real(inv(z));
        reg_b(i,1)=imag(inv(z)); 
        reg_v(i,1)=V(i,1)/120;
  
        reg_r(i,3)= real(z);
        reg_x(i,3)= imag(z);
        reg_g(i,3)= real(inv(z));
        reg_b(i,3)= imag(inv(z));
        reg_v(i,3)= V(i,3)/120;
    case 1
        reg_r(i,1)=real(z);
        reg_x(i,1)=imag(z);
        reg_g(i,1)=real(inv(z));
        reg_b(i,1)=imag(inv(z)); 
        reg_v(i,1)=V(i,1)/120;
  
    case 2
        reg_r(i,2)= real(z); 
        reg_x(i,2)= imag(z);
        reg_g(i,2)= real(inv(z)); 
        reg_b(i,2)= imag(inv(z));
        reg_v(i,2)= V(i,2)/120;
    case 3
        reg_r(i,3)= real(z);
        reg_x(i,3)= imag(z);
        reg_g(i,3)= real(inv(z));
        reg_b(i,3)= imag(inv(z));
        reg_v(i,3)= V(i,3)/120;
end
 reg_tap(i,[1 2 3])=[reg(i,19) reg(i,20) reg(i,21)];
end
end


for k=1:nl
    a=find(conf(:,1) == lin(k,5));
    if lin(k,5)~=0
        Z(k,:)=[conf(a,2) conf(a,3:(length(conf)))*lin(k,4)];
    else
        Z(k,:)=[conf(a,2) ((conf(a,3:(length(conf(1,:))))*milekm)*Sbase/((lin(k,6))))*1e5];
        end
end


 
  %Criação da Ybus
phase=Z(:,1)'; %Fases ABC ou 123 de cada linha
Z=[Z(:,2:length(Z(1,:)))] ;
r=[Z(:,find(confN=='Raa')-2) Z(:,find(confN=='Rab')-2) Z(:,find(confN=='Rac')-2) Z(:,find(confN=='Rbb')-2) Z(:,find(confN=='Rbc')-2) Z(:,find(confN=='Rcc')-2)]./((lin(:,6)).^2/Sbase);
x=[Z(:,find(confN=='Xaa')-2) Z(:,find(confN=='Xab')-2) Z(:,find(confN=='Xac')-2) Z(:,find(confN=='Xbb')-2) Z(:,find(confN=='Xbc')-2) Z(:,find(confN=='Xcc')-2)]./((lin(:,6)).^2/Sbase);
bs=[Z(:,find(confN=='Baa')-2) Z(:,find(confN=='Bab')-2) Z(:,find(confN=='Bac')-2) Z(:,find(confN=='Bbb')-2) Z(:,find(confN=='Bbc')-2) Z(:,find(confN=='Bcc')-2)].*((lin(:,6)).^2/Sbase)*1e-6;

%Tap transformador TAP fixo 

for i=1:nl
   if lin(i,5) == 0 
    TAPm(i,[1 2 3]) = tapt;  %%Colocar aqui o MÓDULO DO TAP DA FASE A, B E C.
    TAPp(i,[1 2 3]) = faset;
   else
    TAPm(i,[1 2 3]) = [1 1 1];
    TAPp(i,[1 2 3]) = [0 0 0];
   end
end

 %Aumento do número de barras devido aos reguladores;
if reg(1,1)~=0
    barra2=barra;
for j=1:length(reg(:,1))
i=nl+j;
lin_tipo(i)=2;
n_lin_reg(j,1)=i;    %branch do regulador
b_de=(find(de==reg_de(j)));
b_para=(find(para==reg_para(j)));
b=intersect(b_de,b_para);
if isempty(b)
    bf=find(de==reg_para(j));
    bt=find(para==reg_de(j));
    b=intersect(bf,bt);
end
if de(b)==reg_dec(j)
    de(b)=reg_nb(j);
else
    para(b)=reg_nb(j);
end
    de(i)=reg_dec(j);
    para(i)=reg_nb(j);
    phase(i)=phase(b);
    r(i,[1 4 6])=reg_r(j,:).*reg_ph(j,1) + (reg_ph(j,1)==0)*vmin;
    x(i,[1 4 6])=reg_x(j,:).*reg_ph(j,:);
    bs(i,[1 2 3 4 5 6])=[0 0 0 0 0 0];
    TAPm(i,:)=1+0.00625*reg(j,[19 20 21]);
    TAPp(i,:)=[0 0 0];
    
   barra2= [barra2 ; barra(reg_parac(j),1:length(barra(1,:))-1) 1 ];

   
end
end

barra2([reg_nb],:)=[reg_nb reg_nb barra2([reg_nb],3:length(barra2(1,:)))];
barra=barra2;


for a=1:length(reg_nb)
    for b = 1:3
     med(length(med(:,1))+[1 2],:)=[(length(med(:,1))+1) 0 20 b reg_nb(a) reg_nb(a);(length(med(:,1))+1) 0 21 b reg_nb(a) reg_nb(a)];
    % med(length(med(:,1))+[1 2],:)=[(length(med(:,1))+1) 0 1 b reg_nb(a) reg_nb(a);(length(med(:,1))+1) 0 2 b reg_nb(a) reg_nb(a)];
    end
end
%Tipo igual a 20: representa uma medida fictícia, devido ao regulador;
%%
%Mudança das linhas em função do regulador
for i=1:length(reg(:,1))
    a(i)=find(lin(:,1)==reg(i,2) & lin(:,2)==reg(i,3));
    lin_reg(i,:)=lin(a(i),:);
end
lin2=[lin;lin_reg];

lin2(:,[1 2])=[de' para'];
for i=1:length(barra(:,1))
    for j=1:length(de)
        if de(j)==barra(i,2)
            lin2(j,1)=barra(i,1);
        end
        if para(j)==barra(i,2)
            lin2(j,2)=barra(i,1);
        end
    end
end
lin=lin2;
    



n=nb+length(reg(:,1));
Ybus=sparse([],[],[],n*3,n*3);
% Colocando os Elementos Shunt na matriz
for i=1:length(shunt(:,1))
    for j=1:nb
        if shunt(i,1) == barra(j,1)
    Ybus([1 2 3]+(j-1)*3,[1 2 3]+(j-1)*3)=diag([(shunt(i,2)) (shunt(i,3)) (shunt(i,4))])*1i/((Sbase*10^-3)/3);
        end
    end
end
vmax=1e12;

nl=nl+length(reg(:,1));
a=find(lin(:,5)== 0);
lin_tipo(a) = 2;

ybus.lin_tipo=lin_tipo;
ybus.phase=phase;
ybus.r=r;
ybus.x=x;
ybus.bs=bs;
ybus.TAPm=TAPm;
ybus.TAPp=TAPp;
ybus.phase=phase;
ybus.de=de;
ybus.para=para;
for i=1:nl
    if lin_tipo(i) ~= 2
        switch phase(i)
            case 1
                z = diag([r(i,1)+x(i,1)*1i vmax vmax]);
                b = diag([bs(i,1)*1i 0 0]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
            case 2
                z = diag([vmax r(i,4)+x(i,4)*1i vmax]);
                b = diag([0 bs(i,4)*1i 0]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
            case 3
                z = diag([vmax vmax r(i,6)+x(i,6)*1i]);
                b = diag([0 0 bs(i,6)*1i]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
            case 12
                z=diag([0 0 vmax]);
                z(1:2,1:2) = [r(i,1:2)+x(i,1:2)*1i; r(i,2)+x(i,2)*1i r(i,4)+x(i,4)*1i];
                b=diag([0 0 0]);
                b(1:2,1:2) = [bs(i,1:2)*1i; bs(i,2)*1i bs(i,4)*1i];
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
            case 13
                z=diag([0 vmax 0]);
                z([1 3],[1 3]) = [r(i,1)+x(i,1)*1i r(i,3)+x(i,3)*1i; r(i,3)+x(i,3)*1i r(i,6)+x(i,6)*1i];
                b=diag([0 0 0]);
                b([1 3],[1 3]) = [bs(i,1)*1i bs(i,3)*1i; bs(i,3)*1i bs(i,6)*1i];
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
            case 23
                z=diag([vmax 0 0]);
                z(2:3,2:3) = [r(i,4:5)+x(i,4:5)*1i; r(i,5:6)+x(i,5:6)*1i];
                b=diag([0 0 0]);
                b(2:3,2:3) = [bs(i,4:5)*1i; bs(i,5:6)*1i];
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
            case 123
                z = [r(i,1:3)+x(i,1:3)*1i; r(i,2)+x(i,2)*1i  r(i,4:5)+x(i,4:5)*1i; r(i,3)+x(i,3)*1i r(i,5:6)+x(i,5:6)*1i];
                b = [bs(i,1:3)*1i; bs(i,2)*1i  bs(i,4:5)*1i; bs(i,3)*1i bs(i,5:6)*1i];
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z)+b/2;     % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(z);         % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(z);         % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(z)+b/2;     % Y22
        end
    end
    if   lin_tipo(i) == 2
        switch phase(i)
            case 1
                z = diag([r(i,1)+x(i,1)*1i vmax vmax]);
                tap=diag([TAPm(i,1)*exp(1i*TAPp(i,1)*pi/180) 1 1]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
            case 2
                z = diag([vmax r(i,4)+x(i,4)*1i vmax]);
                tap=diag([1 TAPm(i,2)*exp(1i*TAPp(i,2)*pi/180) 1]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
            case 3
                z = diag([vmax vmax r(i,6)+x(i,6)*1i]);
                tap=diag([1 1 TAPm(i,3)*exp(1i*TAPp(i,3)*pi/180)]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
            case 12
                z = diag([r(i,1)+x(i,1)*1i r(i,4)+x(i,4)*1i vmax]);
                tap=diag([TAPm(i,[1 2]).*exp(1i.*TAPp(i,[1 2])*pi/180) 1]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
            case 13
                z = diag([r(i,1)+x(i,1)*1i vmax r(i,6)+x(i,6)*1i]);
                tap=diag([1 1 1]);
                tap([1 3],[1 3])=diag(TAPm(i,[1 3]).*exp(1i.*TAPp(i,[1 3])*pi/180));
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
            case 23
                z = diag([vmax r(i,4)+x(i,4)*1i r(i,6)+x(i,6)*1i]);
                tap=diag([ 1 TAPm(i,[2 3]).*exp(1i.*TAPp(i,[2 3])*pi/180)]);
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
            case 123
                z = diag([r(i,1)+x(i,1)*1i r(i,4)+x(i,4)*1i r(i,6)+x(i,6)*1i]);
                tap=diag(TAPm(i,[1 2 3]).*exp(1i.*TAPp(i,[1 2 3])*pi/180));
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(de(i)-1)*3)+inv(z);             % Y11
                Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(de(i)-1)*3,[1 2 3]+(para(i)-1)*3)-inv(tap.*z);        % Y12
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(de(i)-1)*3)-inv(conj(tap).*z);  % Y21
                Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)=Ybus([1 2 3]+(para(i)-1)*3,[1 2 3]+(para(i)-1)*3)+inv(abs(tap).^2.*z);% Y22
                
        end
    end
end





%estrutura para dados do transformador regulador
re.r=reg_r;
re.x=reg_x;
re.g=reg_g;
re.b=reg_b;
re.mph=reg_mph;
re.ph=reg_ph;
re.de= reg_de; %barra de do regulador
re.para= reg_para; %barra para do regulador
re.dec=reg_dec; %barra de - controle
re.parac=reg_parac;  %barra para - controle (barra que está sendo controlada)
re.nb=reg_nb;  %aumento do numero de barras  (barras extras)
re.tap=reg_tap;
re.band= reg_band;
re.V = reg(:,16:18)./reg(:,8);
re.n_lin_reg=n_lin_reg;
re.state=reg(:,23);

end