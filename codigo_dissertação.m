clear
clc
load('gab_anti_horario.mat') %arquivo que contem os dados das coordenadas e dos elementos de fronteira
format short

%%Construcao da parabola do perfil de velocidade na entrada do Rio das Mortes
%
%%%----------------------Determinar a velocidade media----------------------%%% 
%
L1=230; %largura do rio onde foram feitas as medicoes
a1=-(1.1413/((L1/2)^2)); %coeficiente da f1(x)
syms x
f1(x)=a1*x*(x-L1); %funcao que determina o perfil de velocidade onde foram feitas as medicoes
Vm1=(1/L1)*(int (f1(x),x,0,230)); %calcular a area sob a curva para determinar a velocidade media
%
%%%----------------------Determinar a velocidade media----------------------%%%
%L2 = sqrt((coordenadas(1,1)-coordenadas(693,1))^2 + (coordenadas(1,2)-coordenadas(693,2))^2); 
L2=norm((coordenadas(1,:)-coordenadas(693,:)),2); %largura do local onde estamos trabalhando
c=L1/L2; %obtido da equacao V=Vazao/area
Vm2=Vm1*c; %velocidade media no local onde estamos trabalhando
If2=-(4/L2^3)*(int(x*(x-L2),x,0,L2)); %calcular a area sob a funcao f2(x)
q=Vm2/If2 ; %velocidade maxima do rio no local onde estamos trabalhando
b=-((4*q)/(L2^2)); %coeficiente da f2(x)
f2(x)=b*x*(x-L2); %funcao que indica o perfil de velocidade do local onde estamos trabalhando
%
%% Tranformacao Linear
%%%----------------------Dados de entrada do dominio----------------------%%%
%
x1=coordenadas(1,1); %coordenadasenada do no 1 no eixo x
y1=coordenadas(1,2); %coordenadasenada do no 1 no eixo y
x693=coordenadas(693,1); %coordenadasenada do no 623 no eixo x
y693=coordenadas(693,2); %coordenadasenada do no 623 no eixo y
%%%----------------------Determinando o x-vertice e o y-vertice----------------------
%
m=(y1-y693)/(x1-x693); %coeficiente da reta que passa pelas coordenadasenadas da entrada do dominio
n=((y693-y1)/(x1-x693))*x1 + y1; %coeficiente da reta que passa pelas coordenadasenadas da entrada do dominio 
c1=(x1+x693)/2; %media do comprimento do eixo x
c2=(y1+y693)/2; %%media do comprimento do eixo y
xv=(1/(m^2+1))*(c1+c1*m^2-m*q*sqrt(m^2+1)); %xv do vetor que indica velocidade maxima do rio
xv=double(xv); %determina o valor do xv no formato decimal
yv=(-xv/m)+c2+(c1/m); %%xv do vetor que indica velocidade maxima do rio
yv=double(yv); %determina o valor do yv no formato decimal
nos=1:693;
nos=[nos 1];
%%%----------------------Tranformacao Linear----------------------%%%
%
front693=find(elementos_fronteira(:,1)==693); %procura na primeira coluna do elementos_fronteira todos os elementos que tem seus vertices sob a linha 693  
front693=elementos_fronteira(front693,2:3); %determina os nos dos elementos encontrados
front693=[front693(:,1);front693(end,2)]; %retira os nos repetidos
%
for i=1:length(front693)
    D(i)=norm(coordenadas(front693(i),:)-coordenadas(front693(1),:),2); %distancia de cada no ate o no inicial 
end
x=D; 
y=zeros(size(x)); 
for i=1:length(x)
    [X,Y]=transformacao(x(i),y(i),x1,x693,xv,y1,y693,yv,L2,q); %transformacao linear do eixo x
    ptos_eixox(i,:)=[X Y]; %associa os pontos do eixo x a transformacao linear
    p(i)=b*x(i)*(x(i)-L2); %parabola do perfil de velocidade com i variando nos nos do perfil de entrada
    [X,Y]=transformacao(x(i),p(i),x1,x693,xv,y1,y693,yv,L2,q); %transformacao linear da parabola
    ptos_parabola(i,:)=[X Y]; %associa os pontos da parabola a transformacao linear
end
% figure(2)
% plot ([x693;x1],[y693;y1]) %visualizacao dos pontos extremos
% %plot (xv,yv,'r*')
% %plot(coordenadas(nos,1),coordenadas(nos,2),xv,yv,'m*') %visualizacao do desenho do rio e do ponto que ira definir o xv e yv
% view([0 90]) %visualizar o grafico completo
% axis equal %mesmas unidades de comprimento no grafico
% hold on %gera um grafico em cima do outro
% quiver(ptos_eixox(:,1),ptos_eixox(:,2),ptos_parabola(:,1)-ptos_eixox(:,1),ptos_parabola(: ,2)-ptos_eixox(: ,2),'c') %visualizacao dos vetores, lembrar que o 0 plota o tamanho real do vetor
% axis equal 
vel_entrada=double([ptos_eixox(:,1) ptos_eixox(:,2) ptos_parabola(:,1)-ptos_eixox(:,1) ptos_parabola(: ,2)-ptos_eixox(: ,2)]);
save vel_entrada vel_entrada
function [X,Y]= transformacao(x,y,x1,x693,xv,y1,y693,yv,L2,q)
T1 = [((x1-x693)/L2) ((1/(2*q))*(2*xv-x1-x693));((y1-y693)/L2) ((1/(2*q))*(2*yv-y1-y693))];
T2 = [x693;y693];
a = T1*[x;y]+T2; %matriz da transformacao linear
X = a (1) ; %coordenadasenada do X
Y = a (2) ;  %coordenadasenada do Y


%%%% PLOTAGEM %%%%%

load('gab_anti_horario.mat')
load('vel_entrada.mat')
load('CENARIO1_10.mat')
load('nos_contorno.mat')
%%
%plot do campo de velocidade:
quiver(coordenadas(:,1),coordenadas(:,2),solucao(:,end-1),solucao(:,end))
hold on
plot(coordenadas(nos_rio,1),coordenadas(nos_rio,2),"red")
plot(coordenadas(nos_ilha,1),coordenadas(nos_ilha,2),"red")
axis equal

%%
%plot do perfil de poluicao
trisurf(elementos,coordenadas(:,1),coordenadas(:,2),solupol(:,end),'edgeColor','none','faceColor','interp')
axis equal
view([0 90])

%% plot do perfil de poluicao com varias imagens

nti=size(solupol,2);
p=floor(1+nti/10);
if mod(p,2)==0
    p=p/2;
else
    p=ceil(p/2);
end
k=1;
t=tiledlayout(p,2,'TileSpacing','Compact');
for i=1:10:nti
    ax1=nexttile;
    %ax1=subplot(p,2,k);
    trisurf(elementos,coordenadas(:,1),coordenadas(:,2),solupol(:,i),'edgeColor','none','faceColor','interp')
    axis equal
    view([0 90])
    str=strcat('Iteracao:',string(iteracao(i)));
    title(str)
    ax1.FontSize = 8;
    k=k+1;
end
t.TileSpacing='compact';
t.Padding='compact';

%%
str ='1';
str={str,'862'};

for i=3:size(elementos_fronteira,1)
    str{1,i}=mat2str(elementos_fronteira(i,2));
end

plot(coordenadas(elementos_fronteira(:,2),1),coordenadas(elementos_fronteira(:,2),2),'.')
text(coordenadas(elementos_fronteira(:,2),1),coordenadas(elementos_fronteira(:,2),2),str)


%%
A=zeros(size(elementos,1),1);

for i=1:size(elementos,1)
    a(1)=norm(coordenadas(elementos(i,1),:)-coordenadas(elementos(i,2),:),2);
    a(2)=norm(coordenadas(elementos(i,1),:)-coordenadas(elementos(i,3),:),2);
    a(3)=norm(coordenadas(elementos(i,2),:)-coordenadas(elementos(i,3),:),2);

    v(1)=norm([solucao(elementos(i,1),end-1) solucao(elementos(i,2),end)],2);
    v(2)=norm([solucao(elementos(i,1),end-1) solucao(elementos(i,3),end)],2);
    v(3)=norm([solucao(elementos(i,2),end-1) solucao(elementos(i,3),end)],2);
    A(i,1)=min(a)/max(a);
    A(i,2)=mean(a)*mean(v)/0.20;

end
figure
qm=A(:,1);
qm=sort(qm);
plot(qm)

figure
pe=A(:,2);
pe=sort(pe);
plot(pe)


%%%% PROGRAMA OTIMIZADO %%%%%%%%%

% arquivo malha%
load('gab_anti_horario.mat')
load('nos_contorno.mat')
load('vel_entrada.mat')
ntn=size(coordenadas,1);
%
continuacao=input('Continuacao de outra simulacao? (responda 1 para sim ou 0 para nao) = ');
n_iteracoes=input('Numero de iteracoes = ');

if continuacao==1
    nome_arquivo=input('De o nome do arquivo para continuacao: ');
    load(nome_arquivo)
    v1=solucao(:,size(solucao,2)-1:size(solucao,2));
    pol1=solupol(:,size(solupol,2));
else
    solucao=[];
    iteracao=[];
    solupol=[];
    nome_arquivo=strcat('CENARIOteste_','10');
    save(nome_arquivo,'solucao','iteracao','solupol','-v7.3')
    v1=zeros(ntn,2);
    pol1=zeros(ntn,1);
end
proximo=nome_arquivo(end-1:end);
proximo=str2num(proximo);
proximo=proximo+1;
proximo=num2str(proximo);
%%
phii_phij=[1/12 1/24 1/24; 1/24   1/12   1/24; 1/24  1/24   1/12];
dfipsi_dfipsi=[1/2  0  -1/2;0 0 0; -1/2  0  1/2];
dfijpsi_dfieta=[0  0  0;1/2 0  -1/2;-1/2  0  1/2];% dfcsi_dfieta
dfieta_dfipsi=[0  1/2 -1/2;0 0  0;0  -1/2  1/2];%dpjeta_dficsi
dfijeta_dfieta=[0 0 0;0  1/2  -1/2;0  -1/2  1/2];
%
dphipsij_phii=[1/6  0  -1/6;1/6  0  -1/6;1/6  0  -1/6]; %dphi/dpsi * phii
dphijeta_phii=[ 0  1/6  -1/6;0  1/6  -1/6;0  1/6  -1/6];
%
elem_front1=[1/3 1/6;1/6 1/3]; %phii_phij
elem_front2=[1/2;1/2];% um phi

%Dados de entrada da Matriz%
nte=size(elementos,1);
ntn=size(coordenadas,1);
ntef=size(elementos_fronteira,1);
nf=elementos_fronteira(:,2:3);
nf=nf(:);
nf=sort(nf);
nf=unique(nf);
nf=nf(end);

teta1=24e-12;
teta2=24e-2;
teta3=24e-6;
teta4=24e-6;
 
kk5=24e-5;
kk6=0;
kk7=24e-5;
kk8=24e-5;
kk9=0;
kk10=24e-5;
kk11=24e-5;

band=zeros(ntef,1);

band(13124:13200)=1; %G1
band(12250:12779)=2; %G2
band(6037:6040)=3;   %G3
band(5052:5057)=4;   %G4
band(6796:6854)=5;   %G5
band(1:5051)=6;      %G6
band(5058:6036)=7;   %G7
band(6041:6795)=8;   %G8
band(6855:12249)=9;  %G9
band(12780:13123)=10; %G10
band(13201:14724)=11; %G11

%Parametros fundamentais%
dt=0.05;
Re=100; % numero de reynolds
alpha=0.2; % coeficiente de difusao
sigma=24e-10;% decaimento

%Dados de entrada%

fronteira4=nos_mortes;

%co_pre=fronteira2;
inco_pre=(1:ntn)';
%inco_pre(co_pre)=[];

%
inco_v1=[nf+1:ntn]; %co_pre'];% estou tirando a fronteira de entrada,as laterais e as do disco.
inco_v1=sort(inco_v1);
%inco_v1(cacau)=[];
inco_v1=inco_v1';
%
co_v1=(1:ntn)';
co_v1(inco_v1)=[];

% % cognitas e incognitas%
inco_v2=inco_v1;
co_v2=co_v1;

if isempty(iteracao)==1
    it_inicial=1;
else
    it_inicial=iteracao(end)+1;
end



Mfront=zeros(2,3,nf);
MFIJ=zeros(2,6,nf);
parfor k=1:nf
    MFIJ(:,:,k)=EIFFL_FRONT(elementos_fronteira(k,:));
    Mfront(:,:,k)=trab_matriz_front(k,coordenadas,elementos_fronteira,elem_front1,band,...
        teta1,teta2,teta3,teta4,kk5,kk6,kk7,kk8,kk9,kk10,kk11)
end
% I e J para a fronteira phii_phij robin%
I_front_1=MFIJ(:,1:2,:);  I_front_1=I_front_1(:);
J_front_1=MFIJ(:,3:4,:);  J_front_1=J_front_1(:);

% I e J para a fronteira um_phii Von Neumann%
I_front_2=MFIJ(:,5,:);    I_front_2=I_front_2(:);
J_front_2=MFIJ(:,6,:);    J_front_2=I_front_2(:);

% valores para as condicoes de fronteiras

Vfront1=Mfront(:,1:2,:);  Vfront1=Vfront1(:);
Vfront2=Mfront(:,3,:);    Vfront2=Vfront2(:); %alterei aqui

%Montando as matrizes sparse para as fronteiras%

R=sparse(I_front_1,J_front_1,Vfront1,ntn,ntn); %Sparse para front phi_phij
S=sparse(I_front_2,1,Vfront2,ntn,1); %sparse para front phi
parfor ke=1:nte
    %      COORD=[coordenadas(elementos(ke,1),:);   %coordenadas(elementos(ke,1),2);
    %             coordenadas(elementos(ke,2),:);   %coordenadas(elementos(ke,2),2);
    %             coordenadas(elementos(ke,3),:)];   %coordenadas(elementos(ke,3),2)]
    MMM(:,:,ke)=EIFFL_II_JJ(elementos(ke,:));%armazenei as phi_phij
end
% valores para I e J empilhados%
I=MMM(:,1:3,:); I=I(:);
J=MMM(:,4:6,:); J=J(:);

for it=it_inicial:n_iteracoes
    tic
    if it<=it_inicial+0
        % theta variando%
        MM=zeros(3,12,nte);
        MMM=zeros(3,6,nte);
        WW=zeros(3,3,nte);
      
        
        parfor ke=1:nte
            %      COORD=[coordenadas(elementos(ke,1),:);   %coordenadas(elementos(ke,1),2);
            %             coordenadas(elementos(ke,2),:);   %coordenadas(elementos(ke,2),2);
            %             coordenadas(elementos(ke,3),:)];   %coordenadas(elementos(ke,3),2)]
            
            
            MM(:,:,ke)=MATMP(ke,elementos,coordenadas,phii_phij,dfipsi_dfipsi,...
                dfijpsi_dfieta,dfieta_dfipsi,dfijeta_dfieta,dphipsij_phii,dphijeta_phii);
            
            WW(:,:,ke)=vel_media(ke,elementos,coordenadas,v1,dphipsij_phii,dphijeta_phii); % subrotina
        end
        % valores para I e J empilhados%
        %Valores para as matrizes M e PQ%
        V1=MM(:,1:3,:);   V1=V1(:); % phi_phj
        V2=MM(:,4:6,:);   V2=V2(:); % Gradiente phi_phj
        V3=MM(:,7:9,:);   V3=V3(:); % dphix_phi
        V4=MM(:,10:12,:); V4=V4(:); % dphiy_phi
        %Sparse para a matriz M fruto de phii_phij_global%
        M=sparse(I,J, V1,ntn,ntn); % matriz esparsa para phi_phij
        PQ=sparse(I,J,V2,ntn,ntn); % matriz esparsa para o gradiente
        M1=sparse(I,J,V3,ntn,ntn);% matriz esparsa para dphix_phi
        M2=sparse(I,J,V4,ntn,ntn);% matriz esparsa para dphiy_phi
        
        WW=WW(:);% empilhamento dos dados para a velocidade media nos tres nos
        %length(WW)teste
        FF=sparse(I,J,WW,ntn,ntn); %esparsa para a velocidade media nos tres nos
        % montando o b1_esterisco
        
        b1_star=(1/dt)*M(inco_v1,:)*v1(:,1)-((1/dt)*M(inco_v1,co_v1)+(1/Re)*PQ(inco_v1,co_v1)+ ...
            FF(inco_v1,co_v1))*v1(co_v1,1); %comentado geral
        %montando o b2_asterisco
        b2_star=(1/dt)*M(inco_v2,:)*v1(:,2)-((1/dt)*M(inco_v2,co_v2)+(1/Re)*PQ(inco_v2,co_v2)+ ...
            FF(inco_v2,co_v2))*v1(co_v2,2); %comentado geral
        
        %  %comentado geral
        NS_apr1=(1/dt)*M(inco_v1,inco_v1)+(1/Re)*PQ(inco_v1,inco_v1)+FF(inco_v1,inco_v1);
        NS_apr2=(1/dt)*M(inco_v2,inco_v2)+(1/Re)*PQ(inco_v2,inco_v2)+FF(inco_v2,inco_v2);
        
        % Resolvendo os sistemas para v1 e v2 aproximado
        VE1=NS_apr1\b1_star;
        VE2=NS_apr2\b2_star;
        
        %
        vel_star=zeros(ntn,2);
        %
        vel_star(inco_v1,1)=VE1; % abrange toda superficie e os no da fronteira cog_pressao
        vel_star(inco_v2,2)=VE2; % abrange toda superficie e os no da fronteira cog_pressao
        %     vel_star(cacau,1)=0.432e+02;% km/dia
        %
        vel_star(fronteira4,:)=rio_mortes;
        %vel_star(fronteira4,2)=vel_entrada(:,3);
        %corregos%
        vel_star(nos_antartico,:)=rio_antartico;
        vel_star(nos_juma,:)=rio_juma;
        % vel_star(belavista,1)=0.43e+2;
        % vel_star(capivara,1)=-0.43e+2;
        % vel_star(stereza,1)=-0.43e+2;
        %
        bpre=-(1/dt)*M1(inco_pre,:)*vel_star(:,1)-(1/dt)*M2(inco_pre,:)*vel_star(:,2);
        wp=PQ(inco_pre,inco_pre)\bpre;
        pre=zeros(ntn,1);
        pre(inco_pre)=wp; %calculo da pressao
        
        %
        b1=M(inco_v1,:)*vel_star(:,1)-dt*M1(inco_v1,:)*pre;%comentado geral
        b2=M(inco_v2,:)*vel_star(:,2)-dt*M2(inco_v2,:)*pre;
        
        %resolver o sistema
        
        velocidade1=M(inco_v1,inco_v1)\b1;
        velocidade2=M(inco_v2,inco_v2)\b2;
        % Armazenando na matriz velocity12
        v2=zeros(ntn,2);
        
        %v2(fronteira4,1)=vel_rio;%velocidade na entrada do duto
        v2(inco_v1,1)=velocidade1; % abrange toda superficie e os no da fronteira cog_pressao
        v2(inco_v2,2)=velocidade2;
        
        v2(fronteira4,:)=rio_mortes;
        %v2(fronteira4,1)=ptos_parabola(:,1)-ptos_eixox(:,1);
        %Corregos corrigidas%
        v2(nos_antartico,:)=rio_antartico;
        v2(nos_juma,:)=rio_juma;
        % v2(cacau,1)=-0.43e+2;
        % v2(bacuri,1)=-0.43e+2;
        % v2(belavista,1)=0.43e+2;
        % v2(capivara,1)=-0.43e+2;
        % v2(stereza,1)=-0.43e+2;
        %
        
        % km/dia
        %   v2(cacau,1)=43.2;
        % v2(fronteira4,2)=0.00923;%*dt*it;
        % v2(fronteira4,1)=0.0025;%*dt*it;
        % velocidade no dois passos%
        
    else
        v2=v1;
    end
    VPOL=zeros(3,3,nte);
    parfor el=1:nte
        VPOL(:,:,el)=vel_grad_phii(el,elementos,v1,v2,coordenadas,dphipsij_phii,dphijeta_phii);
    end
    VPOL=VPOL(:);
    G=sparse(I,J,VPOL,ntn,ntn);
    %Equacao da Poluicao%
    GE=M+(dt/2)*(alpha*PQ+sigma*M+R)+(dt/2)*G;
    EL=(M-(dt/2)*(alpha*PQ+sigma*M+R)-(dt/2)*G)*pol1+dt*S;%dt*f
    %
    pol2=GE\EL;
    v1=v2;
    pol1=pol2;
    % tetaco1=tetaco;
    % tetaem1=tetaem;
    %
    it
    
    if mod(it,10)==1
        if size(solupol,2)<50
            solucao=[solucao v2];
            iteracao=[iteracao;it];
            solupol=[solupol  pol2];
            save(nome_arquivo,'solucao','iteracao','solupol','-v7.3')
        else
            solucao=v2;
            iteracao=[iteracao;it];
            solupol=pol2;
            nome_arquivo=strcat('CENARIOteste_',proximo);
            save(nome_arquivo,'solucao','iteracao','solupol','-v7.3')
            proximo=nome_arquivo(end-1:end);
            proximo=str2num(proximo);
            proximo=proximo+1;
            proximo=num2str(proximo);
        end
    end
    toc