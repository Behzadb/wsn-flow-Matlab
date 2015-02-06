function Dist_route(numNodes,grid,receiver)
%optimal path flow in routing
%    minimize    q
%    subject to  sigma_(j member of N_i)(r_ij-r_jr)=S_i (for every i member of V)
%                  0<= r_ij <= R_ij (for every i member of V , for every j member of N_i)
%                   T*sigma_(j member of N_i)*E_ij <= B_i (for every i member of V)
% 
%
%   i
%   j
%   B_i :  initial battery energy of node i
%   S_i : the rate at which information is generated at node i
%   r_ij 
%   E_ij : The energy consumption for sending a packet from node i to node j
%   N_i : The set of nodes connected to node i by a link ******----Derived from "E" by E_ij(E_ij(:,'i')==5,:)
%   V : the set of nodes ******----Derived from "E"
%   L : the set of links ******----Derived from "E"
%   R_ij : the maximum flow that a link from node i to node j can support
%   X_ij transmission rate on each link from node i to a neighboring node j, fixed value
%   Tau_ij : the maximum fraction of time for which a link can transmit, fixed
%   c1 
%   c2 : a deterministic path loss model over an AWGN channel, where the received power decays as 1/d^4
%   distance_i 
%	k_iter : iteration k
%   alpha_k : The step-size at iteration k
%	epsilon : an exponentially decreasing
%   (lambda_i,nu_i) : initial point
%   T the rate of information flow from node i to node j
%    

syms i j B_i S_i r_ij E_ij N_i V L R_ij X_ij Tau_ij c1 c2 distance_i k_iter epsilon lambda_i T alpha_k nu_i


%% ============= parameters================
fieldX=500;
fieldY=300;
%Parameters for grid topology
numNodesXY=round(sqrt(numNodes));
step=10;
ini_energy=5; %Joul
k_iter=10


%% =============Main================
R=calc_R(numNodes,fieldX,fieldY)
[netM,RxTxM]=create_netM(numNodesXY,step,grid,fieldX,fieldY);
figure('Color','w','Position',[100 100 800 500]);
E=printNet(R,netM,fieldX,fieldY,receiver);

E_ij=[E(:,1),E(:,2),c1+c2*(E(:,3)).^4];


%% =============Initialize parameters================
i=[1:numNodes];
j=[1:numNodes];

B_i(1:numNodes)=ini_energy;
S_i=poissrnd(numNodes,1,numNodes);
c2=0.1*c1;
X_ij=E;
X_ij(:,3)=0
r_ij=X_ij;
R_ij=X_ij*Tau_ij;

T(r_ij(i,j))=(B_i(i))/symsum(E_ij(i,j)*r_ij(i,j));
T=min(T(r_ij(i,j)));
q=1/T;
pause;%++++++++++++++++++++++++++----------------------------------------




%% ================solve convex problem===========
% cvx_begin 
% 
% variable q
% minimize q
% 
% cvx_end


%-----------------------------Subgradient of -g iteration ----------------

for k=1:k_iter;
    
    alpha_k=max(0.01,(0.5)/sqrt(k_iter));
 
    % ---------------------Eq 7 -------------
    h(k,i)=q(k)*B_i(i)-symsum(E_ij(E_ij(:,1)==i,:)*r(k,i,j));
    f(k,i)=S_i(i)-symsum(r_ij(k,i,j)-r_ij(k,j,i));
    
     % ---------------------Eq 3 -------------
    lambda_i(k+1,i)=lambda_i(k,i)-alpha_k*h(k,i)%plus
    nu_i(k+1,i)=nu_i(k,i)-alpha_k*f(k,i)
    nu_i(k+1,j)=nu_i(k,j)-alpha_k*f(k,j)
    
    q(k)=min(q^2-q*symsum(lambda_i(k,i)*B_i(i)));
    r_ij(k,i,j)=min(epsilon*(r_ij(i,j))^2+r_ij(i,j)*(lambda_i(k,i)*E_ij(i,j)+nu_i(k,i)-nu_i(k,j)));
end;


%% =============Functions================
function R=calc_R(numNodes,fieldX,fieldY)
        p=sqrt((fieldX^2)+(fieldY^2));
        R=p*sqrt(log10(numNodes)/numNodes);
        
function [netM,RxTxM]=create_netM(numNodesXY,step,grid,fieldX,fieldY)
ID=1;
for i=1:numNodesXY
    for j=1:numNodesXY
        netM(1,ID)=ID;% inicializaec topologie
        RxTxM(1,ID)=ID; % inicializace matice RxTxM
        if grid==1
            x=step*j+50;
            y=step*i+50;
        else
            x=rand*fieldX;
            y=rand*fieldY;
        end
        netM(2,ID)=x;
        netM(3,ID)=y;
        RxTxM(2,ID)=0;
        RxTxM(3,ID)=0;
        ID=ID+1;
        
    end
end

% Neigbor matrix creation
function E=printNet(R,netM,fieldX,fieldY,sink)
    set(gca,'FontSize',8,'YGrid','off')
    xlabel('\it x \rm [m] \rightarrow')
    ylabel('\it y \rm [m] \rightarrow')
    
    plot(netM(2,:),netM(3,:),'ko','MarkerSize',5,'MarkerFaceColor','k');
    hold on
    plot(netM(2,sink),netM(3,sink),'h','MarkerSize',15,'MarkerFaceColor','r');
    hold off
    axis([0 fieldX 0 fieldY]);
    hold all;
    radek=1;
    for j=1:numel(netM(1,:))
        for jTemp=1:numel(netM(1,:))
         X1=netM(2,j);
         Y1=netM(3,j);
         X2=netM(2,jTemp);
         Y2=netM(3,jTemp);
         xSide=abs(X2-X1);
         ySide=abs(Y2-Y1);
         d=sqrt(xSide^2+ySide^2);
         if (d<R)&&(j~=jTemp)
             vertice1=[X1,X2];
             vertice2=[Y1,Y2];
             plot(vertice1,vertice2,'-.k','LineWidth',0.1);
             hold all;
             E(radek,1)=j;
             E(radek,2)=jTemp;
             E(radek,3)=d;
             radek=radek+1;
         end
        end
    end
    v=netM(1,:);
    vv=v';
    s=int2str(vv);
    text(netM(2,:)+1,netM(3,:)+3,s,'FontSize',8,'VerticalAlignment','Baseline');
    hold all;

