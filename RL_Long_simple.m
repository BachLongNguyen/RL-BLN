% clear
clc
clear all

load phydatarate_100users_10BSs.mat
N = 5000;                      % Number of interactions                        t = 1..N
u = zeros(S,M);                % Utility payoff of players                     u(j,i)
unoise = zeros(S,M);
a = zeros(N,M);                % Action of players                             a(t,i)
x = zeros(S,M);                % Probability player i choose action j at t     x(j,i)
avgP = zeros(N,M);             % Average payoff of user i at t                 avgP(t,i)
realP = zeros(N,M);            % Real payoff of user i at t                    realP(t,i)
count = zeros(N,S);            % Number of users connecting to j at t          count(t,j)
countweight = zeros(N,S);
avgcount = zeros(N,S);         % Avg number of users connecting to j at t      avgcount(t,j)
no = zeros(S,M);               % Number of times i connected to j untill t     no(j,t,i)
capacity = zeros(S,M);         % Estimated capacity of j by i                  capacity(j,i)
sumCapacity = zeros(S,M);      % Sum estimated capacity of j by i
number = zeros(S,M);           % Estimated no of user on j by i                number(j,i)
sumNumber = zeros(S,M);        % Sum estimated no of user on j by i
sumRegret = zeros(S,S,M);
switching = zeros(N,M);        % Number of switching of user i at t            switching(t,i)
maxSwitchingT = zeros(1,N);    % Max switching number per user at t            maxSwitchingT(1,t)
avgSwitchingT = zeros(1,N);    % Avg switching number per user at t            avgSwitchingT(1,t)
overhead = zeros(N,1);
action = zeros(M,1);           % Number of available action of user i          action(i,1)
% tai moi thoi diem
sumPayoff = zeros(N,1);
fairness = zeros(N,1);
xisquare = zeros(N,M);
noise = 0.3;
mu = zeros(M,1);     
%=====================Vehicular Edge Computing=============================

%=================Simulation Deployment====================================
simTime=150;
D=3000;
nodesX(1)=3000;
nodesY(1)=1000;
nodesX(2)=nodesX(1)+D;
nodesY(2)=1000;
nodesX(3)=nodesX(2)+D;
nodesY(3)=1000;
nodesX0=6000;
nodesY0=2000;
plot(nodesX0,nodesY0,'k^','MarkerSize',20,'MarkerFaceColor','k')
text(nodesX0,nodesY0,'BS','VerticalAlignment','bottom','HorizontalAlignment','left'); 
hold on
grid on

line([1500,10500],[1040,1040])
line([1500,10500],[1030,1030])
line([1500,10500],[1020,1020])
line([1500,10500],[1010,1010])

line([1500,10500],[990,990])
line([1500,10500],[980,980])
line([1500,10500],[970,970])
line([1500,10500],[960,960])    

posRSU=cell(1,length(nodesX));
for i=1:length(posRSU)
   posRSU{1,i}(1,1)=nodesX(1,i);
   posRSU{1,i}(1,2)=nodesY(1,i);
end

text(nodesX(1),nodesY(1),'RSU1','VerticalAlignment','bottom','HorizontalAlignment','left');
text(nodesX(2),nodesY(2),'RSU2','VerticalAlignment','bottom','HorizontalAlignment','left');
text(nodesX(3),nodesY(3),'RSU3','VerticalAlignment','bottom','HorizontalAlignment','left');
plot(nodesX,nodesY,'b^','MarkerSize',20,'MarkerFaceColor','b') 

meanR=600;
R=poissrnd(meanR);

meanr=300;
r=poissrnd(meanr);

viscircles([nodesX(1),nodesY(1)],R,'LineStyle','--','Color','b')    
viscircles([nodesX(2),nodesY(2)],R,'LineStyle','--','Color','b')
viscircles([nodesX(3),nodesY(3)],R,'LineStyle','--','Color','b')
viscircles([nodesX0,nodesY0],4000,'LineStyle','--','Color','k')   

numReq=7;
posX_target=zeros(1,numReq);
posY_target=zeros(1,numReq);
vTarget=zeros(1,numReq);
initY=nodesY(1);
for i=1:numReq
    % requester 1 - task 1 - RSU 1
    posX_target(1,1)=nodesX(1)-100;
    posY_target(1,1)=initY+25;
    vTarget(1,1)=25;
    %---------------------------------------
    % requester 2 - task 2 - RSU 1
    posX_target(1,2)=nodesX(1)+200;
    posY_target(1,2)=initY-15;
    vTarget(1,2)=29;
    %---------------------------------------
    % requester 3 - task 3 - RSU 2
    posX_target(1,3)=nodesX(2)-50;
    posY_target(1,3)=initY+25;
    vTarget(1,3)=33;
    %---------------------------------------
    % requester 4 - task 4 - RSU 2
    posX_target(1,4)=nodesX(2)+70;
    posY_target(1,4)=initY+15;
    vTarget(1,4)=25;
    %---------------------------------------
    % requester 5 -  task 5 - RSU 2
    posX_target(1,5)=nodesX(2)-130;
    posY_target(1,5)=initY-15;
    vTarget(1,5)=29;
    %---------------------------------------
    % requester 6 -  task 6 - RSU 2
    posX_target(1,6)=nodesX(2)+20;
    posY_target(1,6)=initY-25;
    vTarget(1,6)=33;
    %---------------------------------------
    % requester 7 -  task 7 - RSU 3
    posX_target(1,7)=nodesX(3)-10;
    posY_target(1,7)=initY-25;
    vTarget(1,7)=33;
end
posRqs=cell(1,numReq);
for i=1:length(posRqs)
    posRqs{1,i}(1,1)=posX_target(1,i);
    posRqs{1,i}(1,2)=posY_target(1,i);
end 

vTarget=zeros(1,numReq);
for iii=1:numReq
   if posY_target(1,iii)==initY-35||posY_target(1,iii)==initY+35
       vTarget(1,iii)=25; % m/s
   elseif posY_target(1,iii)==initY-25||posY_target(1,iii)==initY+25
       vTarget(1,iii)=28; % m/s
   elseif posY_target(1,iii)==initY-15||posY_target(1,iii)==initY+15
       vTarget(1,iii)=33; % m/s
   end
end
%-------------------------Vehicular Mobility-------------------------------
hr=cell(1,numReq);
% for i=1:simTime    
    for iii=1:numReq
        delete(hr{1,iii})
        hr{1,iii}=plot(posRqs{1,iii}(1,1),posRqs{1,iii}(1,2),'ro','MarkerSize',10,'MarkerFaceColor','r'); 
        text(posRqs{1,iii}(1,1),posRqs{1,iii}(1,2),num2str(iii),'VerticalAlignment','bottom','HorizontalAlignment','left');
        drawnow;
    end

%     for iii=1:numReq
%        if posRqs{1,iii}(1,2)==initY-35||posRqs{1,iii}(1,2)==initY-25||posRqs{1,iii}(1,2)==initY-15
%            posRqs{1,iii}(1,1)=posRqs{1,iii}(1,1)+vTarget(1,iii)*t;
%            posRqs{1,iii}(1,2)=posRqs{1,iii}(1,2);
%        elseif posRqs{1,iii}(1,2)==initY+35||posRqs{1,iii}(1,2)==initY+25||posRqs{1,iii}(1,2)==initY+15
%            posRqs{1,iii}(1,1)=posRqs{1,iii}(1,1)-vTarget(1,iii)*t;
%            posRqs{1,iii}(1,2)=posRqs{1,iii}(1,2);
%        end
%     end 
% end
Rqs=ones(1,length(posRqs));
RqsCurrRSU=cell(1,length(posRqs)); % the vehicle requests which RSU
for i=1:length(posRqs)  
   if Rqs(1,i)==1
       for ii=1:length(posRSU)
           dist=sqrt((posRqs{1,i}(1,1)-posRSU{1,ii}(1,1))^2+(posRqs{1,i}(1,2)-posRSU{1,ii}(1,2))^2);
           if dist<=R
               RqsCurrRSU{1,i}=[i,ii];
           end
       end
   end
end
%--------------------------------------------------------------------------
nRSU=length(posRSU); % a number of RSUs;
nVehicle=length(posRqs); % number of vehicles = number of tasks
nameAVe=cell(1,nVehicle); % name of each action of each requesting vehicle
for i=1:nVehicle
   for ii=1:nRSU
       nameAVe{1,i}{1,ii}=[i,ii];
   end
   nameAVe{1,i}{1,nRSU+1}=[i,0];
end
aVehicle=zeros(1,nVehicle); % a number of actions in nRSUs
probActVehicle=cell(1,nVehicle); % probabiblity that an RSU choose an action
Tsize=zeros(1,nVehicle); % size of each task requested by each vehicle
Tcpu=zeros(1,nVehicle); % (cycles) CPU cycles required to complete each task of each vehicle
RcpuMax=ones(1,nRSU+1); % maximum CPU cycles in each RSU
RcpuT=ones(1,nRSU+1); % (cycles/s) CPU cycles assigned to each task
u=cell(1,N); % utility of each player during N iterations 
%----------------------------Initialization--------------------------------
for i=1:nVehicle
   action(1,i)=nRSU+1; % number of actions of each requesting vehicle
   probActVehicle{1,i}=zeros(1,action(1,i));
   for ii=1:length(probActVehicle{1,i})
      probActVehicle{1,i}(1,ii)=1/length(probActVehicle{1,i}); 
   end
end
%--------------------------------------------------------------------------
nIteration=5000;
M=nVehicle; 
N=nIteration;
S=nRSU+1;
VeActIter=cell(1,N);
TotalCount=cell(1,N);
for it=1:N   % run iterations
   IterCount=zeros(1,S);
   %---------------------Action choosing-----------------------------------
   for i=1:M % the number of vehicles (players)
       temp=rand;
       ii=1;
       k=probActVehicle{1,i}(1,ii);       
       while (temp>k&&ii<S)
           ii=ii+1;
           k=k+probActVehicle{1,i}(1,ii);
       end
       VeActIter{1,it}{1,i}=nameAVe{1,i}{1,ii};  % player 'i' selects action 'ii' at iteration 'it'
       %----Count how many times the players switch their actions----------
       if it>1 
           if VeActIter{1,it}{1,i}~=VeActIter{1,it-1}{1,i}
               switching(it,i)=switching(it-1,i)+1;
           else
               switching(it,i)=switching(it-1,i);
           end
       end
   end
   %--------Count the number of tasks assigned to each RSU-----------------
   for i = 1:S
       for ii=1:length(VeActIter{1,it})
           if i<=S
               if VeActIter{1,it}{1,ii}(1,2)==i
                   IterCount(1,i)=IterCount(1,i)+1;
               end
           else
               if VeActIter{1,it}{1,ii}(1,2)==0
                   IterCount(1,i)=IterCount(1,i)+1;
               end
           end
       end
   end  
   TotalCount{1,it}=IterCount; 
   %----------------------Utility Learning---------------------------------   
%    for i = 1:M
%        if a(it,i) ~= 0
%            %%%%%%%%%% If user i chooses WiFi
%            if a(it,i) <= S/2
%                %%%%%% Real payoff by connecting to WiFi
%                temp = (a(it,:)==a(it,i)).*datarate(a(it,i),:);
%                u(a(it,i),i) = 1/sum(1./temp(temp~=0));
%                %%%%%% Expected payoff if connecting to Wifi
%                for k = 1:S
%                    if (k ~= a(it,i)) && (datarate(k,i) ~= 0)
%                        if k <= S/2
%                            temp = (a(it,:)==k).*datarate(k,:);
%                            u(k,i) = 1/(sum(1./temp(temp~=0)) + 1/datarate(k,i));
%                %%%%%% Expected payoff if connecting to LTE
%                        else
%                            u(k,i) = datarate(k,i)/(count(it,k)+1);
%                        end
%                    end
%                end
%            %%%%%%%%%% If user i chooses LTE
%            else
%                %%%%%% Real payoff by connecting to LTE
%                u(a(it,i),i) = datarate(a(it,i),i)/count(it,a(it,i));
%                %%%%%% Expected payoff if connecting to Wifi
%                for k = 1:S
%                    if (k ~= a(it,i)) %&& (datarate(k,i) ~= 0)
%                        if k <= S/2
%                            temp = (a(it,:)==k).*datarate(k,:);
%                            u(k,i) = 1/(sum(1./temp(temp~=0)) + 1/datarate(k,i));                           
%                %%%%%% Expected payoff if connecting to LTE
%                        else
%                            u(k,i) = datarate(k,i)/(count(it,k)+1);
%                        end
%                    end
%                end
%            end
%            realP(it,i) = u(a(it,i),i);       % Real payoff of user i at t
%            unoise(:,i) = abs(normrnd(u(:,i),noise.*u(:,i)));
%            if it == 1
%                avgP(it,i) = realP(it,i);
%            else
%                avgP(it,i) = (avgP(it-1,i)*(it-1)+realP(it,i))/it;
%            end 
%            xisquare(it,i) = realP(it,i)^2;
%            mu(i,1) = 2*max(unoise(:,i))*(action(i,1)-1)+1; % converge nhanh hay cham do mu
%        end
%    end
   for i = 1:M % the number of players (vehicles)
       if isempty(VeActIter{1,it}{1,i})
       else
           tempVe=VeActIter{1,it}{1,i}(1,1);
           tempRSU=VeActIter{1,it}{1,i}(1,2);
           if tempRSU==RqsCurrRSU{1,tempVe}(1,2)
               Rx=posRSU{1,tempRSU}(1,1);
               Ry=posRSU{1,tempRSU}(1,2);
               Vx=posRqs{1,tempVe}(1,1);
               Vy=posRqs{1,tempVe}(1,2);
               tempRate=DataRate(Rx,Ry,Vx,Vy);
               tempSumCpu=Tcpu(1,tempVe);
               for ii=1:M
                  if ii~=i
                     if VeActIter{1,it}{1,ii}(1,2)==tempRSU
                         tempSumCpu=tempSumCpu+...
                             Tcpu(1,VeActIter{1,it}{1,ii}(1,1));
                     end
                  end
               end
               % Satisfy the constraints of maximum CPU cycles
               if tempSumCpu>RcpuMax(1,tempRSU)
                   delay=1000000;
                   cost=1000000;
               else
                   delay=Tsize(1,tempVe)/tempRate...
                         +Tcpu(1,tempVe)/RcpuT(1,tempRSU);
                   cost=bwV2I*cCom+Tcpu(1,tempVe)*cComp;
               end      
               % real utility
               u{1,it}{1,i}(1,tempRSU)=-(delay+cost);
               % expected utility
               for ii=1:S
                  if ii~=tempRSU
                      if ii==S
                          Rx=posRSU{1,ii}(1,1);
                          Ry=posRSU{1,ii}(1,2);
                          nHop=1;
                          tempSumCpu=Tcpu(1,tempVe);
                          for iii=1:M
                              if iii~=i
                                 if VeActIter{1,it}{1,iii}(1,2)==ii
                                     tempSumCpu=tempSumCpu+...
                                         Tcpu(1,VeActIter{1,it}{1,iii}(1,1));
                                 end
                              end
                          end
                          % Satisfy the constraints of maximum CPU cycles
                          if tempSumCpu>RcpuMax(1,ii)
                               delay=1000000;
                               cost=1000000;
                          else
                               delay=Tsize(1,tempVe)/tempRate...
                                    +(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                                    +Tcpu(1,tempVe)/RcpuT(1,nRSU+1);  
                               cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                                    +Tcpu(1,tempVe)*cComp0;
                          end
                      else
                          Rx=posRSU{1,ii}(1,1);
                          Ry=posRSU{1,ii}(1,2);
                          nHop=floor(abs(posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,1)...
                                -posRSU{1,ii}(1,1))/D);
                          tempSumCpu=Tcpu(1,tempVe);
                          for iii=1:M
                              if iii~=i
                                 if VeActIter{1,it}{1,iii}(1,2)==ii
                                     tempSumCpu=tempSumCpu+...
                                         Tcpu(1,VeActIter{1,it}{1,iii}(1,1));
                                 end
                              end
                          end
                          % Satisfy the constraints of maximum CPU cycles
                          if tempSumCpu>RcpuMax(1,ii)
                               delay=1000000;
                               cost=1000000;
                          else
                               delay=Tsize(1,tempVe)/tempRate...
                                    +(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                                    +Tcpu(1,tempVe)/RcpuT(1,nRSU+1);  
                               cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                                    +Tcpu(1,tempVe)*cComp;
                          end
                      end
                      u{1,it}{1,i}(1,ii)=-(delay+cost);
                  end
               end
           else
               if tempRSU==0
                   Rx=posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,1);
                   Ry=posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,2);
                   Vx=posRqs{1,tempVe}(1,1);
                   Vy=posRqs{1,tempVe}(1,2);
                   tempRate=DataRate(Rx,Ry,Vx,Vy);
                   nHop=1;
                   tempSumCpu=Tcpu(1,tempVe);
                   for ii=1:M
                        if ii~=i
                            if VeActIter{1,it}{1,ii}(1,2)==tempRSU
                                tempSumCpu=tempSumCpu+...
                                Tcpu(1,VeActIter{1,it}{1,ii}(1,1));
                            end
                        end
                   end
                   % Satisfy the constraints of maximum CPU cycles
                   if tempSumCpu>RcpuMax(1,tempRSU)
                        delay=1000000;
                        cost=1000000;
                   else
                        delay=Tsize(1,tempVe)/tempRate...
                            +(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                            +Tcpu(1,tempVe)/RcpuT(1,tempRSU);
                        cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                            +Tcpu(1,tempVe)*cComp0;
                   end      
                   % real utility
                   u{1,it}{1,i}(1,tempRSU)=-(delay+cost);                 
                   % expected utility
                   for ii=1:S
                      if ii~=S
                          Rx=posRSU{1,ii}(1,1);
                          Ry=posRSU{1,ii}(1,2);
                          nHop=floor(abs(posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,1)...
                                -posRSU{1,ii}(1,1))/D);
                          tempSumCpu=Tcpu(1,tempVe);
                          for iii=1:M
                              if iii~=i
                                 if VeActIter{1,it}{1,iii}(1,2)==ii
                                     tempSumCpu=tempSumCpu+...
                                         Tcpu(1,VeActIter{1,it}{1,iii}(1,1));
                                 end
                              end
                          end
                          % Satisfy the constraints of maximum CPU cycles
                          if tempSumCpu>RcpuMax(1,ii)
                               delay=1000000;
                               cost=1000000;
                          else
                               delay=Tsize(1,tempVe)/tempRate...
                                    +(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                                    +Tcpu(1,tempVe)/RcpuT(1,nRSU+1);  
                               cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                                    +Tcpu(1,tempVe)*cComp;
                          end
                          u{1,it}{1,i}(1,ii)=-(delay+cost);
                      end
                   end
               else
                   Rx=posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,1);
                   Ry=posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,2);
                   Vx=posRqs{1,tempVe}(1,1);
                   Vy=posRqs{1,tempVe}(1,2);
                   tempRate=DataRate(Rx,Ry,Vx,Vy);
                   nHop=floor(abs(posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,1)...
                       -posRSU{1,tempRSU}(1,1))/D);
                   tempSumCpu=Tcpu(1,tempVe);
                   for ii=1:M
                        if ii~=i
                            if VeActIter{1,it}{1,ii}(1,2)==tempRSU
                                tempSumCpu=tempSumCpu+...
                                Tcpu(1,VeActIter{1,it}{1,ii}(1,1));
                            end
                        end
                   end
                   % Satisfy the constraints of maximum CPU cycles
                   if tempSumCpu>RcpuMax(1,tempRSU)
                        delay=1000000;
                        cost=1000000;
                   else
                        delay=Tsize(1,tempVe)/tempRate+(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                            +Tcpu(1,tempVe)/RcpuT(1,tempRSU);
                        cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                            +Tcpu(1,tempVe)*cComp;
                   end      
                   % real utility of player 'i'
                   u{1,it}(1,i)=-(delay+cost);
                   % expected utility
                   for ii=1:S
                      if ii~=S
                          nHop=1;
                          tempSumCpu=Tcpu(1,tempVe);
                          for iii=1:M
                              if iii~=i
                                 if VeActIter{1,it}{1,iii}(1,2)==ii
                                     tempSumCpu=tempSumCpu+...
                                         Tcpu(1,VeActIter{1,it}{1,iii}(1,1));
                                 end
                              end
                          end
                          % Satisfy the constraints of maximum CPU cycles
                          if tempSumCpu>RcpuMax(1,ii)
                               delay=1000000;
                               cost=1000000;
                          else
                               delay=Tsize(1,tempVe)/tempRate...
                                    +(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                                    +Tcpu(1,tempVe)/RcpuT(1,nRSU+1);  
                               cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                                    +Tcpu(1,tempVe)*cComp0;
                          end
                          u{1,it}{1,i}(1,ii)=-(delay+cost);
                      elseif ii~=tempRSU
                          Rx=posRSU{1,ii}(1,1);
                          Ry=posRSU{1,ii}(1,2);
                          nHop=floor(abs(posRSU{1,RqsCurrRSU{1,tempVe}(1,2)}(1,1)...
                                -posRSU{1,ii}(1,1))/D);
                          tempSumCpu=Tcpu(1,tempVe);
                          for iii=1:M
                              if iii~=i
                                 if VeActIter{1,it}{1,iii}(1,2)==ii
                                     tempSumCpu=tempSumCpu+...
                                         Tcpu(1,VeActIter{1,it}{1,iii}(1,1));
                                 end
                              end
                          end
                          % Satisfy the constraints of maximum CPU cycles
                          if tempSumCpu>RcpuMax(1,ii)
                               delay=1000000;
                               cost=1000000;
                          else
                               delay=Tsize(1,tempVe)/tempRate...
                                    +(Tsize(1,tempVe)/bwR+2*alpha*nHop)...
                                    +Tcpu(1,tempVe)/RcpuT(1,nRSU+1);  
                               cost=bwV2I*cCom+Tsize(1,tempVe)*cMig...
                                    +Tcpu(1,tempVe)*cComp;
                          end
                          u{1,it}{1,i}(1,ii)=-(delay+cost);
                      end
                   end                   
               end
           end
       end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%% Regret Matching %%%%%%%%%%%%%%%%%%%%%%%%%% 
   for i = 1:M
       if a(it,i) ~= 0
           x(a(it,i),i) = 0;
           for k = 1:S
               if (k ~= a(it,i)) && (datarate(k,i) ~= 0)                   
                   sumRegret(k,a(it,i),i) = sumRegret(k,a(it,i),i) + unoise(k,i) - unoise(a(it,i),i);
                   x(k,i) = 1/mu(i,1)*max(1/it*sumRegret(k,a(it,i),i),0); % check ky
               end
           end
           x(a(it,i),i) = 1 - sum(x(:,i),1);
       end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%% Fairness -- PoA %%%%%%%%%%%%%%%%%%%%%%%%%%   
   fairness(it,1) = (sum(realP(it,1:M),2))^2/sum(xisquare(it,1:M),2)/M;
   sumPayoff(it,1) = sum(realP(it,1:M),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main--algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
% plot(1:N,avgP(:,:),'-','LineWidth',3)
% ylabel('Average payoff')
% xlabel('Iteration')
% 
% figure(2)
% hold on
% plot(1:N,realP(:,:),'-','LineWidth',3)
% ylabel('Real payoff')
% xlabel('Iteration')
% 
% figure(3)
% hold on
% plot(1:N,count(:,:),'-','LineWidth',3)
% ylabel('Number of users connecting to each networks')
% xlabel('Iteration')
% axis([0 N 0 25])
% 
% figure(4)
% hold on
% plot(1:N,avgcount(:,:),'-','LineWidth',3)
% ylabel('Average number of users connecting to each networks')
% xlabel('Iteration')
% axis([0 N 0 25])
% 
% figure(5)
% hold on
% plot(1:N,sumPayoff(:,1),'-x')
% ylabel('Total payoffs of all users')
% xlabel('Iteration')
% legend('sumPayoff')
% 
% figure(6)
% hold on
% plot(1:N,fairness(:,1),'-x')
% ylabel('System Fairness Index')
% xlabel('Iteration')
% legend('fairness')