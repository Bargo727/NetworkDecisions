clear;
%This program will compute probabilities of what fraction of agents make
%the correct choice after first wave and after equilibration.
NNN = 5:190;
L = length(NNN);

%number of trials
nn = 5; 

A1 = zeros(L,1);
A2 = zeros(L,1);

%Initializing the length of time for which the simulation will run
l = 200000;
t = 0:l;

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
p = exp(.05)/10;  %Adjust these parameters to alter SNR.  
q = 1/10;
s = 1-p-q;

evidence_amounts = [log(p/q), -log(p/q), 0];

%Threshold values
thetap = 10; 
thetam = -10;

for jjj = 1:L
    jjj
%Here we will compute some general statistics for the general clique.
%Number of agents in clique
N = NNN(jjj);
for ii = 1:nn


%We draw the observation of the environment at each time from a true
%multinomial distribution.  Each agent's is independent.
LLR = zeros(l,N);
for jj = 1:N
    obs = mnrnd(1,[p,q,s],l);
    LLR(:,jj) = obs*(evidence_amounts)';
end

EV = cumsum(LLR);
[agent_right, tt_right] = find(EV'>= thetap,1);
[agent_wrong, tt_wrong] = find(EV'<= thetam,1);

TT = [tt_right,tt_wrong];
AA = [agent_right,agent_wrong];

tt = min(TT);
agent_index = find(TT == tt,1);
agent = AA(agent_index);
a1 = 0;
d1 = 0;

if(EV(tt,agent) > 0)
    theta = thetap;
    a1 = 1;
else
    theta = thetam;
    d1 = 1;
end

a = 0;
d = 0;

for jj = 1:N
    if(jj ~= agent)
        EV(tt,jj) = EV(tt,jj) + theta;
        a = a + (EV(tt,jj) >= thetap);
        d = d + (EV(tt,jj) <= thetam);
    end
end

A1(jjj) = A1(jjj)+ (a+a1);

u1 = N - 1 - a;

%Indexing Guide for the small threshold example
    %-2 -1 0 1 2 3
    %     0 1 2 3
    
total_states = thetap+abs(thetam) + 1 - 2; %Include 0 in the states but not the boundary values
stateTrackM = zeros(l,total_states);
stateTrackM(1,:) = zeros(1,total_states);
stateTrackM(1,abs(thetam)) = 1;

stateTrackP = zeros(l,total_states);
stateTrackP(1,:) = zeros(1,total_states);
stateTrackP(1,abs(thetam)) = 1;

% negObsProbP = q;
% posObsProbP = p;
% staticObsProb = 1 - p - q;
% 
% negObsProbM = p;
% posObsProbM = q;

survivalP = zeros(1,l);
survivalM = zeros(1,l);
survivalP(1) = 1;
survivalM(1) = 1;
positiveP = zeros(1,l);
positiveM = zeros(1,l);
positiveP(1) = 0;
positiveM(1) = 0;

for j = 2:l
    stateTrackP(j,:) = q*[stateTrackP(j-1,2:end),0] + p*[0,stateTrackP(j-1,1:end-1)] + s*stateTrackP(j-1,:);
    stateTrackM(j,:) = p*[stateTrackM(j-1,2:end),0] + q*[0,stateTrackM(j-1,1:end-1)] + s*stateTrackM(j-1,:);
    survivalP(j) = sum(stateTrackP(j,:));
    survivalM(j) = sum(stateTrackM(j,:));
    positiveP(j) = sum(stateTrackP(j,abs(thetam)+1:end));
    positiveM(j) = sum(stateTrackM(j,abs(thetam)+1:end));
%     negativeP(j) = sum(stateTrackP(j,1:abs(thetam)-1));
%     negativeM(j) = sum(stateTrackM(j,1:abs(thetam)-1));
%     LP(j) = sum(stateTrackP(j,1:abs(thetam)-1));
%     LM(j) = sum(stateTrackM(j,1:abs(thetam)-1));
end

rP = positiveP./survivalP;
rM = positiveM./survivalM;
RP = positiveP./positiveM;
ag2_social = log(survivalP./survivalM);

if(u1 > 0)
    social = (a - (u1-1))*log(RP);
    for jj = 1:N
        if(EV(tt,jj)<thetap && EV(tt,jj) > thetam)
            EV(tt,jj) = EV(tt,jj) + social(tt);
            if(EV(tt,jj) >= thetap)
                a = a + 1;
            end
        end
    end
    
end

A2(jjj) = A2(jjj) + (a + a1);

end
A1(jjj) = A1(jjj)/(nn*N);
A2(jjj) = A2(jjj)/(nn*N);
end


figure(1)
plot(NNN,A1,'bo','MarkerSize',10 ,'MarkerFaceColor','b')
hold on
plot(NNN,A2,'ro','MarkerSize',10,'MarkerFaceColor','r')
set(gca,'fontsize',20)
xlabel('number of agents'); h = legend('first wave','equilibration end'); set(h,'box','off');
ylabel('fraction of agents correct')
ax = gca;
ax.XTick = [1 100 200];



