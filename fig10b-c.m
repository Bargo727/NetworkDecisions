clear;

NNN = [10 50 100 170];
L = length(NNN);

for jjj = 1:L

%Here we will compute some general statistics for the general clique.

%Number of agents in clique
N = NNN(jjj);

%Initializing the length of time for which the simulation will run
l = 1000;
t = 0:l;

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
p = exp(1)/10;
q = 1/10;
s = 1-p-q;

evidence_amounts = [log(p/q), -log(p/q), 0];

%Threshold values
thetap = 20;
thetam = -20;

%We draw the observation of the environment at each time from a true
%multinomial distribution.  Each agent's is independent.
LLR = zeros(l,N);
for jj = 1:N
    obs = mnrnd(1,[p,q,s],l);
    LLR(:,jj) = obs*(evidence_amounts)';
end

EV = cumsum(LLR);

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

a = (N*rP);
d = (N*rM);

%social = (a-d).*(log(rP./rM));

social = (a-d).*(log(RP));

social_per_agent = social./a;

exact = zeros(l,N-1);
for jj = 1:N-1
    exact(:,jj) = factorial(N-1)/(factorial(jj)*factorial(N-1-jj))*(rP.^jj).*((1-rP).^(N-1-jj));
end

compute = sum(exact(:,ceil(N/2):end),2);

st = survivalP.^N;

figure(1)
% subplot(2,2,1)
% plot(t(1:end-1),1-survivalP,'r-o',t(1:end-1),rP,'b-o','LineWidth',3,'MarkerSize',5,'MarkerFaceColor','b')
% set(gca,'fontsize',20)
% xlabel('time')
% ylabel('positive probability')
% 
% subplot(2,2,2)
% plot(t(1:end-1),exact(:,2:2:end),'g-o','LineWidth',3,'MarkerSize',5,'MarkerFaceColor','b')
% set(gca,'fontsize',20)
% xlabel('time')
% ylabel('probability')

% subplot(2,2,4)
% plot(t(1:end-1),social_per_agent,'b-o','LineWidth',3,'MarkerSize',5,'MarkerFaceColor','b')
% set(gca,'fontsize',20)
% xlabel('time')
% ylabel('social evidence per agent')

subplot(1,3,1)
plot(t(1:end-1),social,'-','LineWidth',3,'MarkerSize',5)
set(gca,'fontsize',20)
xlabel('time')
ylabel('social evidence')
% v = legend('N = 10','N=50','N=100','N=200');
% set(v,'box','off')
hold on
% ax = gca;
% ax.XTick = [t(1),t(end-1)];
% ax.XTickLabel = {'1','25'};
% ax.YTick = [floor(social(2)), ceil(social(end))];

subplot(1,3,2)
plot(t(1:end-1),compute,'-','LineWidth',3,'MarkerSize',5)
set(gca,'fontsize',20)
xlabel('time')
ylabel('probability')
 axis([1 20 0 1])
% ax = gca;
% ax.XTick = [1, 5];
% ax.XTickLabel = {'1','5'};
% ax.YTick = [0 0.5 1];
h = legend('N = 10','N = 50','N = 100','N = 200');
set(h,'box','off')
hold on

subplot(1,3,3)
plot(t(1:end-1),st,'-','LineWidth',3,'MarkerSize',5)
set(gca,'fontsize',20)
xlabel('time')
ylabel('probability')
axis([1 100 0 1])
hold on

end


