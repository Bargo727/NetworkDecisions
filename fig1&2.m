clear;

%Initializing the length of time for which the simulation will run
l = 20;
t = 0:l;

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
evidence_amounts = [1, -1, 0];
p = exp(1)/7;
q = 1/7;
s = 1-p-q;

%Threshold values
thetap = 30;
thetam = -2;

%We draw the observation of the environment at each time from a true
%multinomial distribution.  Each agent's is independent.
obs1 = mnrnd(1,[p,q,s],l);
obs2 = mnrnd(1,[p,q,s],l);

%The information obtained from the observations described above.
llr1 = obs1*(evidence_amounts)';
llr2 = obs2*(evidence_amounts)';

%The sum to obtain the log-likelihood ratio at a time t
ev1 = cumsum(llr1);
ag2_private = cumsum(llr2);

%Initializing the parameter that signifies if agent 1 has reached threshold
%or not 0 corresponds to not done, and 1 corresponds to done
done1 = 0; 

for k = 1:l
    if(done1 == 0)
        if(ev1(k) >= thetap)
            ev1(k) = thetap;
            done1 = 1;
            ag1fire = k;
        elseif(ev1(k) <= thetam)
            ev1(k) = thetam;
            done1 = 1;
            ag1fire = k;
        end
    else
        ev1(k) = ev1(k-1);
    end
        
end

%Now we compute Agent 2's Social Evidence 

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

for j = 2:l
    stateTrackP(j,:) = q*[stateTrackP(j-1,2:end),0] + p*[0,stateTrackP(j-1,1:end-1)] + s*stateTrackP(j-1,:);
    stateTrackM(j,:) = p*[stateTrackM(j-1,2:end),0] + q*[0,stateTrackM(j-1,1:end-1)] + s*stateTrackM(j-1,:);
    survivalP(j) = sum(stateTrackP(j,:));
    survivalM(j) = sum(stateTrackM(j,:));
end

ag2_social = log(survivalP./survivalM);
ag2_social(ag1fire+1:end) = ev1(ag1fire);

ev2 = zeros(l,1);
done2 = 0; %Definition similar to done1
lastSocialEv = 0;

for j = 1:l
   if(done2 == 0)
       ev2(j) = ev2(j) + ag2_private(j) + lastSocialEv;
       if(ev2(j) >= thetap)
           ev2(j) = thetap;
           ag2fire = j;
           done2 = 1;
       elseif(ev2(j) <= thetam)
           ev2(j) = thetam;
           ag2fire = j;
           done2 = 1;
       else
          % ev2(j) = ev2(j) + ag2_social(j) - lastSocialEv;
           lastSocialEv = ag2_social(j);
%            if(ev2(j) >= thetap)
%                 ev2(j) = thetap;
%                 ag2fire = j;
%                 done2 = 1;
%            elseif(ev2(j) <= thetam)
%                 ev2(j) = thetam;
%                 ag2fire = j;
%                 done2 = 1;
%            end
%            
       end
   else
       ev2(j) = ev2(j-1);
   end
end

ev1 = [0; ev1];
ag2_private = ag2_private(1:ag2fire+1);
ag2_social = [0, ag2_social];
ev2 = [0; ev2];

ll = length(ag2_private);
tt = 0:ll-1;


figure(1)
subplot(2,1,1)
plot(t,ev1,'k-o',t,ev2,'r-o','LineWidth',5,'MarkerSize',10)
ax = gca;
set(ax,'fontsize',20)
h = legend('Agent 1','Agent 2');
%xlabel('Time')
ylabel('Evidence','fontsize',25)
set(h,'box','off')
ax.XTick = [];
ax.YTick = [thetam thetap];
ax.YTickLabel = {'\theta','-\theta'};
axis([0 l thetam-.2 thetap+.2])

subplot(2,1,2)
plot(t,ag2_social,'-o','Color',[1 .44 0],'LineWidth',5,'MarkerSize',5)
hold on
plot(tt,ag2_private,'-o','Color',[.75 .1 .57],'LineWidth',5,'MarkerSize',5)
ax = gca;
set(ax,'fontsize',20)
xlabel('Time', 'fontsize',25)
h = legend('Social','Private');
set(h,'box','off')
ax.XTick = [0 l];
ax.YTick = [thetam thetap];
axis([0 l thetam-.2 thetap+.2])
