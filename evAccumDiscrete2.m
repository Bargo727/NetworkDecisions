clear;

%The Code here computes the relative difference percentage for the decision
%times of Agent 2 when it is isolated versus when it is obtains social
%information from Agent 1.

%Initializing the length of time for which the simulation will run
l = 400;
t = 0:l;
trials = 40;
N = 20000;   %Number of runs over which we average
C = 2; 
%change = .12;
P = [.4, .35, .3];

delta = zeros(trials,1);
ttt = zeros(trials,C+1);

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively

for kk = 0:C
evidence_amounts = [1, -1, 0];
% p = P(kk+1);%exp(1)/7-kk*(change/2);
 s = 0.5;
% q = 1-p-s;


%Threshold values
for thetap = 1:1:trials;
dtime = 0; %Temp variable to compute average relative difference in decision times
thetam = -2;
delta(thetap) = abs(thetap/thetam);
countN = 0;

for jj = 1:N
    
p = P(kk+1);
q = 1-p-s;
   
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
ev11 = ag2_private;

%Initializing the parameter that signifies if agent 1 has reached threshold
%or not 0 corresponds to not done, and 1 corresponds to done
done1 = 0;
done11 = 0;

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

for k = 1:l
    if(done11 == 0)
        if(ev11(k) >= thetap)
            ev11(k) = thetap;
            done11 = 1;
            ag11fire = k;
        elseif(ev11(k) <= thetam)
            ev11(k) = thetam;
            done11 = 1;
            ag11fire = k;
        end
    else
        ev11(k) = ev11(k-1);
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

ag2_social(ag1fire:end) = ev1(ag1fire);

ev2 = zeros(l,1);
done2 = 0; %Definition similar to done1
lastSocialEv = 0;

for j = 1:l
   if(done2 == 0)
       ev2(j) = ag2_private(j) + lastSocialEv;
       if(ev2(j) >= thetap)
           ev2(j) = thetap;
           ag2fire = j;
           done2 = 1;
       elseif(ev2(j) <= thetam)
           ev2(j) = thetam;
           ag2fire = j;
           done2 = 1;
       else
           lastSocialEv = ag2_social(j);
           
       end
   else
       ev2(j) = ev2(j-1);
   end
end
if(ev11(end)>=thetap && ev2(end)>=thetap)
    dtime = dtime + ((ag11fire-ag2fire)/ag11fire);
    countN = countN + 1;
end
end
ttt(thetap,kk+1) = (dtime/countN)*100;
end
end

figure(3)
plot(delta,ttt(:,1),'ro','MarkerSize',7,'MarkerFaceColor','r')
hold on
plot(delta,ttt(:,2),'bo','MarkerSize',7,'MarkerFaceColor','b')
plot(delta,ttt(:,3),'ko','MarkerSize',7,'MarkerFaceColor','k')
xlabel('\Delta\theta')
ylabel('% Difference')
set(gca,'fontsize',20)
h = legend('p = 0.4','p = 0.35','p = 0.3');
set(h,'box','off')
