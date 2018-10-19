clear;

%This script will produce a plot of the percentage of trials for which the three agents in the unidirectional line choose the correct
%choice.

ps = [1.5/7, 2/7, 2.5/7, 3/7, 3.5/7, 4/7,4.5/7, 5/7,5.5/7,6/7];
NN = length(ps);
PC = zeros(NN,1);
PCU = zeros(NN,1);
DT = zeros(NN,1);
DTU = zeros(NN,1);

count = 10000;

%Initializing the length of time for which the simulation will run
l = 500;
t = 0:l;

%Initializing parameters

for jj = 1:NN

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
p = ps(jj);
q = 1/7;
s = 1-p-q;
evidence_amounts = [log(p/q), -log(p/q), 0];

%Threshold values
thr = 2;
thetap = thr;
thetam = -thr;

%Probability of getting correct answer
alpha = exp(thr)/(1+exp(thr));
dt = 0;
dtu = 0;
cc = 0;
ccu = 0;
for jjj = 1:count

%We draw the observation of the environment at each time from a true
%multinomial distribution.  Each agent's is independent.
obs1 = mnrnd(1,[p,q,s],l);
obs2 = mnrnd(1,[p,q,s],l);
obs3 = mnrnd(1,[p,q,s],l);

%The information obtained from the observations described above.
llr1 = obs1*(evidence_amounts)';
llr2 = obs2*(evidence_amounts)';
llr3 = obs3*(evidence_amounts)';

%The sum to obtain the log-likelihood ratio at a time t
ev1 = cumsum(llr1);
ag2_private = cumsum(llr2);
ag3_private = cumsum(llr3);

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

done22 = 0;
for k = 1:l
    if(done22 == 0)
        if(ag2_private(k) >= thetap)
            ag2_private(k) = thetap;
            done22 = 1;
            ag22fire = k;
        elseif(ag2_private(k) <= thetam)
            ag2_private(k) = thetam;
            done22 = 1;
            ag22fire = k;
        end
    else
        ag2_private(k) = ag2_private(k-1);
    end
        
end

done33 = 0;
for k = 1:l
    if(done33 == 0)
        if(ag3_private(k) >= thetap)
            ag3_private(k) = thetap;
            done33 = 1;
            ag33fire = k;
        elseif(ag3_private(k) <= thetam)
            ag3_private(k) = thetam;
            done33 = 1;
            ag33fire = k;
        end
    else
        ag3_private(k) = ag3_private(k-1);
    end
        
end

%Now we look at Agent 2's evidence and add in the social evidence if
%appropriate.

done2 = 0;
ag2_social = zeros(l,1);
ag2_social(ag1fire+1:end) = sign(ev1(end))*thr;

ev2 = ag2_private + ag2_social;
for k = 1:l
    if(done2 == 0)
        if(ev2(k) >= thetap)
            ev2(k) = thetap;
            done2 = 1;
            ag2fire = k;
        elseif(ev2(k) <= thetam)
            ev2(k) = thetam;
            done2 = 1;
            ag2fire = k;
        end
    else
        ev2(k) = ev2(k-1);
    end
end

%Now we compute Agent 3's Social Evidence 

%Indexing Guide for the small threshold example
    %-2 -1 0 1 2 3
    %     0 1 2 3
 
 ag3_social = zeros(l,1);
 
 total_states = thetap+abs(thetam) + 1 - 2; %Include 0 in the states but not the boundary values
 stateTrackM = zeros(ag2fire-1,total_states);
 stateTrackM(1,:) = zeros(1,total_states);
 stateTrackM(1,abs(thetam)) = 1;

stateTrackP = zeros(ag2fire-1,total_states);
stateTrackP(1,:) = zeros(1,total_states);
stateTrackP(1,abs(thetam)) = 1;
 
survivalP = zeros(1,ag2fire-1);
survivalM = zeros(1,ag2fire-1);
survivalP(1) = 1;
survivalM(1) = 1;

for j = 2:ag2fire-1
    stateTrackP(j,:) = q*[stateTrackP(j-1,2:end),0] + p*[0,stateTrackP(j-1,1:end-1)] + s*stateTrackP(j-1,:);
    stateTrackM(j,:) = p*[stateTrackM(j-1,2:end),0] + q*[0,stateTrackM(j-1,1:end-1)] + s*stateTrackM(j-1,:);
end

survivalP(2:end) = sum(stateTrackP(2:end,abs(thetam):end),2);
 survivalM(2:end) = sum(stateTrackM(2:end,abs(thetam):end),2);
 idk = survivalP./survivalM;


ag3social = zeros(l,1);
ag3social(ag2fire+1:end) = sign(ev2(end))*thr + (ev2(end)>0)*log(idk(end)) + (ev2(end)<0)*log(1/idk(end));

ev3 = ag3_private + ag3social;
done3 = 0; 

for k = 1:l
    if(done3 == 0)
        if(ev3(k) >= thetap)
            ev3(k) = thetap;
            done3 = 1;
            ag3fire = k;
        elseif(ev3(k) <= thetam)
            ev3(k) = thetam;
            done3 = 1;
            ag3fire = k;
        end
    else
        ev3(k) = ev3(k-1);
    end
end
dt = dt + max([ag1fire,ag2fire,ag3fire]);
dtu = dtu + max([ag1fire,ag22fire,ag33fire]);
% if((ev1(end) > 0 && ev2(end)) > 0 || (ev1(end) > 0 && ev3(end)>0)||(ev3(end) > 0 && ev2(end)>0))
%     cc = cc + 1;
% end
% if(ev1(end) > 0 && ag2_private(end) > 0 || ev1(end) > 0 && ag3_private(end)>0||ag3_private(end) > 0 && ag2_private(end)>0)
%     ccu = ccu + 1;
% end
if(ev1(end) > 0 && ev2(end)>0 && ev3(end)>0)
    cc = cc+1;
end
if(ev1(end)>0 && ag2_private(end) >0 && ag3_private(end)>0)
    ccu = ccu + 1;
end
PC(jj) = cc/count;
PCU(jj) = ccu/count;
DT(jj) = dt/count;
DTU(jj) = dtu/count;
end
end


figure(1)
subplot(1,2,1)
plot(ps/q,PC,'bo','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(ps/q,PCU,'ko','MarkerSize',10,'MarkerFaceColor','k')
set(gca,'fontsize',20)
h = legend('Coupled','Uncoupled');
set(h,'box','off')
xlabel('p/q')
ylabel('Percentage of Correct Trials')

subplot(1,2,2)
plot(ps/q,DT,'bo','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(ps/q,DTU,'ko','MarkerSize',10,'MarkerFaceColor','k')
set(gca,'fontsize',20)
ylabel('Decision Time')


