clear;

%The end plot here shows the percentage of trials for which both agents
%chose the correct choice for the cases where they are coupled and the cases
%where they are uncoupled.

%Initializing the length of time for which the simulation will run
l = 400;
t = 0:l;
trials = 20;
N = 2000;
C = 2;
%change = .12;
P = [.4 .35 .3];

%vector to store d\theta 
delta = zeros(trials,1);

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
prc = zeros(trials,C+1);
pru = zeros(trials,C+1);
for kk = 0:C
evidence_amounts = [1, -1, 0];
s = 0.5;





%Threshold values
for thetap = 1:trials;
thetam = -2;
delta(thetap) = abs(thetap/thetam);
%delta = [delta;thetap-thetam];
countc = 0;
countu = 0;

for jj = 1:N
    
    if(mod(jj,2)==0)
        p = P(kk+1);
        q = 1-p-s;
    else
        q = P(kk+1);
        p = 1-q-s;
    end

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

%if(mod(jj,2)==0)
    ag2_social = log(survivalP./survivalM);
% else
%     ag2_social = log(survivalM./survivalP);
% end
ag2_social(ag1fire:end) = ev1(ag1fire);

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
           lastSocialEv = ag2_social(j);            
       end
   else
       ev2(j) = ev2(j-1);
   end
end

if(mod(jj,2)==0)
    if(ev1(ag1fire) >= thetap && ev2(ag2fire) >= thetap)
        countc = countc + 1;
    end

    if(ev1(ag1fire) >= thetap && ev11(ag11fire) >= thetap)
        countu = countu + 1;
    end
else
    if(ev1(ag1fire) <= thetam && ev2(ag2fire) <= thetam)
        countc = countc + 1;
    end

    if(ev1(ag1fire) <= thetam && ev11(ag11fire) <= thetam)
        countu = countu + 1;
    end
end

end
% prc = [prc;100*countc/N];
% pru = [pru;100*countu/N];
prc(thetap,kk+1) = 100*(countc/N);
pru(thetap,kk+1) = 100*(countu/N);
end
end

Ph = P./(P + (.5-P));

Prr = ((1./(2*Ph)).*(1-sqrt(1-4*Ph.*(1-Ph)))).^2;

LL = length(delta);

prr1 = .5*(2-Prr(1))*ones(LL,1);
prr2 = .5*(2-Prr(2))*ones(LL,1);
prr3 = .5*(2-Prr(3))*ones(LL,1);



figure(3)
subplot(3,1,1)
plot(delta,prc(:,1),'co','MarkerSize',8,'MarkerFaceColor','c')
hold on
plot(delta,pru(:,1),'ko','MarkerSize',8,'MarkerFaceColor','k')
plot(delta,100*prr1,'r--','LineWidth',3)
set(gca,'fontsize',20)
h = legend('Coupled','Uncoupled');
set(h,'box','off')
axis([0 max(delta) 70 100] )
ax = gca;
ax.XTick = [];
ax.YTick = [70 100];
title('p = 0.4')

subplot(3,1,2)
plot(delta,prc(:,2),'co','MarkerSize',8,'MarkerFaceColor','c')
hold on
plot(delta,pru(:,2),'ko','MarkerSize',8,'MarkerFaceColor','k')
plot(delta,100*prr2,'r--','LineWidth',3)
set(gca,'fontsize',20)
axis([0 max(delta) 50 100])
ax = gca;
ax.XTick = [];
ax.YTick = [50 100];
title('p = 0.35')
ylabel('% Both Agents Correct')

subplot(3,1,3)
plot(delta,prc(:,3),'co','MarkerSize',8,'MarkerFaceColor','c')
hold on
plot(delta,pru(:,3),'ko','MarkerSize',8,'MarkerFaceColor','k')
plot(delta,100*prr3,'r--','LineWidth',3)
set(gca,'fontsize',20)
axis([0 max(delta) 40 100])
ax = gca;
ax.XTick = [1,max(delta)];
ax.YTick = [50 100];
title('p = 0.3')
xlabel('|\theta_{+}/\theta_{-}|')
