clear;

%Initializing the length of time for which the simulation will run
l = 60;
t1 = 0:l;

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
p = exp(1)/6;
q = 1/6;
s = 1-p-q;

evidence_amounts = [log(p/q), -log(p/q), 0];

%Threshold values
thr = 3;
thetap = thr;
thetam = -thr;

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
ag1_private = cumsum(llr1);
ag2_private = cumsum(llr2);
ag3_private = cumsum(llr3);

done1 = 0;
done2 = 0;
done3 = 0;

done = [done1;done2;done3];
ag1fire = 0;
ag2fire = 0;
ag3fire = 0;

for k = 1:l
    if(done == 0)
        if(ag1_private(k) >= thetap)
            ag1_private(k) = thetap;
            done1 = 1;
            ag1fire = k;
        elseif(ag1_private(k) <= thetam)
            ag1_private(k) = thetam;
            done1 = 1;
            ag1fire = k;
        end
        if(ag2_private(k) >= thetap)
            ag2_private(k) = thetap;
            done2 = 1;
            ag2fire = k;
        elseif(ag2_private(k) <= thetam)
            ag2_private(k) = thetam;
            done2 = 1;
            ag2fire = k;
        end
        if(ag3_private(k) >= thetap)
            ag3_private(k) = thetap;
            done3 = 1;
            ag3fire = k;
        elseif(ag3_private(k) <= thetam)
            ag3_private(k) = thetam;
            done3 = 1;
            ag3fire = k;
        end
    else
        if(done1 == 1)
            ag1_private(k) = ag1_private(k-1);
        end
        if(done2 == 1)
            ag2_private(k) = ag2_private(k-1);
        end
        if(done3 == 1)
            ag3_private(k) = ag3_private(k-1);
        end
   
    end
    done = [done1;done2;done3];
end

agents = [ag1_private, ag2_private, ag3_private];
agfire = [ag1fire;ag2fire;ag3fire];
valid_inds = (agfire>0)';
first_fire = agfire(valid_inds);
boss = agents(:,valid_inds);
decider = find(valid_inds>0,1);

%Now we determine which of three cases we obtained.  The three cases are:
%(1) both remaining agents leaned toward the first agent's decision and
%therefore immediately make the same decision as the first agent when they
%receive the information from the original decider, (2) one of the
%remaining agents agrees with the first decider and the other does not, and
%(3) the remaining agents do not agree with the first. 

state_d = sign(boss(first_fire));
if(done1  == 1)
    state_1 = sign(agents(first_fire+1,2));
    state_2 = sign(agents(first_fire+1,3));
elseif(done2 == 1)
    state_1 = sign(agents(first_fire+1,1));
    state_2 = sign(agents(first_fire+1,3));
else
    state_1 = sign(agents(first_fire+1,1));
    state_2 = sign(agents(first_fire+1,2));
end

if(state_d == state_1 && state_d == state_2)
    situation = 1;
elseif((state_d == state_1 && state_d ~= state_2) || (state_d == state_2 && state_d ~= state_1))
    situation = 2;
else
    situation = 3;
end

%Now that the independent accumulation of evidence phase is over, we will
%begin accounting for social information exchanged between the agents.
%First, from the first decision, there is a kick of size theta_+ or theta_-
%to the remaining agents.

decision_kick = sign(boss(first_fire))*thr;

for k = first_fire+1:l
    if(done1 == 1)
        ag2_private(k) = ag2_private(k) + decision_kick;
        ag3_private(k) = ag3_private(k) + decision_kick;
        if(ag2_private(k) >= thetap)
            ag2_private(k) = thetap;
        elseif(ag2_private(k) <= thetam)
            ag2_private(k) = thetam;
        end
        if(ag3_private(k) >= thetap)
            ag3_private(k) = thetap;
        elseif(ag3_private(k) <= thetam)
            ag3_private(k) = thetam;
        end
    elseif(done2 == 1)
        ag1_private(k) = ag1_private(k) + decision_kick;
        ag3_private(k) = ag3_private(k) + decision_kick;
        if(ag1_private(k) >= thetap)
            ag1_private(k) = thetap;
        elseif(ag1_private(k) <= thetam)
            ag1_private(k) = thetam;
        end
        if(ag3_private(k) >= thetap)
            ag3_private(k) = thetap;
        elseif(ag3_private(k) <= thetam)
            ag3_private(k) = thetam;
        end
    else
        ag1_private(k) = ag1_private(k) + decision_kick;
        ag2_private(k) = ag2_private(k) + decision_kick;
        if(ag1_private(k) >= thetap)
            ag1_private(k) = thetap;
        elseif(ag1_private(k) <= thetam)
            ag1_private(k) = thetam;
        end
        if(ag2_private(k) >= thetap)
            ag2_private(k) = thetap;
        elseif(ag2_private(k) <= thetam)
            ag2_private(k) = thetam;
        end
    end
    
    
end

z = zeros(length(t1),1);


figure(1)
plot(t1,[0;ag1_private]','k-o',t1,[0;ag2_private]','r-o',t1,[0;ag3_private]','b-o',t1,z,'k--','LineWidth',4,'MarkerSize',7)
set(gca,'fontsize',20)
ax = gca;
hh = legend('Agent 1','Agent 2','Agent 3');
set(hh,'box','off')
ax.XTick = [0 first_fire];
ax.YTick = [thetam 0 thetap];
%axis([0 1.3*first_fire thetam thetap])
