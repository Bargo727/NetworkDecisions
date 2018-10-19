clear;

%The Code here computes the first passage time density for coupled versus
%uncoupled agent two decision times.

%Initializing the length of time for which the simulation will run
l = 500;
t = 0:l;
N = 120000;   %Number of runs over which we average


ttt = zeros(N,2);

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively

evidence_amounts = [1, -1, 0];
p = exp(1)/7;
q = 1/7;
s = 1-p-q;

%Threshold values
thetap = 10;
thetam = -2;

for jj = 1:N
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

ttt(jj,1) = ag11fire;

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

ttt(jj,2) = ag2fire;


end

[H1,C1] = hist(ttt(:,1),100);
[H2,C2] = hist(ttt(:,2),100);

figure(3)
plot(C1,H1,'k',C2,H2,'b','LineWidth',4)
set(gca,'fontsize',20)
ax = gca;
ax.YTick = [];
h = legend('Uncoupled','Coupled');
set(h,'box','off')
