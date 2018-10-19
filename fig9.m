clear;

NN = 10000;
thrs = [1.5/7 2/7 3/7 4/7 5/7 6/7];

PC = zeros(length(thrs),1);
PCU = zeros(length(thrs),1);
DT = zeros(length(thrs),1);
DTU = zeros(length(thrs),1);

for jjj = 1:length(thrs)
    count = 0;
    countu = 0;
    timecount = 0;
    timecountu = 0;

for jj = 1:NN

%Initializing the length of time for which the simulation will run
l = 100;
t1 = 0:l;

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
p = thrs(jjj);
q = 1/7;
s = 1-p-q;

evidence_amounts = [log(p/q), -log(p/q), 0];

%Threshold values
thr = 2;
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

ag11_private = cumsum(llr1);
ag22_private = cumsum(llr2);
ag33_private = cumsum(llr3);

done1 = 0;
done2 = 0;
done3 = 0;

done11 = 0;
done22 = 0;
done33 = 0;

done = [done1;done2;done3];
ag1fire = 0;
ag2fire = 0;
ag3fire = 0;

ag11fire = 0;
ag22fire = 0;
ag33fire = 0;

for k = 1:l
    if(done11 == 0)
        if(ag11_private(k) >= thetap)
            ag11_private(k) = thetap;
            done11 = 1;
            ag11fire = k;
        elseif(ag11_private(k) <= thetam)
            ag11_private(k) = thetam;
            done11 = 1;
            ag11fire = k;
        end
    else
        ag11_private(k) = ag11_private(k-1);
    end
    
    if(done22 == 0)
        if(ag22_private(k) >= thetap)
            ag22_private(k) = thetap;
            done22 = 1;
            ag22fire = k;
        elseif(ag22_private(k) <= thetam)
            ag22_private(k) = thetam;
            done22 = 1;
            ag22fire = k;
        end
    else
        ag22_private(k) = ag22_private(k-1);
    end
    
    if(done33 == 0)
        if(ag33_private(k) >= thetap)
            ag33_private(k) = thetap;
            done33 = 1;
            ag33fire = k;
        elseif(ag33_private(k) <= thetam)
            ag33_private(k) = thetam;
            done33 = 1;
            ag33fire = k;
        end
    else
        ag33_private(k) = ag33_private(k-1);
    end
end

if(ag11_private(end) >0 && ag22_private(end) >0 && ag33_private(end) >0)
    countu = countu + 1;
end

timecountu = timecountu + max([ag11fire,ag22fire,ag33fire]);

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
decider = find(valid_inds>0);
social = zeros(l,1);



%Now we determine which of three cases we obtained.  The three cases are:
%(1) both remaining agents leaned toward the first agent's decision and
%therefore immediately make the same decision as the first agent when they
%receive the information from the original decider, (2) one of the
%remaining agents agrees with the first decider and the other does not, and
%(3) the remaining agents do not agree with the first. 

if(decider == 1)
    social(ag1fire+1:end) = sign(ag1_private(end))*thr;
    ev1 = ag1_private;
    ev2 = ag2_private + social;
    ev3 = ag3_private + social;
    done2 = 0;
    done3 = 0;
    done = [done2;done3];
for k = ag1fire+1:l
    if(done == 0)
        if(ev2(k) >= thetap)
            ev2(k) = thetap;
            done2 = 1;
            ag2fire = k;
        elseif(ev2(k) <= thetam)
            ev2(k) = thetam;
            done2 = 1;
            ag2fire = k;
        end
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
        if(done2 == 1)
            ev2(k) = ev2(k-1);
        end
        if(done3 == 1)
            ev3(k) = ev3(k-1);
        end
   
    end
    done = [done2;done3];
end
       agents = [ev2, ev3];
       agfire = [ag2fire;ag3fire];
       valid_inds = (agfire>0)';
       first_fire = agfire(valid_inds);
       boss = agents(:,valid_inds);
       decider1 = find(valid_inds>0,1);
       social = zeros(l,1);
       if(decider1 == 1)
           social(ag2fire+1:end) = sign(ev2(end))*thr;
           ev3 = ev3 + social;
           done3 = 0;
                    for k = ag2fire+1:l
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
       elseif(decider1 == 2)
           social(ag3fire+1:end) = sign(ev3(end))*thr;
           ev2 = ev2 + social;
           done2 = 0;
           
                    for k = ag3fire+1:l
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
       end
elseif(decider == 2)
    social(ag2fire+1:end) = sign(ag1_private(end))*thr;
    ev1 = ag1_private + social;
    ev2 = ag2_private;
    ev3 = ag3_private + social;
    done1 = 0;
    done3 = 0;
    done = [done1;done3];
for k = ag2fire+1:l
    if(done == 0)
        if(ev1(k) >= thetap)
            ev1(k) = thetap;
            done1 = 1;
            ag2fire = k;
        elseif(ev1(k) <= thetam)
            ev2(k) = thetam;
            done1 = 1;
            ag2fire = k;
        end
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
        if(done1 == 1)
            ev1(k) = ev1(k-1);
        end
        if(done3 == 1)
            ev3(k) = ev3(k-1);
        end
   
    end
    done = [done1;done3];
end
       agents = [ev1, ev3];
       agfire = [ag1fire;ag3fire];
       valid_inds = (agfire>0)';
       first_fire = agfire(valid_inds);
       boss = agents(:,valid_inds);
       decider1 = find(valid_inds>0,1);
       social = zeros(l,1);
       if(decider1 == 1)
           social(ag1fire+1:end) = sign(ev1(end))*thr;
           ev3 = ev3 + social;
           done3 = 0;
                    for k = ag1fire+1:l
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
       elseif(decider1 == 2)
           social(ag3fire+1:end) = sign(ev3(end))*thr;
           ev1 = ev1 + social;
           done1 = 0;
           
                    for k = ag3fire+1:l
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
       end
elseif(decider == 3)
    social(ag3fire+1:end) = sign(ag1_private(end))*thr;
    ev2 = ag2_private + social;
    ev1 = ag1_private + social;
    ev3 = ag3_private;
    done2 = 0;
    done1 = 0;
    done = [done1;done2];
for k = ag3fire+1:l
    if(done == 0)
        if(ev1(k) >= thetap)
            ev1(k) = thetap;
            done1 = 1;
            ag1fire = k;
        elseif(ev1(k) <= thetam)
            ev2(k) = thetam;
            done1 = 1;
            ag1fire = k;
        end
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
        if(done1 == 1)
            ev1(k) = ev1(k-1);
        end
        if(done2 == 1)
            ev2(k) = ev2(k-1);
        end
   
    end
    done = [done1;done2];
end
       agents = [ev2, ev1];
       agfire = [ag2fire;ag1fire];
       valid_inds = (agfire>0)';
       first_fire = agfire(valid_inds);
       boss = agents(:,valid_inds);
       decider1 = find(valid_inds>0,1);
       social = zeros(l,1);
       if(decider1 == 1)
           social(ag2fire+1:end) = sign(ev2(end))*thr;
           ev1 = ev1 + social;
           done1 = 0;
                    for k = ag1fire+1:l
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
       elseif(decider1 == 2)
           social(ag1fire+1:end) = sign(ev1(end))*thr;
           ev2 = ev2 + social;
           done2 = 0;
           
                    for k = ag1fire+1:l
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
       end
end

if(ev1(end) >0 && ev2(end) >0 && ev3(end) >0)
    count = count + 1;
end

timecount = timecount + max([ag1fire,ag2fire,ag3fire]);


end
PC(jjj) = count/NN;
PCU(jjj) = countu/NN;
DT(jjj) = timecount/NN;
DTU(jjj) = timecountu/NN;
end
ev1 = [0;ev1];
ev2 = [0;ev2];
ev3 = [0;ev3];


% plot(t1,ev1','r',t1,ev2','b',t1,ev3','k','LineWidth',5)

figure(1)
subplot(1,2,1)
plot(thrs/q,PC,'bo','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(thrs/q,PCU,'ko','MarkerSize',10,'MarkerFaceColor','k')
set(gca,'fontsize',20)
ylabel('percent of trials correct')
h = legend('Coupled','Uncoupled');
set(h,'box','off')

subplot(1,2,2)
plot(thrs/q,DT,'bo','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(thrs/q,DTU,'ko','MarkerSize',10,'MarkerFaceColor','k')
set(gca,'fontsize',20)
ylabel('decision time')


