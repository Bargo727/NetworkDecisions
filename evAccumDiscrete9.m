clear;

%Here we simulate the recurrent two-agent network.

%Initializing the length of time for which the simulation will run
l = 20;
t = 0:l;

%Initializing parameters

%We allow jumps of unit "information" with probabilities p, q, and s,
%respectively
p = exp(1)/7;
q = 1/7;
s = 1-p-q;
evidence_amounts = [log(p/q), log(q/p), 0];

%Threshold values
thetap = 2;
thetam = -1;

%We draw the observation of the environment at each time from a true
%multinomial distribution.  Each agent's is independent.
obs1 = mnrnd(1,[p,q,s],l);
obs2 = mnrnd(1,[p,q,s],l);

%The information obtained from the observations described above.
llr1 = obs1*(evidence_amounts)';
llr2 = obs2*(evidence_amounts)';

%The sum to obtain the log-likelihood ratio at a time t
ev1 = cumsum(llr1);
ev2 = cumsum(llr2);

%Initializing the parameter that signifies if agent 1 has reached threshold
%or not 0 corresponds to not done, and 1 corresponds to done

done1 = 0;
done2 = 0;
social = 0;

equal1 = [];
equal2 = [];

tol = 10E-4;

for k = 1:l
    k
    if(ev1(k) < thetap && ev2(k) < thetap && ev1(k) > thetam && ev2(k) > thetam)
        
        
        %Equilibration after first observation
        if(k == 1)
            
            %Indexing Guide for the small threshold example
            %-2 -1 0 1 2 3
             %   0 1 2 3
             
            total_states = thetap+abs(thetam) + 1 - 2; %Include 0 in the states but not the boundary values
            
            stateTrackM = zeros(l,total_states);
            stateTrackM(1,abs(thetam)) = 1;

            stateTrackP = zeros(l,total_states);
            stateTrackP(1,abs(thetam)) = 1;
            
            stateTrackP(2,:) = q*[stateTrackP(1,2:end),0] + p*[0,stateTrackP(1,1:end-1)] + s*stateTrackP(1,:);
            stateTrackM(2,:) = p*[stateTrackM(1,2:end),0] + q*[0,stateTrackM(1,1:end-1)] + s*stateTrackM(1,:);
            
            survivalP = sum(stateTrackP(2,:));
            survivalM = sum(stateTrackM(2,:));
            social = log(survivalP./survivalM);
            
            equal1 = [equal1;social];
            
            check = 0;
            
            while(check == 0)
                nev1 = ev1(k) + social;
                nev2 = ev2(k) + social;
                if(nev1 >= thetap)
                    check = 1;
                    ev1(k:end) = thetap;
                    ev2(k:end) = ev2(k:end) + thetap;
                    done2 = 0;
                    for kk = 1:l
                        if(done2 == 0)
                            if(ev2(kk) >= thetap)
                                ev2(kk) = thetap;
                                done2 = 1;
                                %ag1fire = kk;
                            elseif(ev2(kk) <= thetam)
                                ev2(kk) = thetam;
                                done2 = 1;
                                %ag1fire = kk;
                            end
                        else
                            ev2(kk) = ev2(kk-1);
                        end

                    end
                elseif(nev2 >= thetap)
                    check = 1;
                    ev2(k:end) = thetap;
                    ev1(k:end) = ev1(k:end) + thetap;
                    done1 = 0;
                    for kk = 1:l
                        if(done1 == 0)
                            if(ev1(kk) >= thetap)
                                ev1(kk) = thetap;
                                done1 = 1;
                                %ag1fire = kk;
                            elseif(ev1(kk) <= thetam)
                                ev1(kk) = thetam;
                                done1 = 1;
                                %ag1fire = kk;
                            end
                        else
                            ev1(kk) = ev1(kk-1);
                        end

                    end
                else
                    survivalP = sum(stateTrackP(2,1:end-floor(social)));
                    survivalM = sum(stateTrackM(2,1:end-floor(social)));
                    new_social = log(survivalP./survivalM);
                    if(new_social == social)
                        check = 1;
                    else
                        equal1 = [equal1;new_social];
                    end
                    social = new_social;
                end
            end
                if(done1 ==0 && done2 ==0)
                    ev1(k) = ev1(k) + social;
                    ev2(k) = ev2(k) + social;
                end
                
        %Equilibration after an arbitrary observation
        elseif(k > 1)
            stateTrackP(k+1,:) = q*[stateTrackP(k,2:end),0] + p*[0,stateTrackP(k,1:end-1)] + s*stateTrackP(k,:);
            stateTrackM(k+1,:) = p*[stateTrackM(k,2:end),0] + q*[0,stateTrackM(k,1:end-1)] + s*stateTrackM(k,:);
            
            survivalP = sum(stateTrackP(k+1,1:end-floor(social)));
            survivalM = sum(stateTrackM(k+1,1:end-floor(social)));
            social = log(survivalP./survivalM);
            
            if(k == 2)
                equal2 = [equal2;social];
            end
            
            check = 0;
            
            while(check == 0)
                nev1 = ev1(k) + social;
                nev2 = ev2(k) + social;
                if(nev1 >= thetap)
                    check = 1;
                    ev1(k:end) = thetap;
                    ev2(k:end) = ev2(k:end) + thetap;
                    done2 = 0;
                    for kk = 1:l
                        if(done2 == 0)
                            if(ev2(kk) >= thetap)
                                ev2(kk) = thetap;
                                done2 = 1;
                                %ag1fire = kk;
                            elseif(ev2(kk) <= thetam)
                                ev2(kk) = thetam;
                                done2 = 1;
                                %ag1fire = kk;
                            end
                        else
                            ev2(kk) = ev2(kk-1);
                        end

                    end
                elseif(nev2 >= thetap)
                    check = 1;
                    ev2(k:end) = thetap;
                    ev1(k:end) = ev1(k:end) + thetap;
                    done1 = 0;
                    for kk = 1:l
                        if(done1 == 0)
                            if(ev1(kk) >= thetap)
                                ev1(kk) = thetap;
                                done1 = 1;
                                %ag1fire = kk;
                            elseif(ev1(kk) <= thetam)
                                ev1(kk) = thetam;
                                done1 = 1;
                                %ag1fire = kk;
                            end
                        else
                            ev1(kk) = ev1(kk-1);
                        end

                    end
                else
                    survivalP = sum(stateTrackP(k+1,1:end-floor(social)));
                    survivalM = sum(stateTrackM(k+1,1:end-floor(social)));
                    new_social = log(survivalP./survivalM);
                    
                    
                    if(abs(new_social-social)<tol)
                        check = 1;
                    else
                        if(k ==2)
                            equal2 = [equal2;new_social];
                        end
                    end
                    social = new_social;
                end
            end
                if(done1 ==0 && done2 ==0)
                    ev1(k) = ev1(k) + social;
                    ev2(k) = ev2(k) + social;
                end
        end
        
    elseif(ev1(k) >= thetap)
        ev1(k:end) = thetap;
        ev2(k:end) = ev2(k:end) + thetap;
        done2 = 0;
        for kk = k:l
            if(done2 == 0)
                if(ev2(kk) >= thetap)
                    ev2(kk) = thetap;
                    done2 = 1;
                elseif(ev2(kk) <= thetam)
                    ev2(kk) = thetam;
                    done2 = 1;
                end
            else
                ev2(kk) = ev2(kk-1);
            end
        end
    elseif(ev1(k) <= thetam)
        ev1(k:end) = thetam;
        ev2(k:end) = ev2(k:end) + thetam;
        done2 = 0;
        for kk = k:l
            if(done2 == 0)
                if(ev2(kk) >= thetap)
                    ev2(kk) = thetap;
                    done2 = 1;
                elseif(ev2(kk) <= thetam)
                    ev2(kk) = thetam;
                    done2 = 1;
                end
            else
                ev2(kk) = ev2(kk-1);
            end
        end
    elseif(ev2(k) >= thetap)
        ev2(k:end) = thetap;
        ev1(k:end) = ev1(k:end) + thetap;
        done2 = 0;
        for kk = k:l
            if(done2 == 0)
                if(ev1(kk) >= thetap)
                    ev1(kk) = thetap;
                    done2 = 1;
                elseif(ev2(kk) <= thetam)
                    ev1(kk) = thetam;
                    done2 = 1;
                end
            else
                ev1(kk) = ev1(kk-1);
            end
        end
    elseif(ev2(k) <= thetam)
        ev2(k:end) = thetam;
        ev1(k:end) = ev1(k:end) + thetam;
        done2 = 0;
        for kk = k:l
            if(done2 == 0)
                if(ev1(kk) >= thetap)
                    ev1(kk) = thetap;
                    done2 = 1;
                elseif(ev1(kk) <= thetam)
                    ev1(kk) = thetam;
                    done2 = 1;
                end
            else
                ev1(kk) = ev1(kk-1);
            end
        end
    end
end

ev1 = [0;ev1];
ev2 = [0;ev2];
v11 = equal1(end)*ones(5,1);
v12 = equal2(end)*ones(5,1);
equal1 = [0;equal1;v11];
equal2 = [0;equal2;v12];

L1 = length(equal1);
L2 = length(equal2);

E1 = 1:L1;
E2 = 1:L2;
    


figure(1)
%subplot(2,1,1)
plot(t',ev1,'k-o',t',ev2,'r-o','LineWidth',5,'MarkerSize',10)
ax = gca;
set(ax,'fontsize',20)
h = legend('agent 1','agent 2');
%xlabel('Time')
ylabel('evidence','fontsize',25)
set(h,'box','off')
ax.XTick = [];
ax.YTick = [thetam thetap];
ax.YTickLabel = {thetam,thetap};
axis([0 l thetam-.2 thetap+.2])

% figure(2)
% plot(E1',equal1,'b-o',E2',equal2,'g-o','LineWidth',5,'MarkerSize',10)
% ax = gca;
% set(ax,'fontsize',20)
% h = legend('equilibration 1','equilibration 2');
% ylabel('equilibration evidence')
% xlabel('equilibration steps')
% set(h,'box','off')

% 
% subplot(2,1,2)
% plot(t,ag2_social,'-o','Color',[1 .44 0],'LineWidth',5,'MarkerSize',5)
% hold on
% plot(tt,ev2,'-o','Color',[.75 .1 .57],'LineWidth',5,'MarkerSize',5)
% ax = gca;
% set(ax,'fontsize',20)
% xlabel('Time', 'fontsize',25)
% h = legend('Social','Private');
% set(h,'box','off')
% ax.XTick = [0 l];
% ax.YTick = [thetam thetap];
% axis([0 l thetam-.2 thetap+.2])
