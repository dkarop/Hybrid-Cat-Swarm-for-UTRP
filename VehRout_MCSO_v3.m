function [FinalRoutes, TXY, elapsed_time, dzero, done, dtwo, dunsat, FinalFitness, ActiveCats, ATT] = VehRout_MCSO_v3(arguments)
global TT TD nNodes TotAct F1_b1 F1_K1 F1_xm F2_K2 F2_b2 F2_a F3_K3 F3_b3 F4_b4 F4_K4 F4_xm transferPenalty nRoutes...
    SwarmPopulation MaxRouteNodes F2_wa F2_wb F2_wc rid  SMP MR TFV InitMaxLength w1 w2 w3 w4 CDCv RSDv Generations itest elapsed_time dzero done dtwo dunsat ATT

if (nargin==0)
    ttimes = 2;
else
    ttimes = arguments(1,17);
end
  
GlobalBest.Routes(1:nRoutes,1:MaxRouteNodes)=-1;
GlobalBest.Fitness = -100000;
GlobalBest.Length(nRoutes) = 0;
GlobalBest.Fit1 = 0;
GlobalBest.Fit2 = 0;
GlobalBest.Fit3 = 0;
GlobalBest.Fit4 = 0; 
GlobalBest.dzero = 0;
GlobalBest.done = 0;
GlobalBest.dtwo = 0;
GlobalBest.dunsat = 0;
GlobalBest.ATT = 0;
GlobalBest.TT(1:nNodes,1:nNodes) = 0;
GlobalBest.transfers(1:nNodes,1:nNodes) = 0;

%%Test Start
for itest=1:ttimes
    close all
    if (nargin==0) % Modular MR
        MR = 0.4;
        % F4_xm = itest*10;
        % F4_b4 = -(F4_K4/F4_xm); % -K4/xm4 <= b4 <= 0
    else
        if (arguments(1,15) == 1)
            MR = 1 - (itest-1)/(ttimes-1) * (1 - 0.001);
        else
            if (arguments(1,15) == 0)
                MR = 0.001 + (1 - 0.001) * (itest - 1) / (ttimes - 1);
            else
                if (arguments(1,20) == 3)
                    MR = max(rand, 0.01);
                else
                    MR = arguments(1,15);
                end
            end
        end
    end
    tic
    clc

    %% Input Variables
    if (nargin==0)
        nRoutes = 4;
        MaxRouteNodes = 8;
        InitMaxNodes = MaxRouteNodes;
        InitMaxLength = 50;
        KINS = 14;
        SwarmPopulation = 2000;
        InitNodeBias = 0.5;
        F1_K1 = 10;
        F1_xm = 20;
        F1_b1 = -(F1_K1/F1_xm);
        F2_K2 = 10;
        F2_wa = 0.8;
        F2_wb = 0.15;
        F2_wc = 0.05;
        w1 = 2;
        w2 = 1;
        w4 = 1;
         if MaxRouteNodes<6
            w3 = MaxRouteNodes/2000;
         else
            w3 = 1;
         end
        F2_a = 100*F2_wa;
        F2_b2 = 3*F2_K2/F2_a/2;
        F3_K3 = 10;
        F3_b3 = -F3_K3/10; % -K3 <= b3 <= 0
        F4_K4 = 0;
        F4_xm = 20;
        F4_b4 = -(F4_K4/F4_xm); % -K4/xm4 <= b4 <= 0
        transferPenalty = 5;
        Generations = 50;
        SMP = 10;
    else
        nRoutes = arguments(1,1);
        MaxRouteNodes = arguments(1,2);
        InitMaxNodes = MaxRouteNodes;
        InitMaxLength = arguments(1,3);
        KINS = arguments(1,4);
        SwarmPopulation = arguments(1,5);
        InitNodeBias = arguments(1,6);
        F1_K1 = arguments(1,7);
        F1_b1 = -(F1_K1/F1_xm); % -K1/xm <= b1 <= 0
        F1_xm = arguments(1,8);
        F2_K2 = arguments(1,9);
        F3_K3 = arguments(1,10);
        F3_b3 = -F3_K3/10; % -K3 <= b3 <= 0
        transferPenalty = arguments(1,11);
        Generations = arguments(1,12);
        F4_K4 = arguments(1,13);
        F4_xm = arguments(1,14);
        F4_b4 = -(F4_K4/F4_xm); % -K4/xm4 <= b4 <= 0
        SMP = arguments(1,16);
        TFV = arguments(1,18);
        w1 = 1;
        w2 = 1;        
        w4 = 1;
        F2_wa = 0.8;
        F2_wb = 0.15;
        F2_wc = 0.05;
        if MaxRouteNodes<6
            w3 = MaxRouteNodes/2000;
        else
            w3 = 1;
        end
        F2_a = 100*F2_wa;
        F2_b2 = 3*F2_K2/F2_a/2;
        CDCv = arguments(1,19);
        RSDv = arguments(1,20);
    end

    TD = readmatrix('td1.txt');
    TT = readmatrix('tt1.txt');
    TXY = readmatrix('Coordinates1.txt');
    nNodes = length(TD);


    Swarm(SwarmPopulation).Routes(1:nRoutes,1:MaxRouteNodes) = -1; %#ok<AGROW>
    Swarm(SwarmPopulation).Length(1:nRoutes,1) = 0; %#ok<AGROW>
    Swarm(SwarmPopulation).Fitness = -100000; %#ok<AGROW>


    [S, ~]=FloydSPR(TT, nNodes);

    links = zeros(1,nNodes);
    nodeBias = zeros(1,nNodes);
    singlenodes = 0;
    for i=1:nNodes
        for j=1:nNodes
            if(TT(i,j)>0)
                links(i) = links(i) + 1;
            end
        end
        if (links(i)>1)
            nodeBias(i) = 1;
        else
            nodeBias(i) = InitNodeBias;
            singlenodes = singlenodes + 1;
        end
    end
    startnodeBias = nodeBias;
    ActLev=sum(TD);
    TotAct = sum(ActLev);
    [sortActLev(:,1),sortActLev(:,2)] = sort(ActLev,'descend');
    PN = zeros(nRoutes,2);
    INS = zeros(nRoutes,2);
    for i=1:KINS
        INS(i,1) = sortActLev(i,2);
        INS(i,2) = sortActLev(i,1);
    end
    j = 1;
    for i=1:nNodes
        if(nodeBias(i) ~= 1)
            INS(KINS+j,1) = i;
            INS(KINS+j,2) = ActLev(i)*10;
            j = j + 1;
        end
    end
    sumProb=0;
    KINS = length(INS);
    for i=1:KINS
        INS(i,2) = 100*(INS(i,2)+sumProb)/sum(INS(1:KINS ,2));
        sumProb=sumProb+INS(i,2);
    end
    INSstart=INS;
    for iFelis=1:SwarmPopulation
        nodeBias = startnodeBias;
        Swarm(iFelis).Fitness=-100000;
        if (iFelis>1); INS=INSstart; end
        Swarm(iFelis).Routes(1:nRoutes,1:MaxRouteNodes)=-1;
        Swarm(iFelis).Length(1:nRoutes,1) = 0;
        startnodes(1:nRoutes,1) = 0;
        Included(1:nNodes) = 0;
        for i=1:nRoutes
            chooseNode = 100*rand;
            for j=1:KINS
                if (chooseNode<INS(j,2))
                    Swarm(iFelis).Routes(i,1)=INS(j,1);
                    Included(INS(j,1)) = 1;
                    if (j<KINS)
                        INS(j:KINS-1,:)=INS(j+1:KINS,:);
                    end
                    if (i==1); INS(KINS,:)=0; end
                    INS(:,2)=0;
                    sumProb=0;
                    for k=1:KINS
                        if (INS(k,1)>0)
                            INS(k,2) = ActLev(INS(k,1))+sumProb;
                            sumProb=sumProb+ActLev(INS(k,1));
                        end
                    end
                    INS(:,2) = 100*INS(:,2)/sumProb;
                    break
                end
            end
        end
        tempnodeBias = nodeBias;
        VNS = zeros(nNodes,2);
        Cycle=0;
        startnodes(:,1) = Swarm(iFelis).Routes(:,1);
        temproute = zeros(MaxRouteNodes,1);
        reverseroute = zeros(MaxRouteNodes,1);
        for i=1:nRoutes
            for rtest = 1:10
                temproute(1) = startnodes(i,1);
                reversed = 0;
                nodeBias = tempnodeBias;
                PN(i,1)= startnodes(i,1);
                PN(i,2)=1;
                temproute(2:MaxRouteNodes) = -1;
                templength = 0;
                m=1;
                while (m<=InitMaxNodes && Swarm(iFelis).Length(i,1) < InitMaxLength)
                    l=0;
                    sumProb=0;
                    VNS(1:nNodes,:) = 0;
                    for j=1:nNodes
                        if (TT(PN(i,1),j)>0)
                            for k=1:PN(i,2)
                                if (j==temproute(k))
                                    Cycle=1;
                                    break
                                end
                            end
                            if (Cycle==0)
                                l=l+1;
                                VNS(l,1)=j;
                                VNS(l,2)=nodeBias(j)*ActLev(j)+sumProb;
                                sumProb=sumProb+nodeBias(j)*ActLev(j);
                            end
                            Cycle=0;
                        end
                    end
                    VNS(:,2) = 100*VNS(:,2)/sumProb;
                    if (sumProb==0)
                        VNS(1,2)=100;
                    end
                    if (isempty(find(temproute(:) < 0, 1)))
                        break;
                    end
                    chooseNode = 100*rand;
                    if (VNS(1,1)==0)
                        if (reversed == 0)
                            reverseroute(:) = -1;
                            for inode = 1:PN(i,2)
                                reverseroute(PN(i,2)-inode+1) = temproute(inode);
                            end
                            temproute(:) = reverseroute(:);
                            PN(i,1)=temproute(PN(i,2));
                            reversed = 1;
                            continue;
                        else
                            break;
                        end
                    end
                    for j=1:nNodes
                        if (chooseNode<VNS(j,2))
                            temproute(PN(i,2)+1) = VNS(j,1);
                            Included(VNS(j,1)) = 1;
                            templength = templength + TT(PN(i,1),VNS(j,1));
                            PN(i,1)=VNS(j,1);
                            PN(i,2)=PN(i,2)+1;
                            nodeBias(VNS(j,1)) = 0.01*nodeBias(VNS(j,1));
                            break;
                        end
                    end
                    m=m+1;
                end
                if (length(find(temproute < 0)) < length(find(Swarm(iFelis).Routes(i,:) < 0)))
                    Swarm(iFelis).Routes(i,1:MaxRouteNodes) = temproute(1:MaxRouteNodes);
                    Swarm(iFelis).Length(i,1) = templength;
                    rr=100000*rand; if rr<100000*MR; Swarm(iFelis).Flag=0; else; Swarm(iFelis).Flag=1; end
                    if (isempty(find(Swarm(iFelis).Routes(i,:) < 0, 1))); break; end
                end
            end
        end
    end
    rid = fopen('iterations.dat','w');
    iteration = 0;

    [Swarm, ~, GlobalBest] = Evaluate(Swarm, Swarm, S, GlobalBest, iteration);
    %Report(GlobalBest, nRoutes, MaxRouteNodes, nNodes, S, itest);


    if (nargin==0)  % Modular Fitness Threshold
        TFV = 0.5;
    else
        if (arguments(1,18) == 1)
            TFV = 1 - (itest-1)/(ttimes-1) * (1 - 0.3);
        else
            if (arguments(1,18) == 0)
                TFV = 0.3 + (1 - 0.1) * (itest - 1) / (ttimes - 1);
            else
                TFV = arguments(1,18);
            end
        end
    end

    FT = TFV * GlobalBest.Fitness;   % Fitness Threshold

    % Filter members based on fitness
    filteredSwarmIndex = [Swarm.Fitness] >= FT;
    filteredSwarm = Swarm(filteredSwarmIndex);

    % Calculate the number of ones and zeros in the Flag field
    num_ones = sum([filteredSwarm.Flag] == 1);
    num_zeros = length(filteredSwarm) - num_ones;

    % Calculate the number of Active Cats
    ActC(itest)=length(filteredSwarm);%#ok<AGROW>

    fprintf([newline, newline, '\nIter: %d, MR: %.2f, FitΤhresh.: %.2f, GlobalBestFit: %.2f\n', ...
        'Cats %d: %d Seeking, %d Tracing\n'], itest, MR, TFV, GlobalBest.Fitness, length(filteredSwarm), num_ones, num_zeros);

    tcounter = 0;  % Tacing Hit/Miss Counter
    scounter = 0;  % Seeking Hit/Miss Counter

    % Main CSO Loop
    for iteration = 1:Generations
          tempSwarm = filteredSwarm; % Working Cats

        if (nargin==0) % Modular CDCv RSDv
            CDCv = 0.3;
            RSDv = 1;
        else
            if (arguments(1,19) == 1)
                CDCv = 1 - (iteration)/(Generations) * (1 - 0.1);
            else
                if (arguments(1,19) == 0)
                    CDCv = 0.1 + (1 - 0.1) * (iteration) / (Generations);
                else
                    if (arguments(1,19) == 3)
                        CDCv = max(rand, 0.1);
                    else
                        CDCv = arguments(1,19);
                    end
                end
            end
            if (arguments(1,20) == 1)
                RSDv = 1 - (iteration)/(Generations) * (1 - 0.1);
            else
                if (arguments(1,20) == 0)
                    RSDv = 0.1 + (1 - 0.1) * (iteration) / (Generations);
                else
                    if (arguments(1,20) == 3)
                        RSDv = max(rand, 0.1);
                    else
                        RSDv = arguments(1,20);
                    end
                end
            end
        end

        CDC = round(nRoutes * CDCv);
        fprintf('CDCv: %.2f, RSDv: %.2f, F2_a: %.2f, F2_b2: %.2f, F4_xm: %.2f, MR: %.2f\n', CDCv, RSDv, F2_a, F2_b2, F4_xm, MR);


        % Initialize parallel pool (choose an appropriate number of workers)
        %parpool('local', 4);

        for iFelis = 1:numel(tempSwarm)
            if tempSwarm(iFelis).Flag == 1 % Seeking Mode
                smpSwarm = repmat(tempSwarm(iFelis), 1, SMP);

                for s = 1:SMP
                    % fprintf("\ns = %d", s);
                    cRs = randperm(nRoutes, CDC);
                    cR2s = randperm(nRoutes, CDC);
                    for idx = 1:CDC
                        swap = 0;
                        counter = 0;
                        cR = cRs(idx);
                        cR2 = cR2s(idx); 
                        while(swap == 0)
                            invalid = 0;                           
                            k=3;
                            while (k <= MaxRouteNodes-1)
                                l=2;
                                while (l <= MaxRouteNodes-1)
                                    if (smpSwarm(s).Routes(cR, k)==GlobalBest.Routes(cR2, l) && smpSwarm(s).Routes(cR, k)>0)
                                        for i=1:k-1
                                            for j=l+1:MaxRouteNodes
                                                if(smpSwarm(s).Routes(cR,i)==GlobalBest.Routes(cR2, j))
                                                    invalid=1;
                                                    break;
                                                end
                                            end
                                            if (invalid == 1); break; end
                                        end
                                        if (invalid==0)
                                            scounter = scounter + 1;
                                            RSD = round(min([MaxRouteNodes - k, MaxRouteNodes - l]) * RSDv);
                                            smpSwarm(s).Routes(cR, k:k+RSD)=GlobalBest.Routes(cR2, l:l+RSD);
                                            smpSwarm(s).Routes(cR, k+RSD+1:MaxRouteNodes)=-1;
                                            swap = 1;
                                            break;
                                        end
                                    end
                                    l = l + 1;
                                end
                                if (swap == 1); break; end
                                k = k + 1;
                            end
                            counter = counter + 1;
                            if counter == 1 && swap == 0
                                 % Flipping smpSwarm
                                routes1 = smpSwarm(s).Routes(cR, :);
                                routes1(routes1 == -1) = [];
                                routes1 = flip(routes1);
                                smpSwarm(s).Routes(cR, :) = [routes1, -1 * ones(1, MaxRouteNodes - numel(routes1))];                                
                            end
                            if counter == 2 && swap == 0
                                n = 1;
                                routes1 = smpSwarm(s).Routes(cR, :);
                                routes2 = GlobalBest.Routes(cR2, :);
                                routes1(routes1 == -1) = [];
                                routes1 = flip(routes1);
                                routes2(routes2 == -1) = [];
                                 % Trying to fill more nodes to original smpSwarm
                                if numel(routes1) < numel(routes2)
                                    cleanedRoute = [];
                                    while n <= numel(routes2)
                                        if TT(routes1(end), routes2(n)) > 0
                                            temp = [routes1(end), routes2(n:numel(routes2))];
                                            % Iterate through the route in reverse order
                                            for o = numel(temp):-1:1
                                                noden = temp(o);
                                                % Check if the node has not been visited yet
                                                if ~ismember(noden, cleanedRoute)
                                                    % Add the node to the cleaned route
                                                    LeanCat = [noden, cleanedRoute];
                                                end
                                            end
                                            isDuplicate = any(diff(sort(LeanCat(LeanCat ~= -1))) == 0);
                                            if isDuplicate == true; break; end
                                            smpSwarm(s).Routes(cR, :) = [LeanCat, -1 * ones(1, MaxRouteNodes - numel(LeanCat))];
                                            break
                                        end
                                        n = n + 1;
                                        if n == numel(routes2); break; end
                                    end
                                else
                                    break;
                                end
                            end
                            if counter == 3 && swap == 0
                                n = 1;
                                routes1 = smpSwarm(s).Routes(cR, :);
                                routes2 = GlobalBest.Routes(cR2, :);
                                routes1(routes1 == -1) = [];
                                routes1 = flip(routes1);
                                routes2(routes2 == -1) = [];
                                 % Trying once more to fill more nodes to flipped smpSwarm
                                if numel(routes1) <= numel(routes2)
                                    cleanedRoute = [];
                                    while n < numel(routes2)
                                        if TT(routes1(end), routes2(n)) > 0
                                            temp = [routes1(end), routes2(n:numel(routes2))];
                                            % Iterate through the route in reverse order
                                            for o = numel(temp):-1:1
                                                noden = temp(o);
                                                % Check if the node has not been visited yet
                                                if ~ismember(noden, cleanedRoute)
                                                    % Add the node to the cleaned route
                                                    LeanCat = [noden, cleanedRoute];
                                                end
                                            end
                                            isDuplicate = any(diff(sort(LeanCat(LeanCat ~= -1))) == 0);
                                            if isDuplicate == true; break; end
                                            smpSwarm(s).Routes(cR, :) = [LeanCat, -1 * ones(1, MaxRouteNodes - numel(LeanCat))];
                                            break
                                        end
                                        n = n + 1;
                                        if n == numel(routes2); break; end
                                    end
                                else
                                    break;
                                end
                            end

                            if (swap == 1 || counter == 4); break; end
                        end  %While swap == 0
                    end
                    for i = 1:nRoutes
                        route = smpSwarm(s).Routes(i, :);
                        route(route == -1) = [];
                        newRoute = HillClimb(route, TT);
                        newRoute2 = [newRoute, -1 * ones(1, MaxRouteNodes - numel(newRoute))];%#ok<AGROW>
                        smpSwarm(s).Routes(i, :) = newRoute2;
                    end
                end
                [tempSwarm(iFelis), GlobalBest] = EvaluateAgent(smpSwarm, tempSwarm(iFelis), S, GlobalBest);
            else % Tracing Mode
                swap = 0;
                cnt = 0;

                while swap == 0
                    cR = ceil(nRoutes * rand);
                    cR2 = ceil(nRoutes * rand);
                     % Route exchange smpSwarm with GlobalBest
                    if isequal(tempSwarm(iFelis).Routes(cR, :), GlobalBest.Routes(cR2, :))
                        break;
                    else
                        tcounter = tcounter + 1;
                        tempSwarm(iFelis).Routes(cR, :) = GlobalBest.Routes(cR2, :);
                        swap = 1;
                    end

                    if swap == 1
                        break;
                    end

                    cnt = cnt + 1;

                    if cnt == MaxRouteNodes
                        break;
                    end
                end
            end
        end
        [filteredSwarm, ~, GlobalBest] = Evaluate(tempSwarm, filteredSwarm, S, GlobalBest, iteration);
        FF(itest,iteration) = GlobalBest.Fitness;%#ok<AGROW>
    end
    %delete(gcp);
    fclose(rid);
    GlobalBest.Routes
    FinalRoutes = zeros(nRoutes, MaxRouteNodes, 2);
    for i=1:nRoutes
        for j=1:MaxRouteNodes
            if (GlobalBest.Routes(i,j)>0)
                FinalRoutes(i,j,1) = TXY(GlobalBest.Routes(i,j),1);
                FinalRoutes(i,j,2) = TXY(GlobalBest.Routes(i,j),2);
            else
                FinalRoutes(i,j,1) = FinalRoutes(i,j-1,1);
                FinalRoutes(i,j,2) = FinalRoutes(i,j-1,2);
            end
        end
    end
  Report(GlobalBest, nRoutes, MaxRouteNodes, nNodes, S, itest);
end %for itest

tempFF = FF(:);
FinalFitness = [0; tempFF];
ActivCats = ActC(:);
ActiveCats = [SwarmPopulation; ActivCats];
fprintf('\n\nSeeking Hits: %d, Seeking Misses: %d\n', round((scounter/(SMP*CDC))), round((Generations*SwarmPopulation-(scounter/(SMP*CDC)))));
fprintf('Tracing Hits: %d, Tracing Misses: %d\n\n', round(tcounter), round((Generations*SwarmPopulation-tcounter)));

Results_Output(itest, nRoutes);
if (nargin==0)
    figure;
    plotRoutes(nRoutes, FinalRoutes, TXY)
    clear
end
end %Main Function end

function [Felis, GlobalBest] = EvaluateAgent(smpSwarm, Felis, S, GlobalBest)
global TT TD nNodes TotAct F1_b1 F1_K1 F1_xm F2_K2 F2_b2 F2_a F3_K3 F3_b3 F4_b4 F4_K4 F4_xm transferPenalty nRoutes...
    MaxRouteNodes F2_wa F2_wb F2_wc w1 w2 w3 w4 rid
Fitness = zeros(1,length(smpSwarm));
transfersRoutes = cell(nNodes,nNodes);
MinRoutes = zeros(MaxRouteNodes,nRoutes);
tempTT = zeros(nNodes,nNodes);
transfers = zeros(nNodes,nNodes);
trans=zeros(1,nRoutes);
for iFelis=1:length(smpSwarm)
    ATT = 0;
    transfersRoutes(:,:) = {zeros(1,nRoutes)};
    tempTT(1:nNodes,1:nNodes) = 0;
    Fit1_num = 0; Fit1_denom = 0; dzero = 0; done = 0; dtwo = 0; dunsat = 0;
    smpSwarm(iFelis).Length(:) = 0;
    for i=1:nRoutes
        LeanCat = smpSwarm(iFelis).Routes(i,:);
        isDuplicate = any(diff(sort(LeanCat(LeanCat ~= -1))) == 0);
        if isDuplicate == true; fprintf("\nAgent with Duplicate Nodes Found!\n"); pause(1); break; end
        for j=2:MaxRouteNodes
            if (smpSwarm(iFelis).Routes(i,j)>0)
                smpSwarm(iFelis).Length(i) = smpSwarm(iFelis).Length(i) + ...
                    TT(smpSwarm(iFelis).Routes(i,j-1),smpSwarm(iFelis).Routes(i,j));
                x = smpSwarm(iFelis).Routes(i,j);
                y = smpSwarm(iFelis).Routes(i,j-1);
                tempTT(x,y) = TT(x,y);
                tempTT(y,x) = TT(y,x);
                transfersRoutes{x,y}(i) = 1;
                transfersRoutes{y,x}(i) = 1;
            else
                break;
            end
        end
    end
    [tempS, tempP]=FloydSPR(tempTT, nNodes);
    for i = 1:nNodes-1
        for j = i+1:nNodes
            MinRoutes(1:MaxRouteNodes,1:nRoutes) = 0;
            if (tempS(i,j) == inf)
                transfers(i,j) = -1;
                transfers(j,i) = transfers(i,j);
                continue;
            end
            index = 1;
            previousk = i;
            k = tempP(i, j);
            if (k==j)
                continue;
            end
            while(previousk ~= j)
                MinRoutes(index,:) = transfersRoutes{previousk,k}(:);
                index = index + 1;
                previousk = k;
                k = tempP(k,j);
            end
            stop = 1;
            RouteEnd = sum(MinRoutes(1,:));
            while (RouteEnd ~= 0)
                trans(1:nRoutes)=0;
                for k = 1:nRoutes
                    for l = stop:MaxRouteNodes
                        trans(k) = trans(k) + MinRoutes(l,k);
                        if(MinRoutes(l,k)==0)
                            break;
                        end
                    end
                end
                stop = stop + max(trans);
                if (stop > MaxRouteNodes)
                    transfers(i,j) = -1;
                    transfers(j,i) = transfers(i,j);
                    break;
                end
                RouteEnd = sum(MinRoutes(stop,:));
                if (RouteEnd ~= 0)
                    transfers(i,j) = transfers(i,j) + 1;
                    transfers(j,i) = transfers(i,j);
                end
            end
            tempS(i,j) = tempS(i,j) + transferPenalty*transfers(i,j);
            tempS(j,i) = tempS(i,j);
        end
    end
    for i = 1:nNodes-1
        for j = i+1:nNodes
            x = tempS(i,j) - S(i,j);
            if (x<=F1_xm)
                tempf = -(F1_b1/F1_xm + F1_K1/F1_xm^2)*x^2 + F1_b1*x + F1_K1;
                Fit1_num = Fit1_num + tempf*TD(i,j);
                if (tempf>0)
                    Fit1_denom = Fit1_denom + TD(i,j);
                end
            end
            switch transfers(i,j)
                case 0
                    dzero = dzero + 2*100*TD(i,j)/TotAct;
                case 1
                    done = done + 2*100*TD(i,j)/TotAct;
                case 2
                    dtwo = dtwo + 2*100*TD(i,j)/TotAct;
                otherwise
                    dunsat = dunsat + 2*100*TD(i,j)/TotAct;
            end
            ATT = ATT + 2*TD(i,j)*tempS(i,j);
        end
    end
    Fit1 = Fit1_num/Fit1_denom;
    dt = F2_wa*dzero + F2_wb*done + F2_wc*dtwo;
    Fit2 = (F2_K2 - F2_b2*F2_a)/F2_a/F2_a*dt^2 + F2_b2*dt;
    Fit3 = -(F3_b3 + F3_K3)*dunsat^2 + F3_b3*dunsat + F3_K3;
    lentot = sum(smpSwarm(iFelis).Length(:)) - F4_xm;
    Fit4 = -(F4_b4/F4_xm + F4_K4/F4_xm^2)*lentot^2 + F4_b4*lentot + F4_K4;
    if (Fit4 > F4_K4); Fit4 = F4_K4; end
    Fitness(iFelis) = w1*Fit1 + w2*Fit2 + w3*Fit3 + w4*Fit4;
    smpSwarm(iFelis).Fitness = Fitness(iFelis);
    if (Fitness(iFelis) > Felis.Fitness)
        Felis.Fitness = Fitness(iFelis);
        Felis.Length(:) = smpSwarm(iFelis).Length(:);
        Felis.Routes(1:nRoutes,1:MaxRouteNodes) = smpSwarm(iFelis).Routes(1:nRoutes,1:MaxRouteNodes);
    end
    if (Fitness(iFelis) > GlobalBest.Fitness)
        GlobalBest.Fitness = Fitness(iFelis);
        GlobalBest.Length(:) = smpSwarm(iFelis).Length(:);
        GlobalBest.Fit1 = Fit1;
        GlobalBest.Fit2 = Fit2;
        GlobalBest.Fit3 = Fit3;
        GlobalBest.Fit4 = Fit4;
        GlobalBest.ATT = ATT;
        GlobalBest.dzero = dzero;
        GlobalBest.done = done;
        GlobalBest.dtwo = dtwo;
        GlobalBest.dunsat = dunsat;
        GlobalBest.Routes(1:nRoutes,1:MaxRouteNodes) = smpSwarm(iFelis).Routes(1:nRoutes,1:MaxRouteNodes);
        GlobalBest.TT = tempS;
        GlobalBest.transfers = transfers;
        disp(newline);
    end
end
% [best, ibest] = max(Fitness);
% fprintf('\n%i\t\t%d\t\t%d\n', iteration, best, GlobalBest.Fitness);
%fprintf(rid, '%e\t%e\t%e\t%e\t%e\t%e\n', GlobalBest.Fitness, GlobalBest.dzero, GlobalBest.done, GlobalBest.dtwo, GlobalBest.ATT, sum(GlobalBest.Length(:)));
end




function [Swarm, ibest, GlobalBest] = Evaluate(tempSwarm, Swarm, S, GlobalBest, iteration)
%fprintf("\nEvaluate_summoned\n")
global TT TD nNodes TotAct F1_b1 F1_K1 F1_xm F2_K2 F2_b2 F2_a F3_K3 F3_b3 F4_b4 F4_K4 F4_xm transferPenalty nRoutes...
    SwarmPopulation MaxRouteNodes F2_wa F2_wb F2_wc rid w1 w2 w3 w4
Fitness = zeros(1,SwarmPopulation);
transfersRoutes = cell(nNodes,nNodes);
MinRoutes = zeros(MaxRouteNodes,nRoutes);
tempTT = zeros(nNodes,nNodes);
transfers = zeros(nNodes,nNodes);
trans=zeros(1,nRoutes);
for iFelis=1:length(Swarm)
    ATT = 0;
    transfersRoutes(:,:) = {zeros(1,nRoutes)};
    tempTT(1:nNodes,1:nNodes) = 0;
    transfers(1:nNodes,1:nNodes) = 0;
    Fit1_num = 0; Fit1_denom = 0; dzero = 0; done = 0; dtwo = 0; dunsat = 0;
    tempSwarm(iFelis).Length(:) = 0;
    for i=1:nRoutes
        LeanCat = tempSwarm(iFelis).Routes(i,:);
        isDuplicate = any(diff(sort(LeanCat(LeanCat ~= -1))) == 0);
        if isDuplicate == true; fprintf("\nSwarm Member with Duplicate Nodes Found!\n"); pause(1); break; end
        for j=2:MaxRouteNodes            
            if (tempSwarm(iFelis).Routes(i,j)>0)
                tempSwarm(iFelis).Length(i) = tempSwarm(iFelis).Length(i) + ...
                    TT(tempSwarm(iFelis).Routes(i,j-1),tempSwarm(iFelis).Routes(i,j));
                x = tempSwarm(iFelis).Routes(i,j);
                y = tempSwarm(iFelis).Routes(i,j-1);
                tempTT(x,y) = TT(x,y);
                tempTT(y,x) = TT(y,x);
                transfersRoutes{x,y}(i) = 1;
                transfersRoutes{y,x}(i) = 1;
            else
                break;
            end
        end
    end
    [tempS, tempP]=FloydSPR(tempTT, nNodes);
    for i = 1:nNodes-1
        for j = i+1:nNodes
            MinRoutes(1:MaxRouteNodes,1:nRoutes) = 0;
            if (tempS(i,j) == inf)
                transfers(i,j) = -1;
                transfers(j,i) = transfers(i,j);
                continue;
            end
            index = 1;
            previousk = i;
            k = tempP(i, j);
            if (k==j)
                continue;
            end
            while(previousk ~= j)
                MinRoutes(index,:) = transfersRoutes{previousk,k}(:);
                index = index + 1;
                previousk = k;
                k = tempP(k,j);
            end
            stop = 1;
            RouteEnd = sum(MinRoutes(1,:));
            while (RouteEnd ~= 0)
                trans(1:nRoutes)=0;
                for k = 1:nRoutes
                    for l = stop:MaxRouteNodes
                        trans(k) = trans(k) + MinRoutes(l,k);
                        if(MinRoutes(l,k)==0)
                            break;
                        end
                    end
                end
                stop = stop + max(trans);
                if (stop > MaxRouteNodes)
                    transfers(i,j) = -1;
                    transfers(j,i) = transfers(i,j);
                    break;
                end
                RouteEnd = sum(MinRoutes(stop,:));
                if (RouteEnd ~= 0)
                    transfers(i,j) = transfers(i,j) + 1;
                    transfers(j,i) = transfers(i,j);
                end
            end
            tempS(i,j) = tempS(i,j) + transferPenalty*transfers(i,j);
            tempS(j,i) = tempS(i,j);
        end
    end
    for i = 1:nNodes-1
        for j = i+1:nNodes
            x = tempS(i,j) - S(i,j);
            if (x<=F1_xm)
                tempf = -(F1_b1/F1_xm + F1_K1/F1_xm^2)*x^2 + F1_b1*x + F1_K1;
                Fit1_num = Fit1_num + tempf*TD(i,j);
                if (tempf>0)
                    Fit1_denom = Fit1_denom + TD(i,j);
                end
            end
            switch transfers(i,j)
                case 0
                    dzero = dzero + 2*100*TD(i,j)/TotAct;
                case 1
                    done = done + 2*100*TD(i,j)/TotAct;
                case 2
                    dtwo = dtwo + 2*100*TD(i,j)/TotAct;
                otherwise
                    dunsat = dunsat + 2*100*TD(i,j)/TotAct;
            end
            ATT = ATT + 2*TD(i,j)*tempS(i,j);
        end
    end
    Fit1 = Fit1_num/Fit1_denom;
    dt = F2_wa*dzero + F2_wb*done + F2_wc*dtwo;
    Fit2 = (F2_K2 - F2_b2*F2_a)/F2_a/F2_a*dt^2 + F2_b2*dt;
    Fit3 = -(F3_b3 + F3_K3)*dunsat^2 + F3_b3*dunsat + F3_K3;
    lentot = sum(tempSwarm(iFelis).Length(:)) - F4_xm;
    Fit4 = -(F4_b4/F4_xm + F4_K4/F4_xm^2)*lentot^2 + F4_b4*lentot + F4_K4;
    if (Fit4 > F4_K4); Fit4 = F4_K4; end
    Fitness(iFelis) = w1*Fit1 + w2*Fit2 + w3*Fit3 + w4*Fit4;
    ATT = ATT/TotAct;
    tempSwarm(iFelis).Fitness = Fitness(iFelis);
    if (Fitness(iFelis) > Swarm(iFelis).Fitness)
        Swarm(iFelis).Fitness = Fitness(iFelis);
        Swarm(iFelis).Length(:) = tempSwarm(iFelis).Length(:);
        Swarm(iFelis).Routes(1:nRoutes,1:MaxRouteNodes) = tempSwarm(iFelis).Routes(1:nRoutes,1:MaxRouteNodes);
    end
    if (Fitness(iFelis) > GlobalBest.Fitness)
        GlobalBest.Fitness = Fitness(iFelis);
        GlobalBest.Length(:) = tempSwarm(iFelis).Length(:);
        GlobalBest.Fit1 = Fit1;
        GlobalBest.Fit2 = Fit2;
        GlobalBest.Fit3 = Fit3;
        GlobalBest.Fit4 = Fit4;
        GlobalBest.ATT = ATT;
        GlobalBest.dzero = dzero;
        GlobalBest.done = done;
        GlobalBest.dtwo = dtwo;
        GlobalBest.dunsat = dunsat;
        GlobalBest.Routes(1:nRoutes,1:MaxRouteNodes) = tempSwarm(iFelis).Routes(1:nRoutes,1:MaxRouteNodes);
        GlobalBest.TT = tempS;
        GlobalBest.transfers = transfers;
        disp(newline);
    end
end
[best, ibest] = max(Fitness);
fprintf('\n%i\t\t%d\t\t%d\n', iteration, best, GlobalBest.Fitness);
fprintf(rid, '%e\t%e\t%e\t%e\t%e\t%e\n', GlobalBest.Fitness, GlobalBest.dzero, GlobalBest.done, GlobalBest.dtwo, GlobalBest.ATT, sum(GlobalBest.Length(:)));
end


function [S, P]=FloydSPR(AdjMax, N)
P=-1*ones(N,N);
S=AdjMax;
for k=1:N-1
    for i=k+1:N
        if (S(i,k)<=0)
            S(i,k)=inf;
            S(k,i)=inf;
        else
            P(i,k)=k;
            P(k,i)=i;
        end
    end
end
for k=1:N
    for i=1:N
        if (i ~= k)
            for j=1:N
                if (i ~= j && k ~= j)
                    if S(i,k)==inf
                        continue;
                    end
                    if S(k,j)==inf
                        continue;
                    end
                    if S(i,j)>S(i,k)+S(k,j)
                        if P(i,k)==-1
                            P(i,j)=k;
                        else
                            P(i,j)=P(i,k);
                        end
                        S(i,j)=S(i,k)+S(k,j);
                    end
                end
            end
        end
    end
end
end


function Report(GlobalBest, nRoutes, MaxRouteNodes, nNodes, S, itest)
results_file = strcat('Results_', int2str(nRoutes), '_', int2str(itest), '.dat');
fid = fopen(results_file,'w');
fprintf(fid, '\n%e', GlobalBest.Fitness);
fprintf(fid, '\n\n%e\t\t%e\t\t%e\t\t%e', GlobalBest.Fit1, GlobalBest.Fit2, GlobalBest.Fit3, GlobalBest.Fit4);
fprintf(fid, '\n\n%e\t\t%e\t\t%e\t\t%e', GlobalBest.dzero, GlobalBest.done, GlobalBest.dtwo, GlobalBest.dunsat);
fprintf(fid, '\n\n%e min\n\n', GlobalBest.ATT);
for i = 1:nRoutes
    for j = 1:MaxRouteNodes
        fprintf(fid, '%d\t', GlobalBest.Routes(i,j));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n\n');
for i = 1:nNodes
    for j = 1:nNodes
        fprintf(fid, '%d\t', GlobalBest.TT(i,j));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n\n');
for i = 1:nNodes
    for j = 1:nNodes
        fprintf(fid, '%d\t', GlobalBest.transfers(i,j));
    end
    fprintf(fid, '\n');
end
fprintf('\nMethod A');
fprintf('\n%e', GlobalBest.Fitness);
fprintf('\n%e\t\t%e\t\t%e\t\t%e', GlobalBest.Fit1, GlobalBest.Fit2, GlobalBest.Fit3, GlobalBest.Fit4);
fprintf('\n%e\t\t%e\t\t%e\t\t%e', GlobalBest.dzero, GlobalBest.done, GlobalBest.dtwo, GlobalBest.dunsat);
fprintf('\n%e min', GlobalBest.ATT);
EvaluateBest(GlobalBest, S, fid);
fclose(fid);
end



function EvaluateBest(GlobalBest, S, fid)
global TT TD nNodes TotAct F1_b1 F1_K1 F1_xm F2_K2 F2_b2 F2_a F3_K3 F3_b3 F4_b4 F4_K4 F4_xm transferPenalty nRoutes...
    MaxRouteNodes F2_wa F2_wb F2_wc dzero done dtwo dunsat w1 w2 w3 w4 elapsed_time ATT
transfersRoutes = cell(nNodes,nNodes);
transfersRoutes(:,:) = {zeros(1,nRoutes)};
tempTT = zeros(nNodes,nNodes);
tempP=-1*ones(nNodes,nNodes);
testTT = zeros(nNodes,nNodes);
testP=-1*ones(nNodes,nNodes);
ATT = 0;
tempTT(1:nNodes,1:nNodes) = inf;
Fit1_num = 0; Fit1_denom = 0; dzero = 0; done = 0; dtwo = 0; dunsat = 0;
for i=1:nRoutes
    testTT(1:nNodes,1:nNodes) = 0;
    for j=2:MaxRouteNodes
        if (GlobalBest.Routes(i,j)>0)
            x = GlobalBest.Routes(i,j-1);
            y = GlobalBest.Routes(i,j);
            testTT(x,y) = TT(x,y);
            testTT(y,x) = TT(y,x);
            tempTT(x,y) = TT(x,y);
            tempTT(y,x) = TT(y,x);
            k = j-1;
            TTprev = TT(x,y);
            transfersRoutes{x,y}(i) = 1;
            transfersRoutes{y,x}(i) = 1;
            testP(x,y) = y;
            testP(y,x) = x;
            tempP(x,y) = y;
            tempP(y,x) = x;
            while(k-1>0)
                x = GlobalBest.Routes(i,k-1);
                w = GlobalBest.Routes(i,k);
                z = GlobalBest.Routes(i,j-1);
                testTT(x,y) = testTT(x,z) + TTprev;
                testTT(y,x) = testTT(z,x) + TTprev;
                testP(x,y) = w;
                testP(y,x) = z;
                k = k - 1;
                if (testTT(x,y) < tempTT(x,y))
                    tempTT(x,y) = testTT(x,y);
                    tempTT(y,x) = testTT(y,x);
                    tempP(x,y) = testP(x,y);
                    tempP(y,x) = testP(y,x);
                end
            end
        else
            break;
        end
    end
end
[tempS, tempP, transfers]=FloydSPR_2(tempTT, tempP, nNodes, transfersRoutes, nRoutes, MaxRouteNodes, transferPenalty);
for i = 1:nNodes-1
    for j = i+1:nNodes
        x = tempS(i,j) - S(i,j);
        if (x<=F1_xm)
            tempf = -(F1_b1/F1_xm + F1_K1/F1_xm^2)*x^2 + F1_b1*x + F1_K1;
            Fit1_num = Fit1_num + tempf*TD(i,j);
            if (tempf>0)
                Fit1_denom = Fit1_denom + TD(i,j);
            end
        end
        switch transfers(i,j)
            case 0
                dzero = dzero + 2*100*TD(i,j)/TotAct;
            case 1
                done = done + 2*100*TD(i,j)/TotAct;
            case 2
                dtwo = dtwo + 2*100*TD(i,j)/TotAct;
            otherwise
                dunsat = dunsat + 2*100*TD(i,j)/TotAct;
        end
        ATT = ATT + 2*TD(i,j)*tempS(i,j);
    end
end
Fit1 = Fit1_num/Fit1_denom;
dt = F2_wa*dzero + F2_wb*done + F2_wc*dtwo;
Fit2 = (F2_K2 - F2_b2*F2_a)/F2_a/F2_a*dt^2 + F2_b2*dt;
Fit3 = -(F3_b3 + F3_K3)*dunsat^2 + F3_b3*dunsat + F3_K3;
lentot = sum(GlobalBest.Length(:)) - F4_xm;
Fit4 = -(F4_b4/F4_xm + F4_K4/F4_xm^2)*lentot^2 + F4_b4*lentot + F4_K4;
if (Fit4 > F4_K4); Fit4 = F4_K4; end
Fitness = w1*Fit1 + w2*Fit2 + w3*Fit3 + w4*Fit4;
ATT = ATT/TotAct;
elapsed_time = toc;
fprintf('\n===========================================================================\nMethod B');
fprintf('\n%.2f', elapsed_time);
fprintf('\n%.2f', Fitness);
fprintf('\n%.2f\t\t%.2f\t\t%.2f\t\t%.2f', Fit1, Fit2, Fit3, Fit4);
fprintf('\n%.2f\t\t%.2f\t\t%.2f\t\t%.2f', dzero, done, dtwo, dunsat);
fprintf('\n%.2f min', ATT);
fprintf(fid, '\n%e seconds', elapsed_time);
fprintf(fid, '\n%e', Fitness);
fprintf(fid, '\n\n%e\t\t%e\t\t%e\t\t%e', Fit1, Fit2, Fit3, Fit4);
fprintf(fid, '\n\n%e\t\t%e\t\t%e\t\t%e', dzero, done, dtwo, dunsat);
fprintf(fid, '\n\n%e min\n\n', ATT);
fprintf(fid, '\n\n%e\n\n', sum(GlobalBest.Length(:)));
fprintf(fid, '\n\n');
for i = 1:nNodes
    for j = 1:nNodes
        fprintf(fid, '%d\t', tempS(i,j));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n\n');
for i = 1:nNodes
    for j = 1:nNodes
        fprintf(fid, '%d\t', transfers(i,j));
    end
    fprintf(fid, '\n');
end
end



function [S, P, transfers]=FloydSPR_2(AdjMax, P, N, transfersRoutes, nRoutes, MaxRouteNodes, transferPenalty)
MinRoutes = zeros(MaxRouteNodes,nRoutes);
transfers = zeros(N,N);
trans=zeros(1,nRoutes);
S=AdjMax;
for i=1:N
    S(i,i) = 0;
end
for k=1:N
    for i=1:N
        if (i ~= k)
            for j=1:N
                if (i ~= j && k ~= j)
                    if S(i,k)==inf
                        continue;
                    end
                    if S(k,j)==inf
                        continue;
                    end
                    if (S(i,j)<S(i,k)+S(k,j))
                        continue;
                    end
                    MinRoutes(1:MaxRouteNodes,1:nRoutes) = 0;
                    changes = 0;
                    backtrack = 0;
                    index = 1;
                    previous = i;
                    finish = k;
                    test = P(i, finish);
                    if (test ~= -1)
                        next = P(test, finish);
                    end
                    if (next == previous)
                        continue;
                    end
                    while(previous ~= j)
                        MinRoutes(index,:) = transfersRoutes{previous,test}(:);
                        index = index + 1;
                        previous = test;
                        if (test==k)
                            finish=j;
                        end
                        test = P(test,finish);
                        if (test ~= -1)
                            if (test == k)
                                next = P(k, j);
                            else
                                next = P(test, finish);
                            end
                        end
                        if (next == previous)
                            backtrack = 1;
                            break;
                        end
                    end
                    if (backtrack == 1)
                        continue;
                    end
                    stop = 1;
                    RouteEnd = sum(MinRoutes(1,:));
                    while (RouteEnd ~= 0)
                        trans(1:nRoutes)=0;
                        for m = 1:nRoutes
                            for l = stop:MaxRouteNodes
                                trans(m) = trans(m) + MinRoutes(l,m);
                                if(MinRoutes(l,m)==0)
                                    break;
                                end
                            end
                        end
                        stop = stop + max(trans);
                        if (stop > MaxRouteNodes)
                            transfers(i,j) = -1;
                            transfers(j,i) = transfers(i,j);
                            break;
                        end
                        RouteEnd = sum(MinRoutes(stop,:));
                        if (RouteEnd ~= 0)
                            changes = changes + 1;
                        end
                    end
                    current_changes = changes - transfers(i,k) - transfers(k,j);
                    if (S(i,j)>S(i,k)+S(k,j) + transferPenalty*current_changes)
                        if P(i,k)==-1
                            P(i,j)=k;
                        else
                            P(i,j)=P(i,k);
                        end
                        S(i,j)=S(i,k)+S(k,j) + transferPenalty*current_changes;
                        transfers(i,j) = changes;
                    end
                end
            end
        end
    end
end
for i=1:N
    for j=1:N
        if (S(i,j)==inf)
            transfers(i,j)=-1;
        end
    end
end
end

% function newRoute = HillClimb(route, TD)
% newRoute = route;
% 
% % Randomly select two distinct indices to swap
% indices = randperm(length(route), 2);
% idx1 = indices(1);
% idx2 = indices(2);
% 
% % Swap the elements at idx1 and idx2
% temp = newRoute(idx1);
% newRoute(idx1) = newRoute(idx2);
% newRoute(idx2) = temp;
% 
% totalDistance = 0;
% for i = 2:length(route)
%     fromNode = newRoute(i - 1);
%     toNode = newRoute(i);
%     totalDistance = totalDistance + TD(fromNode, toNode);
% end
% % If the new route is better and does not loop, keep it; otherwise, revert to the original route
% if totalDistance < CalculateRouteDistance(route, TD, length(route))
%     for i = 1:length(route)
%         for j = i + 1:length(route)
%             if route(i) == route(j)
%                 break;
%             else
%                 newRoute = route;
%             end
%         end
%     end
% end
% end
% 
% function distance = CalculateRouteDistance(route, TD, RouteNodes)
% distance = 0;
% for i = 2:RouteNodes
%     if route(i) == -1
%         break;
%     end
%     fromNode = route(i - 1);
%     toNode = route(i);
%     distance = distance + TD(fromNode, toNode);
% end
% end

function improvedRoute = HillClimb(route, TT)
global bestNeighborRouteTime currentRouteTime
    % Initialize variables
    currentRoute = route;
    currentRouteTime = calculateRouteTime(currentRoute, TT);

    % Perform Hill Climbing iterations
    maxIterations = 100;
    convergenceThreshold = 1;
    iteration = 1;

     while iteration <= maxIterations
        % Find pairs of nodes that can potentially fold the route
        foldCandidates = findFoldCandidates(currentRoute, TT);

        % If foldable pairs are found, explore them
        if ~isempty(foldCandidates)
            bestNeighborRoute = currentRoute;
            bestNeighborRouteTime = currentRouteTime;

            % Iterate through foldable pairs
            for foldIdx = 1:length(foldCandidates)
                splitIdx = foldCandidates(foldIdx);

                % Split the route into two segments
                segment1 = currentRoute(1:splitIdx);
                segment2 = currentRoute(splitIdx+1:end);

                % Combine segments in reverse order
                neighborRoute = [segment2, segment1];

                % Calculate the time of the neighboring solution
                neighborRouteTime = calculateRouteTime(neighborRoute, TT);

                % Check if the neighboring solution is better
                if neighborRouteTime < bestNeighborRouteTime
                    bestNeighborRoute = neighborRoute;
                    bestNeighborRouteTime = neighborRouteTime;
                end
            end

            % Update current route if a better solution is found
            if bestNeighborRouteTime < currentRouteTime
                currentRoute = bestNeighborRoute;
                currentRouteTime = bestNeighborRouteTime;
            end
        end
        % Check for convergence
        if currentRouteTime - bestNeighborRouteTime < convergenceThreshold
            break; % No significant improvement
        end

        iteration = iteration + 1;
    end

    % Return the improved route
    improvedRoute = currentRoute;
    if improvedRoute ~= route
        fprintf("HC Improved Route!");
    end
end

% Helper function to calculate route time based on TT matrix
function routeTime = calculateRouteTime(route, TT)
global routeTime
for k = 1:(length(route) - 1)
    routeTime = routeTime + TT(route(k), route(k + 1));
end
end

% Helper function to find pairs of nodes that can potentially fold the route
function foldCandidates = findFoldCandidates(route, TT)
    foldCandidates = [];
    for i = 2:length(route) - 1 % Exclude the first and last nodes
        node1 = route(i - 1);
        node2 = route(i + 1);
        if TT(node1, node2) > 0
            foldCandidates = [foldCandidates, i]; %#ok<AGROW>
        end
    end
end




function Results_Output(maxtest, nRoutes)
%clc
%clear
results_file_out = strcat('Results_', int2str(nRoutes), '.dat');
fid = fopen(results_file_out,'w');
for test = 1:maxtest
    results_file_in = strcat('Results_', int2str(nRoutes), '_', int2str(test), '.dat');
    fid2 = fopen(results_file_in,'r');
    %%
    Method_A = fscanf(fid2, '%e');
    for i = 1:37+nRoutes
        line = fgetl(fid2); %#ok<NASGU>
    end
    elapsed_time = fscanf(fid2, '%e');
    line = fgetl(fid2); %#ok<NASGU>
    %%
    Method_B = fscanf(fid2, '%e');
    line = fgetl(fid2); %#ok<NASGU>
    ltot = fscanf(fid2, '%e');
    for i = 1:length(Method_A)
        fprintf(fid, '%e\t', Method_A(i));
    end
    for i = 1:length(Method_B)
        fprintf(fid, '%e\t', Method_B(i));
    end
    %%
    fprintf(fid, '%e\t', ltot(1));
    fprintf(fid, '%e\t', elapsed_time);
    fprintf(fid, '\n');
    fclose(fid2);
end
fclose('all');
clear
end

function plotRoutes(nRoutes, FinalRoutes, TXY)
hold 'on'
marker={'--r','--b','--y','--g','--m','--c','--k','-.r'}; %[.49 1 .63]
for i =1:nRoutes
    plot(FinalRoutes(i,:,1),FinalRoutes(i,:,2), marker{i}, 'LineWidth',1.0+(nRoutes-i)/2);
    leg(i,:)=['Route ',num2str(i)]; %#ok<AGROW>
end
plot(TXY(:,1),TXY(:,2), 'rs', 'MarkerSize',10);
leg(i+1,:)='Nodes  ';
xlim([0 100])
ylim([0 100])
xlabel('X')
ylabel('Y')
legend(leg, 'Location', 'best')
grid
hold 'off'
end