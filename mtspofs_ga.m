function varargout = mtspofs_ga(dmat,salesmen,min_tour,pop_size,num_iter,show_prog,show_res)
% MTSPOFS_GA Fixed Start Open Multiple Traveling Salesmen Problem (M-TSP) Genetic Algorithm (GA)
%   Finds a (near) optimal solution to a variation of the "open" M-TSP by
%   setting up a GA to search for the shortest route (least distance needed
%   for each salesman to travel from the start location to unique individual
%   cities without returning to the starting location)
%
% Summary:
%     1. Each salesman starts at the first point, but travels to a unique
%        set of cities after that (and none of them close their loops by
%        returning to their starting points)
%     2. Except for the first, each city is visited by exactly one salesman
%
% Note: The Fixed Start is taken to be the first XY point
%
% Input:
%     DMAT (float) is an NxN matrix of city-to-city distances or costs
%     SALESMEN (scalar integer) is the number of salesmen to visit the cities
%     MIN_TOUR (scalar integer) is the minimum tour length for any of the
%         salesmen, NOT including the start point
%     POP_SIZE (scalar integer) is the size of the population (should be divisible by 8)
%     NUM_ITER (scalar integer) is the number of desired iterations for the algorithm to run
%     SHOW_PROG (scalar logical) shows the GA progress if true
%     SHOW_RES (scalar logical) shows the GA results if true
% DMAT (float)是一个城市到城市距离或成本的NxN矩阵
% SALESMEN (scalar integer)有多少售货员要访问城市
% MIN_TOUR (scalar integer)每个售货员的最短行程是多少，不包括起点
% POP_SIZE(标量整数)是种群的大小(应该能被8整除)
% NUM_ITER(标量整数)是运行算法所需的迭代次数
% 如果为真，SHOW_PROG(标量逻辑)显示GA的进程
% 如果为真，SHOW_RES(标量逻辑)显示GA结果

% Output:
%     OPT_RTE (integer array) is the best route found by the algorithm
%     OPT_BRK (integer array) is the list of route break points (these specify the indices
%         into the route used to obtain the individual salesman routes)
%     MIN_DIST (scalar float) is the total distance traveled by the salesmen
% OPT_RTE(整数数组)是该算法找到的最佳路由
% OPT_BRK(整数数组)是路由断点的列表(这些指定了索引)转换为用于获取单个销售路线的路线)
% MIN_DIST(标量浮点数)是销售人员走过的总距离
% Process Inputs and Initialize Defaults
nargs = 8;
for k = nargin:nargs-1
    switch k
        case 0
            xy = 10*rand(40,2);
        case 1
            N = size(xy,1);
            a = meshgrid(1:N);
            dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),N,N);
        case 2
            salesmen = 5;
        case 3
            min_tour = 2;
        case 4
            pop_size = 80;
        case 5
            num_iter = 5e3;
        case 6
            show_prog = 1;
        case 7
            show_res = 1;
        otherwise
    end
end
% Verify Inputs
  [~,N] = size(dmat);
  n = N - 1; % Separate Start City

% Sanity Checks
salesmen = max(1,min(n,round(real(salesmen(1)))));
min_tour = max(1,min(floor(n/salesmen),round(real(min_tour(1)))));
pop_size = max(8,8*ceil(pop_size(1)/8));
num_iter = max(1,round(real(num_iter(1))));
show_prog = logical(show_prog(1));
show_res = logical(show_res(1));

% Initializations for Route Break Point Selection
num_brks = salesmen-1;
dof = n - min_tour*salesmen;          % degrees of freedom
addto = ones(1,dof+1);
for k = 2:num_brks
    addto = cumsum(addto);
end
cum_prob = cumsum(addto)/sum(addto);

% Initialize the Populations
pop_rte = zeros(pop_size,n);          % population of routes
pop_brk = zeros(pop_size,num_brks);   % population of breaks
for k = 1:pop_size
    pop_rte(k,:) = randperm(n)+1;
    pop_brk(k,:) = randbreaks();
end

% Run the GA
global_min = Inf;
total_dist = zeros(1,pop_size);
dist_history = zeros(1,num_iter);
tmp_pop_rte = zeros(8,n);
tmp_pop_brk = zeros(8,num_brks);
new_pop_rte = zeros(pop_size,n);
new_pop_brk = zeros(pop_size,num_brks);


for iter = 1:num_iter
    % Evaluate Members of the Population
    for p = 1:pop_size
        d = 0;
        p_rte = pop_rte(p,:);
        p_brk = pop_brk(p,:);
        rng = [[1 p_brk+1];[p_brk n]]';
        for s = 1:salesmen
            d = d + dmat(1,p_rte(rng(s,1))); % Add Start Distance
            for k = rng(s,1):rng(s,2)-1
                d = d + dmat(p_rte(k),p_rte(k+1));
            end
        end
        total_dist(p) = d;
    end

    % Find the Best Route in the Population
    [min_dist,index] = min(total_dist);
    dist_history(iter) = min_dist;
    if min_dist < global_min
        global_min = min_dist;
        opt_rte = pop_rte(index,:);
        opt_brk = pop_brk(index,:);
        
        if show_prog
            % Plot the Best Route
            for s = 1:salesmen
                sprintf('Total Distance = %1.4f, Iteration = %d',min_dist,iter);
            end
        end
    end
    
    % Genetic Algorithm Operators
    rand_grouping = randperm(pop_size);
    for p = 8:8:pop_size
        rtes = pop_rte(rand_grouping(p-7:p),:);
        brks = pop_brk(rand_grouping(p-7:p),:);
        dists = total_dist(rand_grouping(p-7:p));
        [~,idx] = min(dists);
        best_of_8_rte = rtes(idx,:);
        best_of_8_brk = brks(idx,:);
        rte_ins_pts = sort(ceil(n*rand(1,2)));
        I = rte_ins_pts(1);
        J = rte_ins_pts(2);
        for k = 1:8 % Generate New Solutions
            tmp_pop_rte(k,:) = best_of_8_rte;
            tmp_pop_brk(k,:) = best_of_8_brk;
            switch k
                case 2 % Flip
                    tmp_pop_rte(k,I:J) = fliplr(tmp_pop_rte(k,I:J));
                case 3 % Swap
                    tmp_pop_rte(k,[I J]) = tmp_pop_rte(k,[J I]);
                case 4 % Slide
                    tmp_pop_rte(k,I:J) = tmp_pop_rte(k,[I+1:J I]);
                case 5 % Modify Breaks
                    tmp_pop_brk(k,:) = randbreaks();
                case 6 % Flip, Modify Breaks
                    tmp_pop_rte(k,I:J) = fliplr(tmp_pop_rte(k,I:J));
                    tmp_pop_brk(k,:) = randbreaks();
                case 7 % Swap, Modify Breaks
                    tmp_pop_rte(k,[I J]) = tmp_pop_rte(k,[J I]);
                    tmp_pop_brk(k,:) = randbreaks();
                case 8 % Slide, Modify Breaks
                    tmp_pop_rte(k,I:J) = tmp_pop_rte(k,[I+1:J I]);
                    tmp_pop_brk(k,:) = randbreaks();
                otherwise % Do Nothing
            end
        end
        new_pop_rte(p-7:p,:) = tmp_pop_rte;
        new_pop_brk(p-7:p,:) = tmp_pop_brk;
    end
    pop_rte = new_pop_rte;
    pop_brk = new_pop_brk;
end

% if show_res
% % Plots
%     figure('Name','MTSPOFS_GA | Results','Numbertitle','off');
% 
%     subplot(1,2,1);
%     imagesc(dmat([1 opt_rte],[1 opt_rte]));
%     title('Distance Matrix');
%     subplot(1,2,2);
%     plot(dist_history,'b','LineWidth',2);
%     title('Best Solution History');
%     set(gca,'XLim',[0 num_iter+1],'YLim',[0 1.1*max([1 dist_history])]);
% end

% Return Outputs
if nargout
    varargout{1} = opt_rte;
    varargout{2} = opt_brk;
    varargout{3} = min_dist;
end

    % Generate Random Set of Break Points
    function breaks = randbreaks()
        if min_tour == 1 % No Constraints on Breaks
            tmp_brks = randperm(n-1);
            breaks = sort(tmp_brks(1:num_brks));
        else % Force Breaks to be at Least the Minimum Tour Length
            num_adjust = find(rand < cum_prob,1)-1;
            spaces = ceil(num_brks*rand(1,num_adjust));
            adjust = zeros(1,num_brks);
            for kk = 1:num_brks
                adjust(kk) = sum(spaces == kk);
            end
            breaks = min_tour*(1:num_brks) + cumsum(adjust);
        end
    end
end
