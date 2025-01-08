function [] = parallelVTNS(file_path,ifc)
% clc
% clear

data_file = file_path;
p = 0.1;
[customer_n, depot_n, car_n, capicity,VRP4E_DATA] = data_trans2VRP4E(data_file,p);
is_high_risk = zeros(customer_n,1);
% 将数据写入新的文件夹
new_dis_name = "/VRP4E/VRP4E";
new_file_name = strcat("VRP4E",data_file);
new_file_name = strcat(new_file_name,"highp");
new_file_name = strcat(new_file_name,num2str(p*10));

%for shenzhen
% [request,~,~,~,~] = read_shenzhen(file_path);
% new_file_name = file_path;
% is_high_risk = zeros(50,1);
% rand_high = randi(50,[1,5]);
% is_high_risk(rand_high) = 1;

%for CVRP
% new_file_name = file_path;
% [request,~,~,~,~] = read_cvrp_dataset(new_file_name);
% v_num = length(request);
% is_high_risk = zeros(v_num,1);

% fileID = fopen(new_file_name,'w+');%'w’表示打开或创建一个新txt文件以输出
% for i = 1:customer_n+depot_n+1
%     fprintf(fileID,'%d %d %d %d\n',VRP4E_DATA(i,:));
% end
% fclose('all');
global v_Pos;
global d_Pos;
global Dis;
global Dis_depot;
% global request;
% data_file = "p01";
% [request,v_Pos,d_Pos,~,~] = Pos_init(new_file_name);
%距离矩阵
[Dis,Dis_depot] = Dis_ini(v_Pos,d_Pos);

% 在使用CVRP数据集时注释掉这行
% is_high_risk = VRP4E_DATA(2:customer_n+1,4);
is_fuzzy_cluster = [1 ,0, 3, 4]; %1 fcm  2 no  3 cscm/scscm  4 kmeans  
% ifc = 2;

cost_selection = [1 0];%1是discost;0是preventioncost
dis_cost = 0;
prevention_cost = 0;
N = 4;%设置循环轮数
ini_sol_fit = zeros(N,2);%比较使用fcm前后的初始解的质量
cost_swarm = cell(N,1);%分别存储所有的解
% cost = zeros(N,1);

pareto_set_1 =  cell(N,1);
pareto_set_2 =  cell(N,1);
co_pareto_set = cell(N,1);
run_mode = [1,0];
%运行模式：
%1：并行单目标优化
%2：融合后加权多目标优化
bs_1 = [];
bs_2 = [];
%%
t0 = cputime;%计算用时
tic
parfor i=1:N
    %并行两个种群进行两个方向的优化
     size(cost_selection);
     size(is_fuzzy_cluster);
     if(cost_selection(mod(i,2)+1))
         clc_fit_1st=@clc_fit;
         clc_fit_2nd=@(p1,p2,p3)clc_prevention_cost(p1,p2,p3,is_high_risk,10);
         [temp_sol,temp_cost,temp_pareto_set,~,bs] = VTNS_parallel(new_file_name,file_path,is_fuzzy_cluster(ifc),[],clc_fit_1st,clc_fit_2nd,is_high_risk,1);
         pareto_set_1{i} = temp_pareto_set;
         bs_1 = bs;
         cost_swarm{i} = temp_sol;
         cost(i) = temp_cost;
     else
         clc_fit_1st=@(p1,p2,p3)clc_prevention_cost(p1,p2,p3,is_high_risk,10);
         clc_fit_2nd=@clc_fit;
         [temp_sol,temp_cost,temp_pareto_set,~,bs] = VTNS_parallel(new_file_name,file_path,is_fuzzy_cluster(ifc),[],clc_fit_1st,clc_fit_2nd,is_high_risk,1);
         pareto_set_2{i} = temp_pareto_set;
         bs_2 = bs;
         cost_swarm{i} = temp_sol;
         cost(i) = temp_cost;
     end
 end
 toc
fprintf("Two direction optimizing done!\n");
delete(gcp('nocreate'))%仅测试时用，正常运行需要注释掉这行
%求解pareto解集
% pareto_file = strcat(new_file_name,"paretoset.xls");
% if(exist(pareto_file,'file'))
%     pareto_set = self_xlsread(pareto_file);
%     delete(pareto_file);
% else
%     pareto_set = [];
% end
pareto_set = [];
for i=1:N
    if(mod(i,2)==0)
        for j = 1:size(pareto_set_1{i},1)
            pareto_set = pareto_update(pareto_set,pareto_set_1{i}(j,1),pareto_set_1{i}(j,2));
        end
    else
        for j = 1:size(pareto_set_2{i},1)
            pareto_set = pareto_update(pareto_set,pareto_set_2{i}(j,2),pareto_set_2{i}(j,1));
        end
    end
end
% %写入解集
% % self_xlswrite(pareto_file,pareto_set);
% % fprintf("Co-optimizing starting!\n");
% %%
% %再进行单目标优化
% %建立第三个种群，用以探索pareto中心部分的解
% delete(gcp('nocreate'))%仅测试时用，正常运行需要注释掉这行
% 
% method = ["fcm", "no", "cscm", "kmeans"];
% save_mat = strcat(file_path, '_', method(ifc), '.mat');
% save_pareto_set = unique(pareto_set,"rows");
% save(save_mat, "save_pareto_set");
% fprintf(method(ifc) + "save done!")
tic
co_cost_swarm = cell(N,1);
parfor i = 1:N
    size(cost_selection);
    size(is_fuzzy_cluster);
    
    if(cost_selection(mod(i,2)+1))
        clc_fit_2nd=@clc_fit;
        clc_fit_1st=@(p1,p2,p3)clc_prevention_cost(p1,p2,p3,is_high_risk,10);
        [~,~,temp_pareto_set] = VTNS_parallel(new_file_name,file_path,is_fuzzy_cluster(ifc),bs_1,clc_fit_1st,clc_fit_2nd,is_high_risk,1);
        t_1 = temp_pareto_set(:,1);
        t_2 = temp_pareto_set(:,2);
        co_pareto_set{i} = [t_2 t_1];
    else
        clc_fit_1st=@clc_fit;
        clc_fit_2nd=@(p1,p2,p3)clc_prevention_cost(p1,p2,p3,is_high_risk,10);
        [~,~,temp_pareto_set] = VTNS_parallel(new_file_name,file_path,is_fuzzy_cluster(ifc),bs_2,clc_fit_1st,clc_fit_2nd,is_high_risk,1);
        co_pareto_set{i} = temp_pareto_set;
    end
end
% toc
fprintf("Co-optimizing done!\n");
delete(gcp('nocreate'))

% 使用rank求最优pareto 前沿
optima_pareto = [];
for nn = 1:N
    tt = unique(co_pareto_set{nn},'rows');
    optima_pareto = [optima_pareto;tt];
end
optima_pareto = unique(optima_pareto,"rows");
% 需要计算结果
pareto_rank_end = FastNonDominatedSorting_Vectorized(optima_pareto);
index_1 = find(pareto_rank_end==1);
optima_pareto_set = optima_pareto(index_1,:);
new_ff = strcat('./new_fcm/',file_path);
new_ff = strcat(new_ff,'_r.mat');
save(new_ff,"optima_pareto_set");

% 只计算用时
% total_time = toc;
% new_ff_t = strcat('./low_FCM/re_time_mode',file_path);
% new_ff_t = strcat(new_ff_t,num2str(ifc));
% new_ff_t = strcat(new_ff_t,'9999.mat');
% save(new_ff_t,'total_time');
%求解联合优化的pareto解集
% co_pareto_file = strcat(new_file_name,"_coparetoset.xls");
% if(exist(co_pareto_file,'file'))
%     pareto_end = self_xlsread(co_pareto_file);
%     delete(co_pareto_file);
% else
%     pareto_end = [];
% end
% for i=1:N
%     for j = 1:length(co_pareto_set{i})
%         pareto_end = pareto_update(pareto_end,co_pareto_set{i}(j,1),co_pareto_set{i}(j,2));
%     end
% end 
% self_xlswrite(co_pareto_file,pareto_end);%
% %end
end
