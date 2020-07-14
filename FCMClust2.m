%%%标准的fcm算法
function [U, center,obj_fcn] = FCMClust2(data2, cluster_n, expo, max_iter, min_impro, display)
% FCMClust.m   采用模糊C均值对数据集data聚为cluster_n类   
% 输入：
%   data2   ---- 数据集
%   cluster_n     ---- 聚类数目
%   expo  ---- 隶属度矩阵U的指数，>1                  (缺省值: 2.0)
%   max_iter  ---- 最大迭代次数                           (缺省值: 100)
%   min_impro  ---- 隶属度最小变化量,迭代终止条件           (缺省值: 1e-5)
%   display  ---- 每次迭代是否输出信息标志                (缺省值: 1)
% 输出：
%   center      ---- 聚类中心
%   U           ---- 隶属度矩阵
%   obj_fcn     ---- 目标函数值


[a b]=size(data2);%图像矩阵的大小
data=data2(:);%转化为列向量
data=double(data);
data_n = size(data, 1); % 求出data的第一维(rows)数,即样本个数
in_n = size(data, 2);   % 求出data的第二维(columns)数，即特征值长度
obj_fcn = zeros(max_iter, 1);	% 初始化输出参数obj_fcn
U = initfcm(cluster_n, data_n);     % 初始化模糊分配矩阵,使U满足列上相加为1,

% Main loop  主要循环
for i = 1:max_iter,
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
	[U, center, obj_fcn(i)] = stepfcm(data, U, cluster_n, expo);
	if display
		fprintf('FCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	end
	% 终止条件判别
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, 
            break;
        end,
	end
end

iter_n = i;	% 实际迭代次数 
obj_fcn(iter_n+1:max_iter) = [];


% 子函数:用于初始化隶属度矩阵U
function U = initfcm(cluster_n, data_n)
% 初始化fcm的隶属度函数矩阵
% 输入:
%   cluster_n   ---- 聚类中心个数
%   data_n      ---- 样本点数
% 输出：
%   U           ---- 初始化的隶属度矩阵
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);%初始化时随机产生一个满足每一列和为1的隶属度矩阵


% 子函数：用于计算每次迭代每一次中产生的新的隶属度矩阵u和聚类中心v
function [U_new, center, obj_fcn] = stepfcm(data, U, cluster_n, expo)
% 模糊C均值聚类时迭代的一步
% 输入：
%   data        ---- nxm矩阵,表示n个样本,每个样本具有m的维特征值
%   U           ---- 隶属度矩阵
%   cluster_n   ---- 标量,表示聚合中心数目,即类别数
%   expo        ---- 隶属度矩阵U的指数                      
% 输出：
%   U_new       ---- 迭代计算出的新的隶属度矩阵
%   center      ---- 迭代计算出的新的聚类中心
%   obj_fcn     ---- 目标函数值
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))'); % 新聚类中心
dist = distfcm(center, data);       % 计算距离矩阵
obj_fcn = sum(sum((dist.^2).*mf));  % 计算目标函数值 
tmp = dist.^(-2/(expo-1));     
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % 计算新的隶属度矩阵

% 子函数：计算样本点到每一个聚类中心的欧氏距离
function out = distfcm(center, data)
% 计算样本点距离聚类中心的距离
% 输入：
%   center     ---- 聚类中心
%   data       ---- 样本点
% 输出：
%   out        ---- 距离
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1), % 对每一个聚类中心
    % 每一次循环求得所有样本点到一个聚类中心的距离
   % out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2),2)',1);
     out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end


