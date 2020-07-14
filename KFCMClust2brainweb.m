function[U, center,obj_fcn] = KFCMClust2brainweb(data_load, cluster_n, expo, max_iter, min_impro, display)
% KFCMClust2brainweb.m   采用结合核函数的模糊C均值对数据集data聚为cluster_n类   
% 输入：
%   data_load   ---- 数据集
%   cluster_n     ---- 聚类数目
%   expo  ---- 隶属度矩阵U的指数，>1                  (缺省值: 2.0)
%   max_iter  ---- 最大迭代次数                           (缺省值: 100)
%   min_impro  ---- 隶属度最小变化量,迭代终止条件           (缺省值: 1e-5)
%   display  ---- 每次迭代是否输出信息标志                (缺省值: 1)
% 输出：
%   center      ---- 聚类中心
%   U           ---- 隶属度矩阵
%   obj_fcn     ---- 目标函数值


data1=data_load;
[a b]=size(data1);
data=data1(:);
data=double(data);
data_n = size(data, 1); % 求出data的第一维(rows)数,即样本个数
in_n = size(data, 2);   % 求出data的第二维(columns)数，即特征值长度

kernel_b=estimateSegma(data);
obj_fcn = zeros(max_iter, 1);	% 初始化输出参数obj_fcn
U = initkfcm(cluster_n, data_n);     % 初始化模糊分配矩阵,使U满足列上相加为1

% 初始化聚类中心：从样本数据点中任意选取cluster_n个样本作为聚类中心。当然，
% 如果采用某些先验知识选取中心或许能够达到加快稳定的效果，但目前不具备这功能
index = randperm(data_n);   % 对样本序数随机排列
center_old = data(index(1:cluster_n),:);  % 选取随机排列的序数的前cluster_n个

% Main loop  主要循环
for i = 1:max_iter,
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
	[U, center, obj_fcn(i)] = stepkfcm(data,U,center_old, expo, kernel_b);
	if display, 
		fprintf('KFCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    center_old = center;    % 用新的聚类中心代替老的聚类中心
	% 终止条件判别
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% 实际迭代次数 
obj_fcn(iter_n+1:max_iter) = [];



%子函数
function U = initkfcm(cluster_n, data_n)
% 初始化fcm的隶属度函数矩阵
% 输入:
%   cluster_n   ---- 聚类中心个数
%   data_n      ---- 样本点数
% 输出：
%   U           ---- 初始化的隶属度矩阵
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);

% 子函数
function [U_new,center_new,obj_fcn] = stepkfcm(data,U,center,expo,kernel_b)
% 模糊C均值聚类时迭代的一步
% 输入：
%   data        ---- nxm矩阵,表示n个样本,每个样本具有m的维特征值
%   U           ---- 隶属度矩阵
%   center      ---- 聚类中心
%   expo        ---- 隶属度矩阵U的指数         
%   kernel_b    ---- 高斯核函数的参数
% 输出：
%   U_new       ---- 迭代计算出的新的隶属度矩阵
%   center_new  ---- 迭代计算出的新的聚类中心
%   obj_fcn     ---- 目标函数值
feature_n = size(data,2);  % 特征维数
cluster_n = size(center,1); % 聚类个数
mf = U.^expo;       % 隶属度矩阵进行指数运算（c行n列)

% 计算新的聚类中心;
KernelMat = gaussKernel(center,data,kernel_b); % 计算高斯核矩阵(c行n列)
num = mf.*KernelMat * data;   
den = sum(mf.*KernelMat,2); 
center_new = num./(den*ones(1,feature_n)); % 计算新的聚类中心(c行p列,c个中心)

% 计算新的隶属度矩阵
kdist = distkfcm(center_new, data, kernel_b);    % 计算距离矩阵
obj_fcn = sum(sum((kdist.^2).*mf));  % 计算目标函数值
tmp = kdist.^(-1/(expo-1));     
U_new = tmp./(ones(cluster_n, 1)*sum(tmp)); 

% 子函数
function out = distkfcm(center, data, kernel_b)
% 计算样本点距离聚类中心的距离
% 输入：
%   center     ---- 聚类中心
%   data       ---- 样本点
% 输出：
%   out        ---- 距离
cluster_n = size(center, 1);
data_n = size(data, 1);
out = zeros(cluster_n, data_n);
for i = 1:cluster_n % 对每个聚类中心 
    vi = center(i,:);
    out(i,:) = 2-2*gaussKernel(vi,data,kernel_b);
end

function segma2=estimateSegma(data)
% 计算高斯距离的参数
% 输入：
%   data           ---- 数据点
%
% 输出：
%   segma2       ---- 参数σ^2
mean_data = mean(data);
segma2 = (data-mean_data)'*(data-mean_data)/size(data,1);

% 子函数
function out = gaussKernel(center,data,kernel_b)
% 高斯核函数计算
% 输入:
%   center      ---- 模糊聚类中心
%   data        ---- 样本数据点
%   kernel_b    ---- 高斯核参数
% 输出：
%   out         ---- 高斯核计算结果

dist = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1), % 对每一个聚类中心
    % 每一次循环求得所有样本点到一个聚类中心的距离
    dist(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end
out = exp(-dist.^2/kernel_b);