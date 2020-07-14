function [U, center,obj_fcn] = MyKFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa)
% MyKFCM_S.m   采用改进的结合邻域信息和核函数的模糊C均值对数据集data聚为cluster_n类   
% 输入：
%   data_load   ---- 数据集
%   cluster_n     ---- 聚类数目
%   expo  ---- 隶属度矩阵U的指数，>1                  (缺省值: 2.0)
%   max_iter  ---- 最大迭代次数                           (缺省值: 100)
%   min_impro  ---- 隶属度最小变化量,迭代终止条件           (缺省值: 1e-5)
%   display  ---- 每次迭代是否输出信息标志                (缺省值: 1)
%   alfa  ---- 邻域信息权重
% 输出：
%   center      ---- 聚类中心
%   U           ---- 隶属度矩阵
%   obj_fcn     ---- 目标函数值
     

data2=data_load;%读入原始图像
[a b]=size(data2);%图像矩阵的大小

img1 = data2;
width = 3;  %局部窗尺寸
delta = (width-1)/2;
for i = delta+1:a-delta
   for j = delta+1:b-delta
%        piexnum = 0;
%        piexvalue = 0;
%        avg_p = sum(sum(data2(i-delta:i+delta,j-delta:j+delta)))/(width*width);
%        std_p = std2(data2(i-delta:i+delta,j-delta:j+delta));
%        avg_gap = sum(sum(abs(data2(i-delta:i+delta,j-delta:j+delta)-ones(width)*data2(i,j))))/(width*width);
%        for ii=i-delta:i+delta
%            for jj=j-delta:j+delta
%                left = abs(data2(ii,jj) - avg_p)+abs(data2(ii,jj)-data2(i,j)); 
%                right = std_p + avg_gap;
%                if left <= right
%                    piexvalue = piexvalue+data2(ii,jj);
%                    piexnum = piexnum + 1;
%                else
%                end
%            end
%        end
%        img1(i,j) = piexvalue/piexnum;
         temp = data2(i-delta:i+delta,j-delta:j+delta);
         temp = sort(temp(:));  
         img1(i,j) = temp((length(temp)+1)/2);
         
   end
end
data1=img1(:);%转化为列向量
data1=double(data1);
data=data2(:);%转化为列向量
data=double(data);
data_n = size(data, 1); % 求出data的第一维(rows)数,即样本个数
in_n = size(data, 2);   % 求出data的第二维(columns)数，即特征值长度


obj_fcn = zeros(max_iter, 1);   % 初始化输出参数obj_fcn     
U = initfcm(cluster_n, data_n);     % 初始化模糊分配矩阵,使U满足列上相加为1, 
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
%center = mf*(data+alfa*data1)./((ones(size(data, 2), 1)*sum(mf'))'*(1+alfa)); %初始聚类中心

%center = [0.9769;108.4991;212.6035;302.3221];

center=kmeans4(data_load);
segma=estimateSegma(data);

% %初始化聚类中心：从样本数据点中任意选取cluster_n个样本作为聚类中心。当然，
% % 如果采用某些先验知识选取中心或许能够达到加快稳定的效果，但目前不具备这功能
% index = randperm(data_n);   % 对样本序数随机排列
% center_old = data(index(1:cluster_n),:);  % 选取随机排列的序数的前cluster_n个

% Main loop  主要循环
for i = 1:max_iter,
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
%     U0 = U;
    [U, center obj_fcn(i)] = stepfcm(data,data1,center,U, cluster_n, expo,alfa,segma);
    if display, 
        fprintf('MyKFCM_S:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    % 终止条件判别
    if i > 1,
        if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, 
            break;
        end,
%         if sum(sum(abs(U-U0)))<1
%             break;
%         end
    end
end

iter_n = i; % 实际迭代次数 
obj_fcn(iter_n+1:max_iter) = [];
     

% 子函数1
function U = initfcm(cluster_n, data_n)
% 初始化fcm的隶属度函数矩阵
% 输入:
%   cluster_n   ---- 聚类中心个数
%   data_n      ---- 样本点数
% 输出：
%   U           ---- 初始化的隶属度矩阵
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);

     
% 子函数2
function [U_new, center, obj_fcn] = stepfcm(data,data1,center,U, cluster_n, expo,alfa,segma)
% 模糊C均值聚类时迭代的一步
% 输入：
%   data        ---- nxm矩阵,表示n个样本,每个样本具有m的维特征值
%   data1        ---- nxm矩阵,表示均值或者中值，均值为FCM_S1算法，中值为FCM_S2算法
%   U           ---- 隶属度矩阵
%   cluster_n   ---- 标量,表示聚合中心数目,即类别数
%   expo        ---- 隶属度矩阵U的指数   
%   alfa        ---- punishment factor
% 输出：
%   U_new       ---- 迭代计算出的新的隶属度矩阵
%   center      ---- 迭代计算出的新的聚类中心
%   obj_fcn     ---- 目标函数值
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
K = Kernel(data, center,segma);        
K1 = Kernel(data1, center,segma); 
%center = mf*(data+alfa*data1)./((ones(size(data, 2), 1)*sum(mf'))'*(1+alfa));
center= sum(mf.*(K.*(ones(cluster_n,1)*data')+alfa*K1.*(ones(cluster_n,1)*data1')),2)./sum(mf.*(K+alfa*K1),2); % 新聚类中心
dist = 2*(1-K);
dist1 = 2*(1-K1);
% dist = distfcm(center, data);       % 计算距离矩阵
% dist1 = distfcm(center,data1);
obj_fcn = sum(sum((dist.^2).*mf))+alfa*sum(sum((dist1.^2).*mf));  % 计算目标函数值 
tmp = (dist.^2+alfa*dist1.^2).^(-1/(expo-1));
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % 计算新的隶属度矩阵

     
function segma2=estimateSegma(data)
% 计算高斯距离的参数
% 输入：
%   data           ---- 数据点
%
% 输出：
%   segma2       ---- 参数σ^2
mean_data = mean(data);
segma2 = (data-mean_data)'*(data-mean_data)/size(data,1);

function d=Kernel(data,center,segma2)
% 计算高斯距离
% 输入：
%   data           ---- 数据点
%   center           ---- 聚类中心
%   segma2         ---- 参数σ
% 输出：
%   d             ---- 高斯距离
d = zeros(size(center, 1), size(data, 1));
data = data';
for k = 1:size(center, 1), % 对每一个聚类中心
    % 每一次循环求得所有样本点到一个聚类中心的距离
    d(k, :) = exp(-(data-center(k)).^2/segma2);
end

