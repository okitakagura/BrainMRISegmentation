function [accuracy, center, U, obj_fcn] = KGFCM_S(data, data1,cluster_n,alfa,beta,options)
% 输入：
%   data        ---- nxm矩阵,表示n个样本,每个样本具有m的维特征值
%   data1        ---- nxm矩阵,表示均值或者中值，均值为FCM_S1算法，中值为FCM_S2算法
%   N_cluster   ---- 标量,表示聚合中心数目,即类别数
%   alfa        ---- punishment factor
%   beta        ---- control factor
%   options     ---- 4x1矩阵，其中
%       options(1):  隶属度矩阵U的指数，>1                  (缺省值: 2.0)
%       options(2):  最大迭代次数                           (缺省值: 100)
%       options(3):  隶属度最小变化量,迭代终止条件           (缺省值: 1e-5)
%       options(4):  每次迭代是否输出信息标志                (缺省值: 1)
% 输出：
%   center      ---- 聚类中心
%   U           ---- 隶属度矩阵
%   obj_fcn     ---- 目标函数值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
close all;
clear all;
clc;

tic

% nii = load_nii( 'sub_TJ001_brain_FLIRT.nii.gz' );  % 装载.nii数据
% img = nii.img;  % 因为这个文件有img和head二个部分，其中img部分是图像数据
% save image.mat img  % 将数据变成mat格式
% load 'image.mat'  % 加载数据
% [n1, n2, n3] = size(img);   % 获取.nii文件的三个维度，一般1、2维是图像维度，第三维是切片
% data_load = imrotate(img(:,:,47),90);
% imshow(data_load,[]);  % 这个是正常显示第45个切片的图像


mark = Mark('brainweb/phantom_1.0mm_normal_crisp.rawb',80);
read = readrawb('brainweb/t1_icbm_normal_1mm_pn5_rf20.rawb',80);
[n1,n2] = size(read);
[row,col] = size(read);
read_new = read;
for i = 1:row
    for j = 1:col
        if mark(i,j) == 0
            read_new(i,j)=0;
        end
    end
end
read_new=imrotate(read_new, 90); 
real_label=imrotate(mark, 90);
real_count1 = 0;
real_count2 = 0;
real_count3 = 0;
real_count4 = 0;
for x=1:n2
    for y=1:n1
        if real_label(x,y) == 0
            real_count1 = real_count1 + 1;
        elseif real_label(x,y) == 1
            real_count2 = real_count2 + 1;
        elseif real_label(x,y) == 2
            real_count3 = real_count3 + 1;
        elseif real_label(x,y) == 3
            real_count4 = real_count4 + 1;
        end
    end
end
data_load = read_new;
imshow(data_load,[]);  % 这个是正常显示第45个切片的图像




data2=data_load;%读入原始图像
options=[2;150;1e-5;1];%指定隶属度矩阵的的模糊指数，算法迭代次数，迭代终止条件等
cluster_n=4; %指定类别数
alfa=0.9;
beta=3;
[a b]=size(data2);%图像矩阵的大小

img1 = data2;
width = 3;  %局部窗尺寸
delta = (width-1)/2;
for i = delta+1:a-delta
   for j = delta+1:b-delta
%         temp = data2(i-delta:i+delta,j-delta:j+delta);
%         temp = sort(temp(:));  
%         img1(i,j) = temp((length(temp)+1)/2);
       img1(i,j) = sum(sum(data2(i-delta:i+delta,j-delta:j+delta)))/(width*width);
   end
end
data1=img1(:);%转化为列向量
data1=double(data1);
data=data2(:);%转化为列向量
data=double(data);
data_n = size(data, 1); % 求出data的第一维(rows)数,即样本个数
in_n = size(data, 2);   % 求出data的第二维(columns)数，即特征值长度
%将options 中的分量分别赋值给四个变量;
expo = options(1);          % 隶属度矩阵U的指数
max_iter = options(2);      % 最大迭代次数 
min_impro = options(3);     % 隶属度最小变化量,迭代终止条件
display = options(4);       % 每次迭代是否输出信息标志 
obj_fcn = zeros(max_iter, 1);   % 初始化输出参数obj_fcn     
U = initfcm(cluster_n, data_n);     % 初始化模糊分配矩阵,使U满足列上相加为1, 
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
center = mf*(data+beta*data1)./((ones(size(data, 2), 1)*sum(mf'))'*(1+beta)); 
segma2=estimateSegma(data);
% Main loop  主要循环
for i = 1:max_iter,
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
%     U0 = U;
    [U, center, obj_fcn(i)] = stepfcm(data,data1, center,U, cluster_n, expo,alfa,beta,segma2);
    if display, 
        fprintf('FCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
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
     

%图像分割
data=data';
wholeG=zeros(size(data));
maxU=max(U);
for k=1:cluster_n
    indexk=(U(k,:)==maxU);
    count{k} = sum(indexk(:));
    Ik = indexk.*data;
    Ik = reshape(Ik,a,b);
    result{k} = Ik;%result{k}记录第Ｋ类中的像素及其位置信息
    wholeG(indexk) = k;%wholeG记录分割后整体图像的信息
end
wholeG=reshape(wholeG,a,b);%将向量转化为a*b大小的矩阵
[label_new,accuracy]=succeed(real_label,cluster_n,wholeG);
%显示分割结果
subplot(1,2,1);imshow(data_load,[]);title('原始图像');
subplot(1,2,2);imshow(wholeG,[]);title('kfcm_s方法分割后的结果');
for i=1:cluster_n
   figure;
   imshow(result{i});
end
count{1}
count{2}
count{3}
count{4}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% 子函数2
function [U_new, center, obj_fcn] = stepfcm(data,data1,center, U, cluster_n, expo,alfa,beta,segma2)
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
K = Kernel(data, center,segma2);        
K1 = Kernel(data1, center,segma2); 
center= sum(mf.*(K.*(ones(cluster_n,1)*data')+beta*K1.*(ones(cluster_n,1)*data1')),2)./sum(mf.*(K+beta*K1),2); % 新聚类中心
dist = 2*(1-K);
dist1 = 2*(1-K1);
alfa_j = min(dist.^2)*alfa;    %alfa_j  is control factor vector
obj_fcn = sum(sum((dist.^2).*mf))+ sum(alfa_j.*sum(U.*(1-U.^(expo-1))))+beta*sum(sum((dist1.^2).*mf));  % 计算目标函数值 
tmp = (dist.^2-ones(cluster_n,1)*alfa_j+beta*dist1.^2).^(-1/(expo-1));
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % 计算新的隶属度矩阵 (5.3)式
% tmp = dist.^(2/(expo-1));
% U_new =ones(cluster_n, size(data,1))./(tmp.*(ones(cluster_n,1)*sum(ones(cluster_n, size(data,1))./tmp)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 子函数3
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
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end

function segma2=estimateSegma(data)
%estimate the parameter σ in Kernel distance
% Input：
%   data           ----samples
%
% Output：
%   segma2       ----the parameter σ^2
mean_data = mean(data);
segma2 = (data-mean_data)'*(data-mean_data)/size(data,1);

function d=Kernel(data,center,segma2)
%compute Kernel distance of vector p and v
% Input：
%   p,v           ----vector p and v
%   segma2         ----the parameter σ^2
% Output：
%   d             ----Kernel distance
d = zeros(size(center, 1), size(data, 1));
data = data';
for k = 1:size(center, 1), % 对每一个聚类中心
    % 每一次循环求得所有样本点到一个聚类中心的距离
    d(k, :) = exp(-(data-center(k)).^2/segma2);
end


