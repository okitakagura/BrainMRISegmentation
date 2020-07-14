close all;
clear all;
clc;

nii = load_nii( 'sub_TJ001_brain_FLIRT.nii.gz' );  % 装载.nii数据
img = nii.img;  % 因为这个文件有img和head二个部分，其中img部分是图像数据
save image.mat img  % 将数据变成mat格式
load 'image.mat'  % 加载数据
[n1, n2, n3] = size(img);   % 获取.nii文件的三个维度，一般1、2维是图像维度，第三维是切片
v_p1 = 0; 
v_p2 = 0; 
v_p3 = 0; 
v_p4 = 0;
n_start = 1;
n_end =79;
for slicenum = n_start:n_end
data_load = imrotate(img(:,:,slicenum),90);
imshow(data_load,[]);  % 这个是正常显示第45个切片的图像

%10??154切片
% mark = Mark('brainweb/phantom_1.0mm_normal_crisp.rawb',100);
% read = readrawb('brainweb/t1_icbm_normal_1mm_pn7_rf20.rawb',100);
% [n1,n2] = size(read);
% [row,col] = size(read);
% read_new = read;
% for i = 1:row
%     for j = 1:col
%         if mark(i,j) == 0
%             read_new(i,j)=0;
%         end
%     end
% end
% read_new=imrotate(read_new, 90); 
% real_label=imrotate(mark, 90);
% real_count1 = 0;
% real_count2 = 0;
% real_count3 = 0;
% real_count4 = 0;
% for x=1:n2
%     for y=1:n1
%         if real_label(x,y) == 0
%             real_count1 = real_count1 + 1;
%         elseif real_label(x,y) == 1
%             real_count2 = real_count2 + 1;
%         elseif real_label(x,y) == 2
%             real_count3 = real_count3 + 1;
%         elseif real_label(x,y) == 3
%             real_count4 = real_count4 + 1;
%         end
%     end
% end
% data_load = read_new;
% imshow(data_load,[]);


data2 = data_load;
[a b]=size(data2);%图像矩阵的大小
data=data2(:);%转化为列向量
data=double(data);
data_n = size(data, 1); % 求出data的第一维(rows)数,即样本个数
in_n = size(data, 2);   % 求出data的第二维(columns)数，即特征值长度

options=[2;150;1e-5;1];%指定隶属度矩阵的的模糊指数，算法迭代次数，迭代终止条件等
cluster_n=4; %指定类别数
expo = options(1);          % 隶属度矩阵U的指数
max_iter = options(2);		% 最大迭代次数 
min_impro = options(3);		% 隶属度最小变化量,迭代终止条件
display = options(4);		% 每次迭代是否输出信息标志 
alfa=1.5;

% [U, center,obj_fcn] = FCMClust2(data2, cluster_n, expo, max_iter, min_impro, display);
%[U, center,obj_fcn] = FCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
%[U, center,obj_fcn] = KFCMClust2brainweb(data_load, cluster_n, expo, max_iter, min_impro, display);
%[U, center,obj_fcn] = KFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
[U, center,obj_fcn] = MyKFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);

fm = U.^expo;
V_pc = sum(fm(:))/data_n;
V_xieben = XieBeniInverted(U, expo, center, data);

%图像分割
data=data';
wholeG=zeros(size(data));
maxU=max(U);
x=[1,2,3,4];
[center,id] = sort(center);
x=x(id);
for k=1:cluster_n
    indexk=(U(k,:)==maxU);
    count{k} = sum(indexk(:));
    Ik = indexk.*data;
    Ik = reshape(Ik,a,b);
    result{k} = Ik;%result{k}记录第Ｋ类中的像素及其位置信息
    wholeG(indexk) = x(k);%wholeG记录分割后整体图像的信息
end

wholeG=reshape(wholeG,a,b);%将向量转化为a*b大小的矩阵
% [label_new,accuracy]=succeed(real_label,cluster_n,wholeG);
% wholeG = reshape(label_new,a,b);
%显示分割结果
subplot(1,2,1);imshow(data_load,[]);title('原始图像');
subplot(1,2,2);imshow(wholeG,[]);title('分割后的结果');
% 
% for i=1:cluster_n
%    figure;
%    imshow(result{x(i)});
% end
%体积计算
s_p2 = count{x(2)} * 0.04;
s_p3 = count{x(3)} * 0.04;
s_p4 = count{x(4)} * 0.04;
if slicenum == n_start || slicenum == n_end
    v_p2 = v_p2+s_p2;
    v_p3 = v_p3+s_p3;
    v_p4 = v_p4+s_p4;
else
    v_p2 = v_p2+2*s_p2;
    v_p3 = v_p3+2*s_p3;
    v_p4 = v_p4+2*s_p4;
end
end
V_p2 = v_p2*0.2/2;
V_p3 = v_p3*0.2/2;
V_p4 = v_p4*0.2/2;
V_p2
V_p3
V_p4


function xieBeni = XieBeniInverted(U, uExpoent, center, data)
%XIEBENIINVERTED Implementation of the measure of cluster validation of the Xie-Beni.

% Calculate Xie-Beni Inverted cluster validation method

dist = sum(sum((U.^uExpoent) .* (distfcm(center, data).^2)));
minimum = intmax;
for i = 1 : size(center, 1)-1
    minimumAux = min(distfcmfp(center(i, :), center(i+1:end, :)).^2);
    if(minimum > minimumAux)
        minimum = minimumAux;
    end
end
xieBeni = dist / (size(data, 1) .* minimum);
end

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
end

function out = distfcmfp(center, data)
%DISTFCMFP Distance measure in fuzzy c-mean clustering with focal point.
%	OUT = DISTFCMFP(CENTER, DATA) calculates the Euclidean distance
%	between each row in CENTER and each row in DATA, and returns a
%	distance matrix OUT of size M by N, where M and N are row
%	dimensions of CENTER and DATA, respectively, and OUT(I, J) is
%	the distance between CENTER(I,:) and DATA(J,:).
out = zeros(size(center, 1), size(data, 1));
% fill the output matrix
if size(center, 2) > 1,
    for k = 1:size(center, 1),
        out(k, :) = sqrt(sum(((data-ones(size(data, 1), 1)*center(k, :)).^2)'));
    end
else	% 1-D data
    for k = 1:size(center, 1),
        out(k, :) = abs(center(k)-data)';
    end
end
end


