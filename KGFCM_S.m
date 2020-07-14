function [accuracy, center, U, obj_fcn] = KGFCM_S(data, data1,cluster_n,alfa,beta,options)
% ���룺
%   data        ---- nxm����,��ʾn������,ÿ����������m��ά����ֵ
%   data1        ---- nxm����,��ʾ��ֵ������ֵ����ֵΪFCM_S1�㷨����ֵΪFCM_S2�㷨
%   N_cluster   ---- ����,��ʾ�ۺ�������Ŀ,�������
%   alfa        ---- punishment factor
%   beta        ---- control factor
%   options     ---- 4x1��������
%       options(1):  �����Ⱦ���U��ָ����>1                  (ȱʡֵ: 2.0)
%       options(2):  ����������                           (ȱʡֵ: 100)
%       options(3):  ��������С�仯��,������ֹ����           (ȱʡֵ: 1e-5)
%       options(4):  ÿ�ε����Ƿ������Ϣ��־                (ȱʡֵ: 1)
% �����
%   center      ---- ��������
%   U           ---- �����Ⱦ���
%   obj_fcn     ---- Ŀ�꺯��ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
close all;
clear all;
clc;

tic

% nii = load_nii( 'sub_TJ001_brain_FLIRT.nii.gz' );  % װ��.nii����
% img = nii.img;  % ��Ϊ����ļ���img��head�������֣�����img������ͼ������
% save image.mat img  % �����ݱ��mat��ʽ
% load 'image.mat'  % ��������
% [n1, n2, n3] = size(img);   % ��ȡ.nii�ļ�������ά�ȣ�һ��1��2ά��ͼ��ά�ȣ�����ά����Ƭ
% data_load = imrotate(img(:,:,47),90);
% imshow(data_load,[]);  % �����������ʾ��45����Ƭ��ͼ��


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
imshow(data_load,[]);  % �����������ʾ��45����Ƭ��ͼ��




data2=data_load;%����ԭʼͼ��
options=[2;150;1e-5;1];%ָ�������Ⱦ���ĵ�ģ��ָ�����㷨����������������ֹ������
cluster_n=4; %ָ�������
alfa=0.9;
beta=3;
[a b]=size(data2);%ͼ�����Ĵ�С

img1 = data2;
width = 3;  %�ֲ����ߴ�
delta = (width-1)/2;
for i = delta+1:a-delta
   for j = delta+1:b-delta
%         temp = data2(i-delta:i+delta,j-delta:j+delta);
%         temp = sort(temp(:));  
%         img1(i,j) = temp((length(temp)+1)/2);
       img1(i,j) = sum(sum(data2(i-delta:i+delta,j-delta:j+delta)))/(width*width);
   end
end
data1=img1(:);%ת��Ϊ������
data1=double(data1);
data=data2(:);%ת��Ϊ������
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����
%��options �еķ����ֱ�ֵ���ĸ�����;
expo = options(1);          % �����Ⱦ���U��ָ��
max_iter = options(2);      % ���������� 
min_impro = options(3);     % ��������С�仯��,������ֹ����
display = options(4);       % ÿ�ε����Ƿ������Ϣ��־ 
obj_fcn = zeros(max_iter, 1);   % ��ʼ���������obj_fcn     
U = initfcm(cluster_n, data_n);     % ��ʼ��ģ���������,ʹU�����������Ϊ1, 
mf = U.^expo;       % �����Ⱦ������ָ��������
center = mf*(data+beta*data1)./((ones(size(data, 2), 1)*sum(mf'))'*(1+beta)); 
segma2=estimateSegma(data);
% Main loop  ��Ҫѭ��
for i = 1:max_iter,
    %�ڵ�k��ѭ���иı��������ceneter,�ͷ��亯��U��������ֵ;
%     U0 = U;
    [U, center, obj_fcn(i)] = stepfcm(data,data1, center,U, cluster_n, expo,alfa,beta,segma2);
    if display, 
        fprintf('FCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    % ��ֹ�����б�
    if i > 1,
        if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, 
            break;
        end,
%         if sum(sum(abs(U-U0)))<1
%             break;
%         end
    end
end
     
iter_n = i; % ʵ�ʵ������� 
obj_fcn(iter_n+1:max_iter) = [];
     

%ͼ��ָ�
data=data';
wholeG=zeros(size(data));
maxU=max(U);
for k=1:cluster_n
    indexk=(U(k,:)==maxU);
    count{k} = sum(indexk(:));
    Ik = indexk.*data;
    Ik = reshape(Ik,a,b);
    result{k} = Ik;%result{k}��¼�ڣ����е����ؼ���λ����Ϣ
    wholeG(indexk) = k;%wholeG��¼�ָ������ͼ�����Ϣ
end
wholeG=reshape(wholeG,a,b);%������ת��Ϊa*b��С�ľ���
[label_new,accuracy]=succeed(real_label,cluster_n,wholeG);
%��ʾ�ָ���
subplot(1,2,1);imshow(data_load,[]);title('ԭʼͼ��');
subplot(1,2,2);imshow(wholeG,[]);title('kfcm_s�����ָ��Ľ��');
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
% �Ӻ���1
function U = initfcm(cluster_n, data_n)
% ��ʼ��fcm�������Ⱥ�������
% ����:
%   cluster_n   ---- �������ĸ���
%   data_n      ---- ��������
% �����
%   U           ---- ��ʼ���������Ⱦ���
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% �Ӻ���2
function [U_new, center, obj_fcn] = stepfcm(data,data1,center, U, cluster_n, expo,alfa,beta,segma2)
% ģ��C��ֵ����ʱ������һ��
% ���룺
%   data        ---- nxm����,��ʾn������,ÿ����������m��ά����ֵ
%   data1        ---- nxm����,��ʾ��ֵ������ֵ����ֵΪFCM_S1�㷨����ֵΪFCM_S2�㷨
%   U           ---- �����Ⱦ���
%   cluster_n   ---- ����,��ʾ�ۺ�������Ŀ,�������
%   expo        ---- �����Ⱦ���U��ָ��   
%   alfa        ---- punishment factor
% �����
%   U_new       ---- ������������µ������Ⱦ���
%   center      ---- ������������µľ�������
%   obj_fcn     ---- Ŀ�꺯��ֵ
mf = U.^expo;       % �����Ⱦ������ָ��������
K = Kernel(data, center,segma2);        
K1 = Kernel(data1, center,segma2); 
center= sum(mf.*(K.*(ones(cluster_n,1)*data')+beta*K1.*(ones(cluster_n,1)*data1')),2)./sum(mf.*(K+beta*K1),2); % �¾�������
dist = 2*(1-K);
dist1 = 2*(1-K1);
alfa_j = min(dist.^2)*alfa;    %alfa_j  is control factor vector
obj_fcn = sum(sum((dist.^2).*mf))+ sum(alfa_j.*sum(U.*(1-U.^(expo-1))))+beta*sum(sum((dist1.^2).*mf));  % ����Ŀ�꺯��ֵ 
tmp = (dist.^2-ones(cluster_n,1)*alfa_j+beta*dist1.^2).^(-1/(expo-1));
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % �����µ������Ⱦ��� (5.3)ʽ
% tmp = dist.^(2/(expo-1));
% U_new =ones(cluster_n, size(data,1))./(tmp.*(ones(cluster_n,1)*sum(ones(cluster_n, size(data,1))./tmp)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �Ӻ���3
function out = distfcm(center, data)
% �������������������ĵľ���
% ���룺
%   center     ---- ��������
%   data       ---- ������
% �����
%   out        ---- ����
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1), % ��ÿһ����������
    % ÿһ��ѭ��������������㵽һ���������ĵľ���
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end

function segma2=estimateSegma(data)
%estimate the parameter �� in Kernel distance
% Input��
%   data           ----samples
%
% Output��
%   segma2       ----the parameter ��^2
mean_data = mean(data);
segma2 = (data-mean_data)'*(data-mean_data)/size(data,1);

function d=Kernel(data,center,segma2)
%compute Kernel distance of vector p and v
% Input��
%   p,v           ----vector p and v
%   segma2         ----the parameter ��^2
% Output��
%   d             ----Kernel distance
d = zeros(size(center, 1), size(data, 1));
data = data';
for k = 1:size(center, 1), % ��ÿһ����������
    % ÿһ��ѭ��������������㵽һ���������ĵľ���
    d(k, :) = exp(-(data-center(k)).^2/segma2);
end


