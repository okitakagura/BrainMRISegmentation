function[U, center,obj_fcn] = KFCMClust2brainweb(data_load, cluster_n, expo, max_iter, min_impro, display)
% KFCMClust2brainweb.m   ���ý�Ϻ˺�����ģ��C��ֵ�����ݼ�data��Ϊcluster_n��   
% ���룺
%   data_load   ---- ���ݼ�
%   cluster_n     ---- ������Ŀ
%   expo  ---- �����Ⱦ���U��ָ����>1                  (ȱʡֵ: 2.0)
%   max_iter  ---- ����������                           (ȱʡֵ: 100)
%   min_impro  ---- ��������С�仯��,������ֹ����           (ȱʡֵ: 1e-5)
%   display  ---- ÿ�ε����Ƿ������Ϣ��־                (ȱʡֵ: 1)
% �����
%   center      ---- ��������
%   U           ---- �����Ⱦ���
%   obj_fcn     ---- Ŀ�꺯��ֵ


data1=data_load;
[a b]=size(data1);
data=data1(:);
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����

kernel_b=estimateSegma(data);
obj_fcn = zeros(max_iter, 1);	% ��ʼ���������obj_fcn
U = initkfcm(cluster_n, data_n);     % ��ʼ��ģ���������,ʹU�����������Ϊ1

% ��ʼ���������ģ����������ݵ�������ѡȡcluster_n��������Ϊ�������ġ���Ȼ��
% �������ĳЩ����֪ʶѡȡ���Ļ����ܹ��ﵽ�ӿ��ȶ���Ч������Ŀǰ���߱��⹦��
index = randperm(data_n);   % �����������������
center_old = data(index(1:cluster_n),:);  % ѡȡ������е�������ǰcluster_n��

% Main loop  ��Ҫѭ��
for i = 1:max_iter,
    %�ڵ�k��ѭ���иı��������ceneter,�ͷ��亯��U��������ֵ;
	[U, center, obj_fcn(i)] = stepkfcm(data,U,center_old, expo, kernel_b);
	if display, 
		fprintf('KFCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    center_old = center;    % ���µľ������Ĵ����ϵľ�������
	% ��ֹ�����б�
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% ʵ�ʵ������� 
obj_fcn(iter_n+1:max_iter) = [];



%�Ӻ���
function U = initkfcm(cluster_n, data_n)
% ��ʼ��fcm�������Ⱥ�������
% ����:
%   cluster_n   ---- �������ĸ���
%   data_n      ---- ��������
% �����
%   U           ---- ��ʼ���������Ⱦ���
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);

% �Ӻ���
function [U_new,center_new,obj_fcn] = stepkfcm(data,U,center,expo,kernel_b)
% ģ��C��ֵ����ʱ������һ��
% ���룺
%   data        ---- nxm����,��ʾn������,ÿ����������m��ά����ֵ
%   U           ---- �����Ⱦ���
%   center      ---- ��������
%   expo        ---- �����Ⱦ���U��ָ��         
%   kernel_b    ---- ��˹�˺����Ĳ���
% �����
%   U_new       ---- ������������µ������Ⱦ���
%   center_new  ---- ������������µľ�������
%   obj_fcn     ---- Ŀ�꺯��ֵ
feature_n = size(data,2);  % ����ά��
cluster_n = size(center,1); % �������
mf = U.^expo;       % �����Ⱦ������ָ�����㣨c��n��)

% �����µľ�������;
KernelMat = gaussKernel(center,data,kernel_b); % �����˹�˾���(c��n��)
num = mf.*KernelMat * data;   
den = sum(mf.*KernelMat,2); 
center_new = num./(den*ones(1,feature_n)); % �����µľ�������(c��p��,c������)

% �����µ������Ⱦ���
kdist = distkfcm(center_new, data, kernel_b);    % ����������
obj_fcn = sum(sum((kdist.^2).*mf));  % ����Ŀ�꺯��ֵ
tmp = kdist.^(-1/(expo-1));     
U_new = tmp./(ones(cluster_n, 1)*sum(tmp)); 

% �Ӻ���
function out = distkfcm(center, data, kernel_b)
% �������������������ĵľ���
% ���룺
%   center     ---- ��������
%   data       ---- ������
% �����
%   out        ---- ����
cluster_n = size(center, 1);
data_n = size(data, 1);
out = zeros(cluster_n, data_n);
for i = 1:cluster_n % ��ÿ���������� 
    vi = center(i,:);
    out(i,:) = 2-2*gaussKernel(vi,data,kernel_b);
end

function segma2=estimateSegma(data)
% �����˹����Ĳ���
% ���룺
%   data           ---- ���ݵ�
%
% �����
%   segma2       ---- ������^2
mean_data = mean(data);
segma2 = (data-mean_data)'*(data-mean_data)/size(data,1);

% �Ӻ���
function out = gaussKernel(center,data,kernel_b)
% ��˹�˺�������
% ����:
%   center      ---- ģ����������
%   data        ---- �������ݵ�
%   kernel_b    ---- ��˹�˲���
% �����
%   out         ---- ��˹�˼�����

dist = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1), % ��ÿһ����������
    % ÿһ��ѭ��������������㵽һ���������ĵľ���
    dist(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end
out = exp(-dist.^2/kernel_b);