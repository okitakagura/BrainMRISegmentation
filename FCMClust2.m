%%%��׼��fcm�㷨
function [U, center,obj_fcn] = FCMClust2(data2, cluster_n, expo, max_iter, min_impro, display)
% FCMClust.m   ����ģ��C��ֵ�����ݼ�data��Ϊcluster_n��   
% ���룺
%   data2   ---- ���ݼ�
%   cluster_n     ---- ������Ŀ
%   expo  ---- �����Ⱦ���U��ָ����>1                  (ȱʡֵ: 2.0)
%   max_iter  ---- ����������                           (ȱʡֵ: 100)
%   min_impro  ---- ��������С�仯��,������ֹ����           (ȱʡֵ: 1e-5)
%   display  ---- ÿ�ε����Ƿ������Ϣ��־                (ȱʡֵ: 1)
% �����
%   center      ---- ��������
%   U           ---- �����Ⱦ���
%   obj_fcn     ---- Ŀ�꺯��ֵ


[a b]=size(data2);%ͼ�����Ĵ�С
data=data2(:);%ת��Ϊ������
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����
obj_fcn = zeros(max_iter, 1);	% ��ʼ���������obj_fcn
U = initfcm(cluster_n, data_n);     % ��ʼ��ģ���������,ʹU�����������Ϊ1,

% Main loop  ��Ҫѭ��
for i = 1:max_iter,
    %�ڵ�k��ѭ���иı��������ceneter,�ͷ��亯��U��������ֵ;
	[U, center, obj_fcn(i)] = stepfcm(data, U, cluster_n, expo);
	if display
		fprintf('FCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	end
	% ��ֹ�����б�
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, 
            break;
        end,
	end
end

iter_n = i;	% ʵ�ʵ������� 
obj_fcn(iter_n+1:max_iter) = [];


% �Ӻ���:���ڳ�ʼ�������Ⱦ���U
function U = initfcm(cluster_n, data_n)
% ��ʼ��fcm�������Ⱥ�������
% ����:
%   cluster_n   ---- �������ĸ���
%   data_n      ---- ��������
% �����
%   U           ---- ��ʼ���������Ⱦ���
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);%��ʼ��ʱ�������һ������ÿһ�к�Ϊ1�������Ⱦ���


% �Ӻ��������ڼ���ÿ�ε���ÿһ���в������µ������Ⱦ���u�;�������v
function [U_new, center, obj_fcn] = stepfcm(data, U, cluster_n, expo)
% ģ��C��ֵ����ʱ������һ��
% ���룺
%   data        ---- nxm����,��ʾn������,ÿ����������m��ά����ֵ
%   U           ---- �����Ⱦ���
%   cluster_n   ---- ����,��ʾ�ۺ�������Ŀ,�������
%   expo        ---- �����Ⱦ���U��ָ��                      
% �����
%   U_new       ---- ������������µ������Ⱦ���
%   center      ---- ������������µľ�������
%   obj_fcn     ---- Ŀ�꺯��ֵ
mf = U.^expo;       % �����Ⱦ������ָ��������
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))'); % �¾�������
dist = distfcm(center, data);       % ����������
obj_fcn = sum(sum((dist.^2).*mf));  % ����Ŀ�꺯��ֵ 
tmp = dist.^(-2/(expo-1));     
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % �����µ������Ⱦ���

% �Ӻ��������������㵽ÿһ���������ĵ�ŷ�Ͼ���
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
   % out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2),2)',1);
     out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end


