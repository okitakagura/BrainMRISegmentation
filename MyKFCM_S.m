function [U, center,obj_fcn] = MyKFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa)
% MyKFCM_S.m   ���øĽ��Ľ��������Ϣ�ͺ˺�����ģ��C��ֵ�����ݼ�data��Ϊcluster_n��   
% ���룺
%   data_load   ---- ���ݼ�
%   cluster_n     ---- ������Ŀ
%   expo  ---- �����Ⱦ���U��ָ����>1                  (ȱʡֵ: 2.0)
%   max_iter  ---- ����������                           (ȱʡֵ: 100)
%   min_impro  ---- ��������С�仯��,������ֹ����           (ȱʡֵ: 1e-5)
%   display  ---- ÿ�ε����Ƿ������Ϣ��־                (ȱʡֵ: 1)
%   alfa  ---- ������ϢȨ��
% �����
%   center      ---- ��������
%   U           ---- �����Ⱦ���
%   obj_fcn     ---- Ŀ�꺯��ֵ
     

data2=data_load;%����ԭʼͼ��
[a b]=size(data2);%ͼ�����Ĵ�С

img1 = data2;
width = 3;  %�ֲ����ߴ�
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
data1=img1(:);%ת��Ϊ������
data1=double(data1);
data=data2(:);%ת��Ϊ������
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����


obj_fcn = zeros(max_iter, 1);   % ��ʼ���������obj_fcn     
U = initfcm(cluster_n, data_n);     % ��ʼ��ģ���������,ʹU�����������Ϊ1, 
mf = U.^expo;       % �����Ⱦ������ָ��������
%center = mf*(data+alfa*data1)./((ones(size(data, 2), 1)*sum(mf'))'*(1+alfa)); %��ʼ��������

%center = [0.9769;108.4991;212.6035;302.3221];

center=kmeans4(data_load);
segma=estimateSegma(data);

% %��ʼ���������ģ����������ݵ�������ѡȡcluster_n��������Ϊ�������ġ���Ȼ��
% % �������ĳЩ����֪ʶѡȡ���Ļ����ܹ��ﵽ�ӿ��ȶ���Ч������Ŀǰ���߱��⹦��
% index = randperm(data_n);   % �����������������
% center_old = data(index(1:cluster_n),:);  % ѡȡ������е�������ǰcluster_n��

% Main loop  ��Ҫѭ��
for i = 1:max_iter,
    %�ڵ�k��ѭ���иı��������ceneter,�ͷ��亯��U��������ֵ;
%     U0 = U;
    [U, center obj_fcn(i)] = stepfcm(data,data1,center,U, cluster_n, expo,alfa,segma);
    if display, 
        fprintf('MyKFCM_S:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
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

     
% �Ӻ���2
function [U_new, center, obj_fcn] = stepfcm(data,data1,center,U, cluster_n, expo,alfa,segma)
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
K = Kernel(data, center,segma);        
K1 = Kernel(data1, center,segma); 
%center = mf*(data+alfa*data1)./((ones(size(data, 2), 1)*sum(mf'))'*(1+alfa));
center= sum(mf.*(K.*(ones(cluster_n,1)*data')+alfa*K1.*(ones(cluster_n,1)*data1')),2)./sum(mf.*(K+alfa*K1),2); % �¾�������
dist = 2*(1-K);
dist1 = 2*(1-K1);
% dist = distfcm(center, data);       % ����������
% dist1 = distfcm(center,data1);
obj_fcn = sum(sum((dist.^2).*mf))+alfa*sum(sum((dist1.^2).*mf));  % ����Ŀ�꺯��ֵ 
tmp = (dist.^2+alfa*dist1.^2).^(-1/(expo-1));
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % �����µ������Ⱦ���

     
function segma2=estimateSegma(data)
% �����˹����Ĳ���
% ���룺
%   data           ---- ���ݵ�
%
% �����
%   segma2       ---- ������^2
mean_data = mean(data);
segma2 = (data-mean_data)'*(data-mean_data)/size(data,1);

function d=Kernel(data,center,segma2)
% �����˹����
% ���룺
%   data           ---- ���ݵ�
%   center           ---- ��������
%   segma2         ---- ������
% �����
%   d             ---- ��˹����
d = zeros(size(center, 1), size(data, 1));
data = data';
for k = 1:size(center, 1), % ��ÿһ����������
    % ÿһ��ѭ��������������㵽һ���������ĵľ���
    d(k, :) = exp(-(data-center(k)).^2/segma2);
end

