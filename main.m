close all;
clear all;
clc;

nii = load_nii( 'sub_TJ001_brain_FLIRT.nii.gz' );  % װ��.nii����
img = nii.img;  % ��Ϊ����ļ���img��head�������֣�����img������ͼ������
save image.mat img  % �����ݱ��mat��ʽ
load 'image.mat'  % ��������
[n1, n2, n3] = size(img);   % ��ȡ.nii�ļ�������ά�ȣ�һ��1��2ά��ͼ��ά�ȣ�����ά����Ƭ
v_p1 = 0; 
v_p2 = 0; 
v_p3 = 0; 
v_p4 = 0;
n_start = 1;
n_end =79;
for slicenum = n_start:n_end
data_load = imrotate(img(:,:,slicenum),90);
imshow(data_load,[]);  % �����������ʾ��45����Ƭ��ͼ��

%10??154��Ƭ
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
[a b]=size(data2);%ͼ�����Ĵ�С
data=data2(:);%ת��Ϊ������
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����

options=[2;150;1e-5;1];%ָ�������Ⱦ���ĵ�ģ��ָ�����㷨����������������ֹ������
cluster_n=4; %ָ�������
expo = options(1);          % �����Ⱦ���U��ָ��
max_iter = options(2);		% ���������� 
min_impro = options(3);		% ��������С�仯��,������ֹ����
display = options(4);		% ÿ�ε����Ƿ������Ϣ��־ 
alfa=1.5;

% [U, center,obj_fcn] = FCMClust2(data2, cluster_n, expo, max_iter, min_impro, display);
%[U, center,obj_fcn] = FCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
%[U, center,obj_fcn] = KFCMClust2brainweb(data_load, cluster_n, expo, max_iter, min_impro, display);
%[U, center,obj_fcn] = KFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
[U, center,obj_fcn] = MyKFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);

fm = U.^expo;
V_pc = sum(fm(:))/data_n;
V_xieben = XieBeniInverted(U, expo, center, data);

%ͼ��ָ�
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
    result{k} = Ik;%result{k}��¼�ڣ����е����ؼ���λ����Ϣ
    wholeG(indexk) = x(k);%wholeG��¼�ָ������ͼ�����Ϣ
end

wholeG=reshape(wholeG,a,b);%������ת��Ϊa*b��С�ľ���
% [label_new,accuracy]=succeed(real_label,cluster_n,wholeG);
% wholeG = reshape(label_new,a,b);
%��ʾ�ָ���
subplot(1,2,1);imshow(data_load,[]);title('ԭʼͼ��');
subplot(1,2,2);imshow(wholeG,[]);title('�ָ��Ľ��');
% 
% for i=1:cluster_n
%    figure;
%    imshow(result{x(i)});
% end
%�������
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


