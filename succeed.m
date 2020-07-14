function [label_new, accuracy]=succeed(real_label,K,result)
%����K���۵��࣬id��ѵ����ľ�������N*1�ľ���
[m,n]=size(result);
id = reshape(result,m*n,1);
real_label = reshape(real_label,m*n,1)+ones(m*n,1);
N=size(id,1);   %��������
p=perms(1:K);   %ȫ���о���
p_col=size(p,1);   %ȫ���е�����
new_label=zeros(N,p_col);   %�����������п���ȡֵ��N*p_col
num=zeros(1,p_col);  %����ʵ������һ���ĸ���
%��ѵ�����ȫ����ΪN*p_col�ľ���ÿһ��Ϊһ�ֿ�����
for i=1:N
    for j=1:p_col
        for k=1:K
            if id(i)==k
                new_label(i,j)=p(j,k);  %1 2 3
            end
        end
    end
end
%����ʵ����ȶԣ����㾫ȷ��
for j=1:p_col
    for i=1:N
        if new_label(i,j)==real_label(i)
                num(j)=num(j)+1;
        end
    end
end
[M,I]=max(num);
accuracy=M/N;
label_new=new_label(:,I);