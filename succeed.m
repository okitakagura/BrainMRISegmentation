function [label_new, accuracy]=succeed(real_label,K,result)
%输入K：聚的类，id：训练后的聚类结果，N*1的矩阵
[m,n]=size(result);
id = reshape(result,m*n,1);
real_label = reshape(real_label,m*n,1)+ones(m*n,1);
N=size(id,1);   %样本个数
p=perms(1:K);   %全排列矩阵
p_col=size(p,1);   %全排列的行数
new_label=zeros(N,p_col);   %聚类结果的所有可能取值，N*p_col
num=zeros(1,p_col);  %与真实聚类结果一样的个数
%将训练结果全排列为N*p_col的矩阵，每一列为一种可能性
for i=1:N
    for j=1:p_col
        for k=1:K
            if id(i)==k
                new_label(i,j)=p(j,k);  %1 2 3
            end
        end
    end
end
%与真实结果比对，计算精确度
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