function g = readrawb(filename, num)
% ����readrawb(filename, num)�еĵ�һ������filename������ȡ��rawb�ļ����ļ������ڶ�������num���ǵڶ����š�
fid = fopen(filename);
% ������ȡ181*217*181�����ݣ���ʱ��temp��һ������Ϊ181*217*181��������
% �Ƚ�rawb�е��������ݴ��ݸ�temp���飬Ȼ��tempreshape��ͼƬ����
temp = fread(fid, 181 * 217 * 181);
% ���԰��������һ��181*217�У�181�е����飬�������Ĵ��룬�����181��ͼƬ�����ݣ�ÿһ�ж�Ӧһ��ͼ��
% ����ͼƬ�����顣ͼƬ��images������ÿһ�б�ʾһ��ͼƬ��
images = reshape(temp, 181 * 217, 181);  
% ��ȡ�����еĵ�num�У��õ�������reshape��ͼƬԭ����������������181*217��
image = images(:, num);
image = reshape(image, 181, 217);
g = image;
fclose(fid);
end