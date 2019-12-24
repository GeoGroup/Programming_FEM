function [ g ] = num_to_g(num,nf)
%This subroutine finds the g vector from num and nf.
% ���ݽڵ�ź����ɶȱ��ȷ����Ԫ���ɶ�����
 
 nod=size(num,1) ;
 nodof=size(nf,1);
 g=zeros(nod*nodof,1);
 for i=1:nod
   k=i*nodof;
   g(k-nodof+1:k)=nf(:,num(i));
 end



end

