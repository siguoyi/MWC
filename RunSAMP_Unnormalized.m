function [ theta, Pos_theta ] = RunSAMP_Unnormalized( y,A,S )
% SAMP alogrithm for MMV/SMV: Y=AX
% This function returns a solution Xr together with the support (Supp)
% The termination criteria are:
%    Stop if more than NumIters iterations
%    Stop if residual norm is greater than ResThreshold
%    Stop if the ratio between residual norm and the norm of the current solution is greater than ResvsSolThreshold
% If SymmetricSupport is true, then the algorithm select pairs of indices every
%      iteration, such that the indices are symettric with respect to solution
%       dimensions (special feature for conjugate symmetric solutions).
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [M,N] = size(A);%���о���AΪM*N����
    theta = zeros(N,1);%�����洢�ָ���theta(������)
    Pos_theta = [];%�������������д洢A��ѡ��������
    r_n = y;%��ʼ���в�(residual)Ϊy
    L = S;%��ʼ������(Size of the finalist in the first stage)
    Stage = 1;%��ʼ��Stage
    IterMax = M;
    NormACols = sqrt(diag(A'*A));
    for ii=1:IterMax%������M��
        %(1)Preliminary Test
        product = A'*r_n;%���о���A������в���ڻ�
        product_1 = sqrt(sum(abs(product).^2,2))./NormACols;
        [val,pos]=sort(abs(product_1),'descend');%��������
        Sk = pos(1:L);%ѡ������L��
        %(2)Make Candidate List
        Ck = union(Pos_theta,Sk);
        %(3)Final Test
        if length(Ck)<=M
            At = A(:,Ck);%��A���⼸����ɾ���At
        else
            theta_ls=0;
            break;
        end
        %y=At*theta��������theta����С���˽�(Least Square)
        theta_ls = (At'*At)^(-1)*At'*y;%��С���˽�
        [val,pos]=sort(abs(theta_ls),'descend');%��������
        F = Ck(pos(1:L));
        %(4)Compute Residue
        %A(:,F)*theta_ls��y��A(:,F)�пռ��ϵ�����ͶӰ
        theta_ls = (A(:,F)'*A(:,F))^(-1)*A(:,F)'*y;
        r_new = y - A(:,F)*theta_ls;%���²в�r
        if norm(r_new)<1e-6%halting condition true 
            Pos_theta = F;
            %r_n = r_new;%����r_n�Ա�������µ�r_n
            break;%quit the iteration
        elseif norm(r_new)>=norm(r_n)%stage switching 
            Stage = Stage + 1;%Update the stage index 
            L = Stage*S;%Update the size of finalist
            if ii == IterMax%���һ��ѭ��
                Pos_theta = F;%����Pos_theta����theta_lsƥ�䣬��ֹ����
            end
            %ii = ii - 1;%��������������
        else
            Pos_theta = F;%Update the finalist Fk
            r_n = r_new;%Update the residue
        end
    end
%     theta(Pos_theta)=theta_ls;%�ָ�����theta
end