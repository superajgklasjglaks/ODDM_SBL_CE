% function [x, index] = OMP(A,b,sparsity)
% %Step 1
% index = []; k = 1; [Am, An] = size(A); r = b; x=zeros(An,1);
% cor = A'*r; 
% while k <= sparsity
%     %Step 2
%     [Rm,ind] = max(abs(cor)); 
%     index = [index ind]; 
%     %Step 3
%     P = A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
%     r = (eye(Am)-P)*b; cor=A'*r;
%     k=k+1;
% end
% %Step 5
% xind = inv(A(:,index)'*A(:,index))*A(:,index)'*b;
% x(index) = xind;
% end

function [x, index, k_v_omp, l_t_omp] = OMP(x_p,x_kp,x_lp,M,N,N_T,M_T,y_T,r_v,r_t,k_max,l_max)
N_v = ceil(2 * k_max / r_v);   % Doppler domain virtual sampling grid
M_t = ceil(l_max / r_t);   % delay domain virtual sampling grid
virtual_size = N_v * M_t;
rho = 1e-2;
c = 1e-4;
d = 1e-4;

%% begin off-gird CE
k_v = zeros(virtual_size,1);
l_t = zeros(virtual_size,1);
[kv_bar,lt_bar] = First_Order_Linear_Approximation(N_v,M_t,k_max,r_v,r_t);
[phi_trunc,phi_t,phi_t_v,phi_t_t] = Gen_measurement_matrix(x_p,x_kp,x_lp,M,N,M_T,N_T,M_t,N_v,kv_bar,lt_bar,k_v,l_t);
A = phi_t + phi_t_v * diag(k_v) + phi_t_t * diag(l_t);

% [x, index] = OMP_inner(y_T, A, 5);

%Step 1
index = []; k = 1; [Am, An] = size(A); r = y_T; x=zeros(An,1); sparsity = 20;
cor = A'*r; 
while k <= sparsity
    %Step 2
    [Rm,ind] = max(abs(cor)); 
    index = [index ind]; 
    %Step 3
    P = A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
    r = (eye(Am)-P)*y_T; cor=A'*r;
    k=k+1;
end
%Step 5
xind = inv(A(:,index)'*A(:,index))*A(:,index)'*y_T;
x(index) = xind;

k_v_omp = kv_bar;
l_t_omp = lt_bar;
end

function [h_i, Omega] = OMP_inner(y,P,S)      
 N=size(y,1);
 K= size(P,2);

 err= 1e-3;  
h_i = zeros(K, 1);  % Initialize the recovery signal
y_i = y - P * h_i;    % residual error

it=0;
stop = 0;
P_s = [];
pos_array=[];
while ~stop                     
           for ntl=1:1:K
                product(ntl)=abs(P(:,ntl)'*y_i);      %  Step 1): calculate the correlation between the (1+(l-1)*Nt)-th  to (l*Nt)-th  columns and the residual signal 
           end
            
             y_i_before=y_i;    
            [val,pos]=max(product);                       %  Step 2): find the value and position of the most matching columns
           
            P_s=[P_s,P(:,pos)];  pos_array= [pos_array,pos];  %  Step 3): update the argment matrix and record the position,
            P(:,pos)=zeros(N,1);                        %           remove the column just selected
            h_s=(P_s'*P_s)^(-1)*P_s'*y; 
%             h_s = P_s\y;                       %  Step 4): solve the LS problem to obatain a new signal estimate
            y_i = y - P_s*h_s;                          %  Step 5): Calculate the new residual
                                              
           it = it + 1;                                  %Iteration counter              
                                                     
    %Check Halting Condition
    if (it >= S  )  ||  norm(y_i) <=err*norm(y) 
      stop = 1;
    end
%     ||  norm(y_i) <=err*norm(y) 
%     ||norm( y_i) >= norm( y_i_before)
end
  h_i(pos_array) = (P_s'*P_s)^(-1)*P_s'*y;    %  Step 7): get the recovered signal
 %End CoSaMP iteration
%  used_iter = it;
  Omega=pos_array;     
end