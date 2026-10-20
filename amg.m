clc;
clear;


%parameters
n = 60;
A1 = gallery('poisson', n);
tol=1e-6;
AMG = amg_setup(A1);
b = ones(size(A1,1),1);
[x,res] = amg_solve(AMG, b, 200,tol);
X = reshape(x, n, n);
figure;
contourf(X, 20);
colorbar;
title('AMG Solution Contours');
figure;
semilogy(res, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
title('AMG Residual History');
xlabel('Cycle');
ylabel('||r||');
grid on;


%strength connection function 
function [S,ST]=connections(A,theta)
% find the connections form matrix A
% 
% S{i}: strong variables that influence node i
% ST{i}: strong variables that are influenced by node i
% 
% S{i}=-aij>=theta*max(-aik)

%rows of A
n=size(A,1);

%matrices to store the nodes combinations 
S=cell(n,1);
ST=cell(n,1);
for i = 1:n
    S{i}  = [];
    ST{i} = [];
end

for i=1:n
    
    %get row 
    [~, cols, vals] = find(A(i, :));
    
    %ignore the aii
    mask = (cols ~= i)& (vals < 0);
    cols = cols(mask);
    vals = vals(mask);
    if isempty(vals), continue; end
 
    %form the -aij,j!=i
    strength = -vals;
    
    %form the max(-aik)
    alpha = max(strength);

    %find the connection -aij>=theta*max(-aik)
    strong_idx = cols(strength >= theta * alpha);
    S{i}=strong_idx;

    for j=S{i}
        ST{j}=[ST{j},i];
    end
end
end




function [C,F]=coarsening(S,ST)
 % Runge-Stuben C/F splitting 
 % lambda(i)=|ST(i)|
 % C points: coarse variables 
 % F points: fine variables 
 % states:0 undecided, 1 coarse , -1 fine
   
 n=length(S);
 state=zeros(n,1); %vector to st ore the states of each node
 lambda=zeros(n,1); %vector that will store the length if each ST{i}

 for i=1:n
    lambda(i)=length(ST{i});
 end

 %choose the points that will become coarse and fine
 while any(state==0)
  %make copy of lambda
  score=lambda; 
  score(state~=0)=-1; %initialize decided values to infinity for 
  % search start

  %get max lambda value
  [~,i]=max(score);
  state(i)=1;  %the i element of S becomes coarse
    
  %the dependent points of the i element becomes fine 
  % for j=S{i}
  for j=S{i}
      if state(j)==0
        state(j)=-1; %convert state to -1 ->fine

         % increase importance of neighbours
          % for k=ST{j}
          for k=ST{j}
              if state(k)==0
                  lambda(k)=lambda(k)+1;
              end
          end
      end
  end
 end
for i = find(state == -1)'
    for j = intersect(S{i}, find(state == -1))'
        % Check if i and j share a C-point
        if isempty(intersect(S{i}, intersect(S{j}, find(state == 1))))
            state(i) = 1; % Convert i to coarse point
            break;
        end
    end
end
C=find(state==1);

F=find(state==-1);

end


%interpolation function
function P=interpolate(A,S,C,F) 
% P=[I;W]
% coarse: P(i,i)=1
% fine: e_i=sum of wij*ej for j in Ci

% wij=-up/down
% up=aij+sum_{m belongs to F_i}(aim*amj/sum_{k belongs to C_i}(amk))
% down=aii+sum_{m belongs to F_i}(aim)

n=size(A,1);
nc=length(C);
P=sparse(n,nc); %P nxnc

%mapping coarse global indexes to P indexes
map=zeros(n,1);
for k=1:nc
    map(C(k))=k;
    P(C(k),k)=1;%for coarse points Pii=1
end

%interpolate wij
for i=F(:)'
    Fi=intersect(S{i},F);  %for m belongs to Fi^s
    Ci=intersect(S{i},C);  %for k belongs to Ci

     if isempty(Ci)
        continue
     end

    % Fi_w=setdiff(F,Fi); %Fi^w =F-Fi^s for m belongs to Fi^w
    % 
    % down=A(i,i)+sum(A(i,Fi_w));
    % Fi_weak=setdiff(F,Fi);
    % down=A(i,i)+sum(A(i,Fi_weak));
    off_diag = find(A(i,:));
    off_diag(off_diag == i) = []; 
    weak_nodes = setdiff(off_diag, S{i});

% Denominator includes A(i,i) plus weak off-diagonals
down = A(i,i) + sum(A(i, weak_nodes));
    if abs(down)<eps
        continue;
    end

    for j=Ci(:)'
         up=A(i,j);
        for m=Fi(:)'
            denom2 = sum(A(m,Ci));
          
            if abs(denom2) > eps
                up = up + A(i,m)*A(m,j)/denom2;
            end
        end 
        % wij=-up/down;
        P(i,map(j))=-up/down; 
    end
end
end



%function that implemnets the galerkin approach 
function Ac=galerkin(A,P)
    Ac=P'*A*P;
%     disp("this is a message");
%     d = full(diag(Ac));
% 
% fprintf("size(Ac) = %d\n", size(Ac,1));
% fprintf("nnz(Ac) = %d\n", nnz(Ac));
% fprintf("min(diag) = %g\n", min(d));
% fprintf("max(diag) = %g\n", max(d));
end

function AMG=amg_setup(A)
    theta=0.5;
     % max_levels = 5;
    AMG.A{1} = A;
    level = 1;
    
    while size(AMG.A{level},1) > 50 && level<5
        Alevel = AMG.A{level};
        fprintf("Level %d: size=%d nnz=%d\n", ...
            level,size(Alevel,1),nnz(Alevel));
        [S, ST] = connections(Alevel, theta);
        [C,F]   = coarsening(S,ST);
        fprintf("Fine/Coarse ratio: F=%d, C=%d\n", numel(F), numel(C));
        P       = interpolate(Alevel,S,C,F);
        Ac      = galerkin(Alevel,P);
        
        fprintf("Level %d: Size = %d -> %d coarse nodes\n",...
            level, size(Alevel,1), length(C));

        AMG.P{level}   = P;
        AMG.A{level+1} = Ac;

        level = level + 1;
    end
    AMG.levels=level;
end


function x=gauss(A,x,b,iter)
    DL=tril(A);
    for i=1:iter
        x=x+DL\(b-A*x);
    end
end


function x = V(AMG, level, x, b)

    A = AMG.A{level};

    % base case: coarsest level
    if level == AMG.levels
        x = A \ b;
        return
    end

    % pre-smoothing
    x = gauss(A, x, b, 4);

    % residual
    r = b - A*x;

    % restrict residual
    P   = AMG.P{level};
    r_c = P' * r;

    % coarse correction
    e_c = zeros(size(r_c));
    e_c = V(AMG, level+1, e_c, r_c);

    % % prolong correction
    % e = P * e_c;
    % 
    % % correct
    % x = x + e;
    x = x + AMG.P{level} * e_c;

    % post-smoothing
    x = gauss(A, x, b, 3);
end


function [x,res] = amg_solve(AMG, b, cycles,tol)
    n = size(AMG.A{1},1);
    x = zeros(n,1);
    bnorm = norm(b);
    res=zeros(cycles,1);
    for k = 1:cycles
        x = V(AMG, 1, x, b);
        r = b - AMG.A{1} * x;
        rnorm = norm(r);
        res(k)=rnorm;
        if mod(k,10)==0
            fprintf('Iteration : %d with Residual: %f\n',k,res(k));
        end
        % stopping criterion
        if rnorm < tol * bnorm
            res(k+1:end)=[];
            fprintf('Converged at cycle %d\n', k);
            fprintf('Residual: %g\n',res(k));
            break
        end
    end
end
