function DEMC_krig2

t=1;T=30000;

N=4;
d=2;

%% new bayesian code
A=[-20 0];
B=[40 15];

prior1 = @(N,d) [unifrnd(A,B); unifrnd(A,B); unifrnd(A,B); unifrnd(A,B)]; %prior
prior= @(kkk) heaviside(-kkk(1)+B(1))*heaviside(kkk(1)-A(1))*heaviside(-kkk(2)+B(2))*heaviside(kkk(2)-A(2));
pdf=@(x) likelihood_krig(x);

% x = nan(T,d); p_x = nan(T,1); % Preallocate memory for chain and density
% xp=nan(T,d);
imim=zeros(T,1);
%%
tic
gamma_RWM = 2.38/sqrt(2 * d); % Calculate jump rate (RWM based)
x = nan(T,d,N); p_x = nan(T,N); % Preallocate memory for chains and density
st=[0.5,0.5;0.5,0.5;0.5,0.5;0.5,0.5];
X = prior1(N,d); p_X = pdf(X); % Create initial population and compute density
X_initial=X;
x(1,1:d,1:N) = reshape(X',1,d,N); p_x(1,1:N) = p_X'; % Store initial position of chain and density
for i = 1:N
    R(i,1:N-1) = setdiff(1:N,i); 
end % R?matrix: ith chain, the index of chains for DE

%% t
for t =1:T, % Dynamic part: Evolution of N chains
g = 1*randsample([gamma_RWM 1],1,true,[0.9 0.1]); % Set gamma: 2.38/sqrt(2d) or 1, 90/10 ratio

for i = 1:N, % Create proposal and accept/reject
Xp=nan(N,d);
Xp(i,:)=[-1 -1];

% while Xp(i,1)<0 | Xp(i,2)<0
[~,draw] = sort(rand(N-1,N)); % Randomly permute [1,...,N?1] N times
r1 = R(i,draw(1,i)); % Derive r1
r2 = R(i,draw(2,i)); % Derive r2; r1 not equal r2 not equal i
Xp(i,1:d) = (X(i,1:d) +g*(X(r1,1:d)-X(r2,1:d))+ 1e-1 * randn(1,d)); % Create ith point with differential evolution
% end

p_Xp(i,1) = pdf(Xp(i,1:d)); % Calculate density of ith proposal
ttt(t,i)=pdf(Xp(i,1:d))-pdf(X(i,1:d));
alpha(t,i) = min(ttt(t,i),0); % Compute Metropolis ratio

idx = alpha(t,i) > log(rand); % Alpha larger than U[0,1] or not?

if prior([Xp(i,1:d)])==0;
idx=0;
end

if idx,
X(i,1:d) = Xp(i,1:d); 
p_X(i,1) = p_Xp(i,1); % True: Accept proposal
imim(t,i)=1;
imim(t,:)
end
end
x(t,1:d,1:N) = reshape(X',1,d,N); p_x(t,1:N) = p_X'; % Add current position and density to chain
end
toc 
%% 
