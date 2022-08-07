function y = costFun(p,c,M)  % clear for me


N = length(p)/M; % 126/6 =21
P = reshape(p, M, N); % make all values related to pg1 in one row and etc


y = 0; % initial cost==0 

for i = 1:N


costs = c(1,1)*P(1,i).^2+c(1,2)*P(1,i)+c(1,3)+...
        c(2,1)*P(2,i).^2+c(2,2)*P(2,i)+c(2,3)+...  % for i=1......,21 it has different fuel power generation corrospond to each load i
        c(3,1)*P(3,i).^2+c(3,3)*P(3,i)+c(3,3);

    y = y + costs;

end

