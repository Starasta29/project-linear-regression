clc
clear all
load proj_fit_01.mat

% we load the data for x and y id from the data set
x1_id = id.X{1,1};
x2_id = id.X{2,1};
y_id = id.Y;

% we load the data for x and y val from the data set
x1_val = val.X{1,1};
x2_val = val.X{2,1};
y_val = val.Y;

% we reshape and transpose y_id to get a column matrix
y_reshape = reshape(y_id,1,[]);
y_reshape_transposed = y_reshape';

% we reshape and transpose y_val to get a column matrix
y_reshape_val = reshape(y_val,1,[]);
y_reshape_transposed_val = y_reshape_val';

% counter for the MSE vector
count = 1;
% we check what is the best n for our data
for n = 1:30
    % initializing pascal matrix and phi matrix dimensions
    X = zeros(n);
    P = abs(pascal(n,1));
    phi_id = zeros(length(x1_id)*length(x2_id),nnz(P));
    
    pos = 1;

    for e1 = 1:length(x1_id) 
        for e2 = 1:length(x2_id)
            X = zeros(n);
            % we construct the polynomial line for one x1id and x2id
            % combination
            for i = 1:n
                p2 = 0;
                for p1 = i-1:-1:0
                    X(i,p2+1) = X(i,p2+1) + x1_id(e1)^(p1) * x2_id(e2)^(p2);
                    if X(i,p2+1) == 0
                        X(i,p2+1) = 's';
                    end
                    p2 = p2 + 1;
                end
            end
            X = transpose(X);
            X = X(:);
            X = transpose(X);
            X = X(X~=0);
            X(X=='s') = 0;
            % we add the polynomial line into the phi matrix
            for i = 1:length(X)
                phi_id(pos,i) = phi_id(pos,i) + X(i);
            end
            pos = pos + 1;
        end
    end

% initializing pascal matrix and phi matrix dimensions
X = zeros(n);
P = abs(pascal(n,1));
phi_val = zeros(length(x1_val)*length(x2_val),nnz(P));
pos = 1;

for e1 = 1:length(x1_val) 
    for e2 = 1:length(x2_val)
        X = zeros(n);
        % we construct the polynomial line for one x1val and x2val
        % combination
        for i = 1:n
            p2 = 0;
            for p1 = i-1:-1:0
                X(i,p2+1) = X(i,p2+1) + x1_val(e1)^(p1) * x2_val(e2)^(p2);
                if X(i,p2+1) == 0
                    X(i,p2+1) = 's';
                end
                p2 = p2 + 1;
            end
        end
        X = transpose(X);
        X = X(:);
        X = transpose(X);
        X = X(X~=0);
        X(X=='s') = 0;
        % we add the polynomial line into the phi matrix
        for i = 1:length(X)
            phi_val(pos,i) = phi_val(pos,i) + X(i);
        end
        pos = pos + 1;
    end
end
    % using the formulas from the course to calculate theta and y_hat
    theta = phi_id \ y_id(:);
    y_hat = phi_id * theta; 
    y_hat_val = phi_val * theta;
    MSE(count) = (1/1681)*sum((y_hat-y_reshape_transposed).^2);
    MSE_val(count) = (1/5041)*sum((y_hat_val-y_reshape_transposed_val).^2);
    count = count + 1;
end

% now we plot for the best n = 23
for n = 1:23
    % initializing pascal matrix and phi matrix dimensions
    X = zeros(n);
    P = abs(pascal(n,1));
    phi_id = zeros(length(x1_id)*length(x2_id),nnz(P));
    pos = 1;
    for e1 = 1:length(x1_id) 
        for e2 = 1:length(x2_id)
            X = zeros(n);
            % we construct the polynomial line for one x1id and x2id
            % combination
            for i = 1:n
                p2 = 0;
                for p1 = i-1:-1:0
                    X(i,p2+1) = X(i,p2+1) + x1_id(e1)^(p1) * x2_id(e2)^(p2);
                    if X(i,p2+1) == 0
                        X(i,p2+1) = 's';
                    end
                    p2 = p2 + 1;
                end
            end
            X = transpose(X);
            X = X(:);
            X = transpose(X);
            X = X(X~=0);
            X(X=='s') = 0;
            % we add the polynomial line into the phi matrix
            for i = 1:length(X)
                phi_id(pos,i) = phi_id(pos,i) + X(i);
            end
            pos = pos + 1;
        end
    end

% initializing pascal matrix and phi matrix dimensions
X = zeros(n);
P = abs(pascal(n,1));
phi_val = zeros(length(x1_val)*length(x2_val),nnz(P));
pos = 1;

for e1 = 1:length(x1_val) 
    for e2 = 1:length(x2_val)
        X = zeros(n);
        % we construct the polynomial line for one x1id and x2id
        % combination
        for i = 1:n
            p2 = 0;
            for p1 = i-1:-1:0
                X(i,p2+1) = X(i,p2+1) + x1_val(e1)^(p1) * x2_val(e2)^(p2);
                if X(i,p2+1) == 0
                    X(i,p2+1) = 's';
                end
                p2 = p2 + 1;
            end
        end
        X = transpose(X);
        X = X(:);
        X = transpose(X);
        X = X(X~=0);
        X(X=='s') = 0;
        % we add the polynomial line into the phi matrix
        for i = 1:length(X)
            phi_val(pos,i) = phi_val(pos,i) + X(i);
        end
        pos = pos + 1;
    end
end
    % using the formulas from the course to calculate theta and y_hat
    theta = phi_id \ y_id(:);
    y_hat = phi_id * theta; 
    y_hat_val = phi_val * theta;
    count = count + 1;
end

% comparing the calculated y_hat with the given y_id for the best n
figure
mesh(reshape(y_hat,[41,41]),'FaceColor','green')
title('Mesh for n = ', 23)
hold
mesh(y_id,'FaceColor','red')

% comparing the calculated y_hat_val with the given y_val for the best n
figure
mesh(reshape(y_hat_val,[71,71]),'FaceColor','green')
title('Mesh for MSE = 0.3275; n = 23')
hold
mesh(y_val,'FaceColor','red')

% comparing the MSE vectors
figure
plot(MSE)
title('MSE_{id} vs MSE_{val}')
hold
plot(MSE_val)
legend('MSE_{id}', 'MSE_{val}')