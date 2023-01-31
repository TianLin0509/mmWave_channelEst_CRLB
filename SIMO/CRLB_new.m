Nr =192;
Nt = 1;
Nrf = 1;
Ns = 1;
% SNR = -20 : 5 : 10;
SNR = 5:5:20;

% Nloop = 1000;
Nloop = 10;
% random observation matrices
Nw = 256; % number of observations
Nv = 10;
% W = exp( 1i*unifrnd(0,2*pi,Nr,Nw));
W = eye(Nr);
theta = unifrnd(-pi/2, pi/2);
alpha = sqrt(1/2) * (randn() + 1j * randn());
Music = zeros(length(SNR));
crlb  = zeros(length(SNR));
for j = 1 : length(SNR)
    fprintf('SNR = %d\n', SNR(j))
    crlb_mse = zeros(1, Nloop);
    music_mse = zeros(1, Nloop);
    for ll = 1 : Nloop
        if ll == 10
            fprintf('迭代10次了\n')
        end
        snr = SNR(j);
        Vn = 1/ 10 ^ (snr / 10);
        V = ones(1, Nv);
        % compute phi
        phi = kron(V', W');
        
        % compute A
        index = (0 : (Nr - 1))';
        a1 = exp(1j * pi * index * sin(theta));
        a2 = a1 * 1j;
        a3 = zeros(Nr,1);
        for m = 1 : Nr
            a3(m) = alpha * (m-1) * 1j * pi * cos(theta) * exp(1j*pi*(m-1)*sin(theta));
        end
        A = [a1 a2 a3];
        %A = a3;
        F = 2 * Nrf / Nr / Vn * real(A'*phi'*phi*A);
        C = F^(-1);
        
        crlb_mse(ll) = C(3,3);
        %crlb_mse(ll) = C(1,1);
        %music
        s = ones(1, Nv);
        h = alpha *a1;
        n = Vn / sqrt(2) * (randn(Nr,Nv) + 1j * randn(Nr,Nv));
        
        y = W' * h * s + W' * n;
%         r = y * y';
%         [U, S, V] = svd(r);
%         NoiseSpace = U(:, 2:end);
%         
%         theta_set =  -pi/2:0.01:pi/2;
%         for i = 1 : length(theta_set)
%             theta_tmp = theta_set(i);
%             x = exp(1j * (0:(Nr-1))'* pi * sin(theta_tmp));
%             tmp =  NoiseSpace' * W' * x;
%             p(i) = norm(tmp, 'fro');
%         end
%         [~, idx] = min(p);        
%         theta_est = theta_set(idx);
        theta_est = MUSIC(y);
        music_mse(ll) = (theta - asin(theta_est))^2;
    end
    Music(j) = mean(music_mse);
    crlb(j) = mean(crlb_mse);
end

figure; %hold on; box on; grid on;
plot(SNR, crlb(:,1), 'r:s', 'LineWidth', 2)
hold on;
plot(SNR, Music, 'g-<', 'LineWidth', 2)
xlabel('SNR')
ylabel('MSE')
set(gca, 'YScale', 'log');
legend({  
        'CRLB',...    
        'MUSIC'},...    
        'Interpreter','latex', 'Box','off'); 