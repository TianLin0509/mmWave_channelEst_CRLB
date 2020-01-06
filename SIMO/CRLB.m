Nr =64;
Nt = 1;
Nrf = 1;
Ns = 1;
SNR = 5:5:10;
for f = 1 : 3
    Nloop = 50000;
    % random observation matrices
    Nw = 10; % number of observations
    Nv = 1;
    W = exp( 1i*unifrnd(0,2*pi,Nr,Nw));
    theta = unifrnd(-pi/2, pi/2);
    alpha = sqrt(1/2) * (randn() + 1j * randn());
    for j = 1 : length(SNR)
        crlb_mse = zeros(1, Nloop);
        music_mse = zeros(1, Nloop);
        for ll = 1 : Nloop
            snr = SNR(j);
            Vn = 1/ 10 ^ (snr / 10);
            % compute phi
            phi = W';
            
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
            r = y * y';
            [U, S, V] = svd(r);
            NoiseSpace = U(:, 2:end);
            
            theta_set =  -pi/2:0.01:pi/2;
            for i = 1 : length(theta_set)
                theta_tmp = theta_set(i);
                x = exp(1j * (0:(Nr-1))'* pi * sin(theta_tmp));
                tmp =  NoiseSpace' * W' * x;
                p(i) = norm(tmp, 'fro');
            end
            [~, idx] = min(p);
            theta_est = theta_set(idx);
            music_mse(ll) = (theta - theta_est)^2;
        end
        Music(j) = mean(music_mse);
        crlb(j) = mean(crlb_mse);
    end
    Music
    crlb
    clear Music crlb
end

semilogy(SNR, crlb, 'k', 'LineWidth', 2)
hold on
semilogy(SNR, Music, 'b', 'LineWidth', 2)
legend('CRLB', 'Music')
xlabel('SNR')
ylabel('MSE')

