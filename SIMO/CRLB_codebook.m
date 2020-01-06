Nr =192;
Nt = 1;
Nrf = 1;
Ns = 1;
SNR = -20:5:10;

Nloop = 1000;
% random observation matrices

Nv = 10;
Ncode = 256;
index = (0 : (Nr - 1))';

% generate codebook
theta_code = -pi/2 : pi/(Ncode-1) : pi/2;
W = exp(1j * pi * index * sin(theta_code));
theta = unifrnd(-pi/2, pi/2);
alpha = sqrt(1/2) * (randn() + 1j * randn());

for j = 1 : length(SNR)
    crlb_mse = zeros(1, Nloop);
    music_mse = zeros(1, Nloop);
    for ll = 1 : Nloop
        snr = SNR(j);
        Vn = 1/ 10 ^ (snr / 10);
        V = ones(1, Nv);
        % compute phi
        phi = kron(V', W');
        
        % compute A
        
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
       [~, idx] = max(mean(y,2));
        music_mse(ll) = (theta - theta_code(idx))^2;
    end
    Music(j) = mean(music_mse);
    crlb(j) = mean(crlb_mse);
end


semilogy(SNR, crlb, 'k', 'LineWidth', 2)
hold on
semilogy(SNR, Music, 'b', 'LineWidth', 2)
legend('CRLB', 'codebook')
xlabel('SNR')
ylabel('MSE')

