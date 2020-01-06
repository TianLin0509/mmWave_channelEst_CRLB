function A = A_ULA(thetaT, thetaR, alpha)

global Nr Nt

L = size(thetaT, 2);

A = zeros(Nr * Nt, 4 * L); %alpha_R, alpha_I, thetar, thetat

for i = 1 : L
    thetat = thetaT(i);
    thetar = thetaR(i);
    a = alpha(i);
    tmp = exp(1j * pi * (0:Nr-1)' * sin(thetar)) * exp(1j * pi * (0:Nr-1)' * sin(thetat))';
    A(:, i) = tmp(:);
    A(:, L + i) = 1j * tmp(:);
    for m = 1 : Nr
        tmp(m) = a * 1j * pi * (m-1) * cos(thetar) *  exp(1j * pi * (m-1) * sin(thetar));
    end
    tmp = tmp * exp(1j * pi * (0:Nr-1)' * sin(thetat))';
    A(:, 2 * L + i) = tmp(:);
    for n = 1 : Nt
        tmp(n) = -a * 1j * pi * (n-1) * cos(thetat) *  exp(-1j * pi * (n-1) * sin(thetat));
    end
    tmp = exp(1j * pi * (0:Nr-1)' * sin(thetar)) * tmp.';
    A(:, 3 * L + i) = tmp(:);
end

if Nt == 1
    A = A(:, 3 * L);
end