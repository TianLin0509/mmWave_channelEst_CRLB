function CRLB = compute_CRLB(A, phi)

global Nrf Nr Vn

F = 2 * Nrf / Nr / Vn * real(A'*phi'*phi*A);
CRLB = F^(-1);