
N = 10;
S = diag(randi([0 1], 1, N));
A=randn(10)+1j*randn(10);
B=randn(10)+1j*randn(10);
H1 = A*S*B;

H2 = zeros(N);
Bt = B.';
for i=1:N
  H2 = H2+S(i,i)*A(:,i)*Bt(:,i).';
end

norm(H1-H2)

R1 = log2(det(eye(N) + H1*H1'));


HH = zeros(N);
for i=1:N
  for j=1:N
    HH = HH + S(i,i)*S(j,j)*(Bt(:,i).'*conj(Bt(:,j)))*A(:,i)*A(:,j)';
  end
end
R2 = log2(det(eye(N) + HH));

norm(R1-R2)


HH = zeros(N);
for k=1:N
    HH = HH + S(k,k)*(Bt(:,k).'*conj(Bt(:,k)))*A(:,k)*A(:,k)';
end
R3 = log2(det(eye(N) + HH));

norm(R1-R3)

mu=10;
H1 = H1/norm(H1);
log2(real(det(eye(N) + 1/mu*H1*H1')))
log2(1 + 1/mu*trace(H1*H1'))
