function [ U_next, Bi_next, V_next ] = Bidiag_Francis_Step_Update_U_V(Ui, Bi, Vi)
  A = Bi;
  
  function [ G ] = Givens(tau11, tau12, tauMM)
      x = [ tau11 - tauMM; tau12 ];
      gamma = x(1) / norm(x);
      sigma = x(2) / norm(x);
      G = [ gamma -sigma; sigma gamma ];
  end

  U = Ui;
  V = Vi;

  [ m, ~ ] = size(A);
  
  % Compute first Givens rotation from T.
  G = Givens(A(1,1)^2, A(1,2) * A(1,1), A(m-1,m)^2 + A(m,m)^2);

  % Apply G to A from right
  coli = A(:,1);
  A(:,1) = coli * G(1,1) + A(:,2) * G(2,1);
  A(:,2) = coli * G(1,2) + A(:,2) * G(2,2);

  %% NEW - Apply Givens rotation to V
  % Apply G to V from the right
  coli = V(:,1);
  V(:,1) = coli * G(1,1) + V(:,2) * G(2,1);
  V(:,2) = coli * G(1,2) + V(:,2) * G(2,2);  

  % Chase the bulge
  i = 2;
  j = 1;
  while i < m
    if i > j
      G = Givens(A(i-1,j), A(i,j), 0);

      % Apply to A from the left
      rowi = A( i-1,: );
      A(i-1,:) = G(1,1) * rowi + G(2,1) * A(i,:);
      A(i,:) = G(1,2) * rowi + G(2,2) * A(i,:);
      
      %% NEW - Apply Givens rotation to U
      % Apply to U from the right
      coli = U( :,j );
      U(:,j) = coli * G(1,1) + U(:,j+1) * G(2,1);
      U(:,j+1) = coli * G(1,2) + U(:,j+1) * G(2,2);
    else
      G = Givens(A(i,j-1), A(i, j), 0);
      
      % Apply to A from the right
      coli = A(:,j-1);
      A(:,j-1) = coli * G(1,1) + A(:,j) * G(2,1);
      A(:,j) = coli * G(1,2) + A(:,j) * G(2,2);

      %% NEW - Apply Givens rotation to V
      % Apply to V from the right
      coli = V(:,j-1);
      V(:,j-1) = coli * G(1,1) + V(:,j) * G(2,1);
      V(:,j) = coli * G(1,2) + V(:,j) * G(2,2);  
    end
    tmp = i;
    i = j;
    j = tmp + 1;
  end

  % Remove the bulge on the last row
  G = Givens(A(m-1,m-1), A(m, m-1), 0);
  
  % Apply to A from the left
  rowi = A(m-1,:);
  A(m-1,:) = G(1,1) * rowi + G(2,1) * A(m,:);
  A(m,:) = G(1,2) * rowi + G(2,2) * A(m,:);

  %% NEW - Apply Givens rotation to U
  % Apply to U from the right
  coli = U(:,m-1);
  U(:,m-1) = coli * G(1,1) + U(:,m) * G(2,1);
  U(:,m) = coli * G(1,2) + U(:,m) * G(2,2);
  
  Bi_next = A;
  U_next = U;
  V_next = V;
    
end