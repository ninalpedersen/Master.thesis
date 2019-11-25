function fastICAparallel (inputdata,outputdata,p,iter) #fastICAparallel("datamat.csv","output.csv",3,50)

  pkg load statistics
  
  datamat = csvread(inputdata);
  datamat = datamat(2:end,2:end);

  iter = 5;
  r = size(datamat,1);
  c = size(datamat,2);

  datac = center(datamat);
  datac = datac';

  #Whitening
  V = (datac * datac')/r;
  [u, s, v] = svd(V);
  D = sqrt(s)^(-1);
  K = D * u';
  K = K(1:p,:);
  datatilde = K * datac;

  winit = normrnd(0,1,[p,p]);
  m = size(datatilde,2);

  f = size(datatilde,2);
  W = winit;
  [u, s, v] = svd(W);
  W = u*s^(-1)* u' *W;
  w = W;
  it = 1;
  while it < iter
    w1 = tanh(W * datatilde) * datatilde'/f;
    xg = 1-(tanh(W * datatilde)).^2;
    w2 = diag(mean(xg,2)) * W;
    w = w1 - w2;
    [u1,s1,v1] = svd(w);
    w = u1 * s1^(-1) * u1' * w;
    W = w;
    it++;
  end

  w = W * K;
  S = w * datac;
  A = w' * inv(w * w');
  datac = datac';
  K = K';
  W = W';
  A = A';
  S = S';
  csvwrite(outputdata,S)
endfunction
