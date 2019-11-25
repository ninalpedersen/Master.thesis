function fastICAdeflation (inputdata,outputdata,p,iter) #fastICAdeflation("datamat.csv","output.csv",3,50)

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
  W = zeros(p,p);

  for i = 1:p
    w = winit(i,:)';
    if i > 1
      tempw = w;
      tempw(:) = 0;
      for j = 1:(i-1)
        tempw = tempw + sum(w .* W(j,:)) .* W(j,:);
      end
      w = w - tempw;
    end
    w = w/sqrt(sum(w.^2));
    it = 1;
    while it < iter
      xg = datatilde .* repmat(tanh(w' * datatilde),[p,1]);
      w1 = mean(xg') - (mean(1 - (tanh(w' * datatilde)).^2) .* w);
      w1 = diag(w1);
      if i > 1
        tempw = w1;
        tempw(:) = 0;
        for k = 1:(i-1)
          tempw = tempw + sum(w1 .* W(k,:)) .* W(k,:);
        end
        w1 = w1 - tempw;
      end
      w1 = w1/sqrt(sum(w1.^2));
      w = w1;
      it++;
    end
    W(i,:) = w;
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
