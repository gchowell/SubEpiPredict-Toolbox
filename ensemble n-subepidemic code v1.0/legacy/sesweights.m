function weights=sesweights(n,alpha1)

weights= alpha1*(1-alpha1).^((1:n)-1)';

weights=weights(end:-1:1);

