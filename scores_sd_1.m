load("/Users/lizzy/ppp_patterns/Data/mn2_all.mat");

size(EH)
size(EWA)

% scale EWA by normalization constant
all(EWA * diag(denom) == EWA_scaled)
all(EWA(:,1)*denom(1) == EWA_scaled(:,1))
all(EWA(:,2)*denom(2) == EWA_scaled(:,2))

% scale varWA by normalization constant
% var(aX) = a^2var(x)
varWA_scaled = varWA * diag(denom.^2);

% scale EWA to sd = 1
std_denom = 1./std(EWA_scaled);

EWA_sd1 = EWA_scaled * diag(std_denom);
std(EWA_sd1)

% scale varWA to sd(EWA) = 1
% var(aX) = a^2var(x)
varWA_sd1 = varWA_scaled * diag(std_denom.^2);

% figure(1);
% histogram(varWA)
% figure(2);
% histogram(varWA_scaled)
% 
% figure(3);
% histogram(EWA)
% figure(4);
% histogram(EWA_scaled)

tiledlayout(2,3)

nexttile
histogram(EWA)
title('EWA')

% Bottom plot
nexttile
histogram(varWA_scaled)
title('EWA scaled')

% Bottom plot
nexttile
histogram(EWA_sd1)
title('EWA sd1')

nexttile
histogram(varWA)
title('varWA')

% Bottom plot
nexttile
histogram(varWA_scaled)
title('varWA scaled')

% Bottom plot
nexttile
histogram(varWA_sd1)
title('varWA sd1')

save("/Users/lizzy/ppp_patterns/Data/mn2_EWA_sd1.mat", 'EWA_sd1');
save("/Users/lizzy/ppp_patterns/Data/mn2_WA_var_sd1.mat", 'varWA_sd1');

save("/Users/lizzy/ppp_patterns/Data/mn2_Walpha.mat", 'alphaW');
save("/Users/lizzy/ppp_patterns/Data/mn2_Wbeta.mat", 'betaW');
save("/Users/lizzy/ppp_patterns/Data/mn2_Aalpha.mat", 'alphaA');
save("/Users/lizzy/ppp_patterns/Data/mn2_Abeta.mat", 'betaA');

save("/Users/lizzy/ppp_patterns/Data/mn2_EWA.mat", 'EWA_scaled');
save("/Users/lizzy/ppp_patterns/Data/mn2_EH.mat", 'EH_scaled');

save("/Users/lizzy/ppp_patterns/Data/mn2_EWA_un.mat", 'EWA');
save("/Users/lizzy/ppp_patterns/Data/mn2_WA_var_un.mat", 'varWA');

