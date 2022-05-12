

%-----------------------Harmonicas--------------------------%

%calcula o número de harmónicas
n_harm = floor((n_freqs-1)*delta_f/f0);
harmonicas = 1:n_harm;

for i=1:n_harm
 harmonicas(i) = espetro(1 + round(i * f0/delta_f));
end