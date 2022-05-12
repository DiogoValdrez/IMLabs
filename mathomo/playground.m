fs = 10000;
nsamples = 1000;
t = 0:1/fs:1.5;
data = cos(2*pi*50*t);
%-----------------------ESTIMATIVA DA FREQUENCIA--------------------------%
%funciona daqui para baixo

%DFT
espectro_c = fft(data, nsamples)/nsamples;
%calcula resolução e número de freqs diferentes positivas
delta_f = fs/nsamples;
nfreqs = floor(nsamples/2)+1;


%espetro real
aux = rem(nsamples, 2);
%extrai o lado direito do espetro e obtém as amplitudes
espetro_lateral = abs(espectro_c(1:nfreqs));
%repõem a energia do espetro
if aux ~= 0
%nº impar de aquisições - nº de frequencias positivas e negativas igual
espetro_lateral(2:end) = 2*espetro_lateral(2:end);
else
%nº par de aquisições - nº de frequencias positivas e negativas diferentes
espetro_lateral(2:(end-1)) = 2*espetro_lateral(2:(end-1));
end

espectro = espetro_lateral;

n = 2*pi/nsamples;
%indice do maximo do especrto
[~,index_M] = max(abs(espectro_c(2:nfreqs)));
%escolhe o vizinho maior
if (abs(espectro_c(index_M-1)) > abs(espectro_c(index_M+1)))
L = index_M-1;
else
L = index_M;
end

%IpDFT
U_1 = real(espectro_c(L));
U_2 = real(espectro_c(L+1));
V_1 = imag(espectro_c(L));
V_2 = imag(espectro_c(L+1));
L=L-1;
K_Opt = ( (V_2-V_1)*sin(n * L) + (U_2-U_1)*cos(n * L) ) / (U_2-U_1);
Z_1 = V_1*(K_Opt - cos(n * L))/sin(n * L) + U_1;
Z_2 = V_2*(K_Opt - cos(n * (L+1)))/sin(n * (L+1)) + U_2;
lamda = (1/n) * acos( (Z_2*cos(n * (L+1)) - Z_1*cos(n * L)) / (Z_2 - Z_1) );
%f0 é a frequência estimada
f0 = lamda * delta_f

%-----------------------Valor Médio--------------------------%

nsamplesperiod = fs/f0;
nperiods = nsamples/nsamplesperiod;
ncompleteperiods = floor(abs(nperiods));
naverage = floor(ncompleteperiods * nsamplesperiod);
media = sum(data)/naverage

%-----------------------Valor Eficaz--------------------------%

data1 = data(1:naverage);
sef = (sum(data1.^2))/naverage;
eficaz = sqrt(sef)


%-----------------------Harmonicas--------------------------%

%calcula o número de harmónicas
n_harm = floor((nfreqs-1)*delta_f/f0);
harmonicas = 1:n_harm;

for i=1:n_harm
 harmonicas(i) = espectro(1 + round(i * f0/delta_f));
end

%-------------------Distorção HArmonica---------------------%
%harmonica fundamental
 a1 = harmonicas(1)^2;
 an = 0;

 %somatorio das restantes harmonicas
 for k = 2:n_harm
 an = an + harmonicas(k)^2;
 end

 %calculo do THD
 dist_harm = 20 *log10(sqrt( an/a1 ) )
 
 %----------------------Desfasagem--------------%
 
 dados1 = cos(2*pi*50*t);
 dados2 = sin(2*pi*50*t);
 
 espectro_c1 = fft(dados1, nsamples)/nsamples;
 espectro_c2 = fft(dados2, nsamples)/nsamples;
 
 
 f_1 = round(f0 / delta_f) + 1;
 f_2 = round(f0 / delta_f) + 1;
 
 S_1 = espectro_c1(f_1);
 S_2 = espectro_c2(f_2);
 
 ang = angle(S_1) - angle(S_2);
 ang = (180*ang)/pi;
 if  ang < -180
     ang = ang + 360;
 end
 ang
 