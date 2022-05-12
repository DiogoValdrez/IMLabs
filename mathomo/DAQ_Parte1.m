% fs = 1000;
% range = 5;
% nsamples = 3000;
% 
% %DAQ
% devlist("ni")
% d = daq("ni");
% d.Rate = fs;
% ch = addinput(d,"dev3..", "ai1..", "Voltage");
% ch.Range = [-range range];
% data = read(d, nsamples);

fs = 10000;
range = 5;
nsamples = 1000;
t = 0:1/fs:1.5;
data = cos(2*pi*50*t);

%--------------------------Estimar a IPDFT-------------------------------%
espectro_compl = fft(data, nsamples)/nsamples;
%calcula resolução e número de freqs diferentes positivas
delta_f = fs/nsamples;
nfreqs = floor(nsamples/2)+1;
n = 2*pi/nsamples;
%indice do maximo do espectro
[~,index_M] = max(abs(espectro_compl(2:nfreqs)));
%escolhe o vizinho maior
if (abs(espectro_compl(index_M-1)) > abs(espectro_compl(index_M+1)))
    L = index_M-1;
else
    L = index_M;
end
%IpDFT
U_1 = real(espectro_compl(L));
U_2 = real(espectro_compl(L+1));
V_1 = imag(espectro_compl(L));
V_2 = imag(espectro_compl(L+1));
L=L-1;
K_Opt = ( (V_2-V_1)*sin(n * L) + (U_2-U_1)*cos(n * L) ) / (U_2-U_1);
Z_1 = V_1*(K_Opt - cos(n * L))/sin(n * L) + U_1;
Z_2 = V_2*(K_Opt - cos(n * (L+1)))/sin(n * (L+1)) + U_2;
lamda = (1/n) * acos( (Z_2*cos(n * (L+1)) - Z_1*cos(n * L)) / (Z_2 - Z_1) );
%f0 é a frequência estimada
f0 = lamda * delta_f

%---------------------Calcular periodos inteiros--------------------------%
nsamplesperiod = fs/f0;
ncompleteperiods = floor(abs(nsamples/nsamplesperiod));
ncorrect = floor(ncompleteperiods * nsamplesperiod);
correctdata = data(1:ncorrect);
%--------------------------Calcular Espectro------------------------------%
espectro_compl = fft(correctdata, ncorrect)/ncorrect;
delta_f = fs/ncorrect;
nfreqs = floor(ncorrect/2)+1;
remain = rem(ncorrect, 2);
lateral_s = abs(espectro_compl(1:nfreqs));

if remain ~= 0
    lateral_s(2:end) = 2*lateral_s(2:end);
else
    lateral_s(2:(end-1)) = 2*lateral_s(2:(end-1));
end

espectro = lateral_s;
%------------------------------Valor Médio--------------------------------%
media = sum(correctdata)/ncorrect;
%-----------------------------Valor Eficaz--------------------------------%
rms = sqrt(sum(correctdata.^2)/ncorrect);
%------------------------------Harmonicas---------------------------------%
nharm = floor((nfreqs-1)*delta_f/f0);
harmonicas = 1:nharm;
for i=1:nharm
    harmonicas(i) = espectro(1+round(i*f0/delta_f));
end
%--------------------------Distorção Harmonica----------------------------%
%harmonica fundamental
a1 = harmonicas(1)^2;
an = 0;
%somatorio das restantes harmonicas
for k = 2:nharm
    an = an + harmonicas(k)^2;
end
%calculo do THD
dist_harm = 20*log10(sqrt(an/a1));
harmonicasdB = 20*log10(harmonicas);  
%----------------------------------PLOT-----------------------------------%
t = 0 : 1 : ncorrect;
datat = t / fs;
datafinal = data(1:length(t));
subplot(2,1,1);
plot(datat, datafinal);

xlabel('Tempo[s]');
ylabel('Amplitude[V]');
title(sprintf('Frequência = %.4f [Hz]; Valor Médio = %.4d [V];\n Valor eficaz = %.4f [V]; Nº de amostras = %d ;\n Fs = %.4f [Hz]; Alcance = [%.1f , %.1f ]\n \\Delta_t = %.4d [s], \\Delta_f = %.4d [Hz]', f0, media, rms, nsamples, fs, -range, range, (1/fs),delta_f )); 
legend('Sinal');
f = 0 : 1 : nfreqs-1;
dataf = f(1:length(f))*delta_f; 
espectrodB = 20*log10(espectro);

subplot(2,1,2);
plot(dataf, espectrodB); xlabel('Frequência[Hz]'); ylabel('Potência[dB]'); legend('Espetro');