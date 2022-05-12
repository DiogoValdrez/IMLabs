% fs = 1000;
% range = 5;
% nsamples = 3000;
% naquis = 10; 
%
% %DAQ
% devlist("ni")
% d = daq("ni");
% d.Rate = fs;
% ch = addinput(d,"dev3..", "ai1..", "Voltage");
% ch.Range = [-range range];

fs = 10000;%--------------em ultimo caso tirar a estimativa do f0 como a parte 3 do diogo
range = 5;
nsamples = 15001;
naquis = 3;
t = 0:1/fs:1.5;


datafreqs = zeros(1, naquis);
data = zeros(naquis, nsamples);

for i = 1 : naquis
    data1 = 100*cos(2*pi*(1/10)*t+47);%read(d, nsamples);
    data(i,:) = data1;
    
    %Estimar a IPDFT
    espectro_compl = fft(data, nsamples)/nsamples;
    %calcula resolução e número de freqs diferentes positivas
    delta_f = fs/nsamples;
    nfreqs = floor(nsamples/2)+1;
    n = 2*pi/nsamples;
    %indice do maximo do espectro
    [~,index_M] = max(abs(espectro_compl(2:nfreqs)));
    %abs(espectro_compl(2:nfreqs))
    %escolhe o vizinho maior
    
    %tentar corrigir bug
    if index_M-1 < 1
        aux = 1;
    else
        aux = index_M-1;
    end
    
    
    if (abs(espectro_compl(aux)) > abs(espectro_compl(index_M+1)))
        L = aux;
    else
        L = index_M;%--------------------------------------+1?
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
    %frequência estimada
    datafreqs(i) = lamda * delta_f;
end

%delete(d)
%clear d

%Valor Médio de frequencia 0
mediaf = sum(datafreqs)/naquis;


%Calcular periodos inteiros
nsamplesperiod = fs/f0;
ncompleteperiods = floor(abs(nsamples/nsamplesperiod));
ncorrect = floor(ncompleteperiods * nsamplesperiod);
correctdata = data(:, 1:ncorrect);
signaldata = correctdata(1,:);
%Calcular Espectro
for i = 1:naquis
    espectro_compl = fft(correctdata(i, :), ncorrect)/ncorrect;
    delta_f_pred = fs/ncorrect;
    nfreqs_pred = floor(ncorrect/2)+1;
    remain = rem(ncorrect, 2);
    lateral_s = abs(espectro_compl(1:nfreqs_pred));

    if remain ~= 0
        lateral_s(2:end) = 2*lateral_s(2:end);
    else
        lateral_s(2:(end-1)) = 2*lateral_s(2:(end-1));
    end
    espectro = lateral_s;  
    if i == 1
        nfreqs = nfreqs_pred;
        delta_f = delta_f_pred;
        espectros = zeros(naquis, nfreqs);
    end
    espectros(i, :) = espectro;
end
espectro = sum(espectros)/naquis;
espectropow = (espectro.^2)/2;
%Valor Médio
media = sum(signaldata)/ncorrect;
%Valor Eficaz
rms = sqrt(sum(signaldata.^2)/ncorrect);
%Ruido Eficaz
aux = round(mediaf/delta_f)+1;

ruido = sum(espectropow([2:aux-1, aux+1:length(espectropow)]));%ruido = sum(espectropow([2:aux-1, aux+1:length(espectropow)]));
ruidoEF = sqrt(ruido);
%ENOB e SINAD
harm1 = round(mediaf/deltaf)+1;
pow1 = 10*log10(espectropow(harm1));
pow2 = 10 * log10(range);
if pow1 > pow2
    pow0 = pow1;
else
    pow0 = pow2;
end
ruido = sum(espectropow([2:harm1-1, harm1+1:length(espectropow)]));
ruidoPOW = 10*log10(ruido);
sinaddB = pow0 - ruidoPOW;
enob = (sinaddB-1.76)/6.02;




t = 0 : 1 : corrected_n;
data_t = t/fs;
aux = data(1,:);
data_f = aux(1:length(t));
subplot(2,1,1);

plot(data_t, data_f);

xlabel('Tempo[s]');
ylabel('Amplitude[V]');
title(sprintf(['Frequência = %.4f [Hz]; Valor Médio = %.4d [V]; Valor Eficaz = %.4f [V]; Ruído Eficaz = %.4f [dB]; Nº eficaz de bits = %.4f; ' ...
    'Nº de amostras = %d; Fs = %d [Hz]; Alcance = [%.2f, %.2f]'], mediaf, media, rms, 20*log10(ruidoEF), enob, nsamples, fs, ...
    -range, range));
legend('Sinal');

f = 0 : 1 : nfreqs-1;
data_f = f(1:length(f) ) * delta_f; 
espectrodB = 20*log10(espectro);

subplot(2,1,2);
plot(data_f, espectrodB); xlabel('Frequência[Hz]'); ylabel('Potência[dB]'); legend('Espetro das aquisições');