% fs = 1000;
% range = 5;
% nsamples = 3000;
% naquis = 10; 
%
% %DAQ
% devlist("ni")
% d = daq("ni");
% d.Rate = fs;
% ch1 = addinput(d,"dev3..", "ai1..", "Voltage");
% ch2 = addinput(d,"dev3..", "ai1..", "Voltage");
% ch1.Range = [-range range];
% ch2.Range = [-range range];
% data = read(d, nsamples);
R = 10;

fs = 10000;
range = 5;
nsamples = 1000;
t = 0:1/fs:1.5;
data = [cos(2*pi*50*t); sin(2*pi*50*t)];

%--------------------------Espectro1-------------------------------%
espectro_compl1 = fft(data(1,:), nsamples)/nsamples;%-------------why data(1,:)?
delta_f1 = fs/nsamples;
nfreqs1 = floor(nsamples/2)+1;
n = 2*pi/nsamples;
%indice do maximo do espectro
[~,index_M] = max(abs(espectro_compl1(2:nfreqs1)));
%escolhe o vizinho maior
if (abs(espectro_compl1(index_M-1)) > abs(espectro_compl1(index_M+1)))
    L = index_M-1;
else
    L = index_M;
end
%IpDFT
U_1 = real(espectro_compl1(L));
U_2 = real(espectro_compl1(L+1));
V_1 = imag(espectro_compl1(L));
V_2 = imag(espectro_compl1(L+1));
L=L-1;
K_Opt = ( (V_2-V_1)*sin(n * L) + (U_2-U_1)*cos(n * L) ) / (U_2-U_1);
Z_1 = V_1*(K_Opt - cos(n * L))/sin(n * L) + U_1;
Z_2 = V_2*(K_Opt - cos(n * (L+1)))/sin(n * (L+1)) + U_2;
lamda = (1/n) * acos( (Z_2*cos(n * (L+1)) - Z_1*cos(n * L)) / (Z_2 - Z_1) );
%f0 é a frequência estimada
f01 = lamda * delta_f1;

nsamplesperiod = fs/f0;
ncompleteperiods = floor(abs(nsamples/nsamplesperiod));
ncorrect1 = floor(ncompleteperiods * nsamplesperiod);
correctdata1 = data(1:ncorrect1);

espectro_compl1 = fft(correctdata1, ncorrect1)/ncorrect1;
delta_f1 = fs/ncorrect1;
nfreqs = floor(ncorrect1/2)+1;
remain = rem(ncorrect1, 2);
lateral_s = abs(espectro_compl1(1:nfreqs));
if remain ~= 0
    lateral_s(2:end) = 2*lateral_s(2:end);
else
    lateral_s(2:(end-1)) = 2*lateral_s(2:(end-1));
end
espectro1 = lateral_s;
%--------------------------Espectro2-------------------------------%

espectro_compl2 = fft(data(2,:), nsamples)/nsamples;%-------------why data(1,:)?
delta_f2 = fs/nsamples;
nfreqs2 = floor(nsamples/2)+1;
n = 2*pi/nsamples;
%indice do maximo do espectro
[~,index_M] = max(abs(espectro_compl2(2:nfreqs2)));
%escolhe o vizinho maior
if (abs(espectro_compl2(index_M-1)) > abs(espectro_compl2(index_M+1)))
    L = index_M-1;
else
    L = index_M;
end
%IpDFT
U_1 = real(espectro_compl2(L));
U_2 = real(espectro_compl2(L+1));
V_1 = imag(espectro_compl2(L));
V_2 = imag(espectro_compl2(L+1));
L=L-1;
K_Opt = ( (V_2-V_1)*sin(n * L) + (U_2-U_1)*cos(n * L) ) / (U_2-U_1);
Z_1 = V_1*(K_Opt - cos(n * L))/sin(n * L) + U_1;
Z_2 = V_2*(K_Opt - cos(n * (L+1)))/sin(n * (L+1)) + U_2;
lamda = (1/n) * acos( (Z_2*cos(n * (L+1)) - Z_1*cos(n * L)) / (Z_2 - Z_1) );
%f0 é a frequência estimada
f02 = lamda * delta_f2;

nsamplesperiod = fs/f0;
ncompleteperiods = floor(abs(nsamples/nsamplesperiod));
ncorrect2 = floor(ncompleteperiods * nsamplesperiod);
correctdata2 = data(1:ncorrect2);

espectro_compl2 = fft(correctdata2, ncorrect2)/ncorrect2;
delta_f2 = fs/ncorrect2;
nfreqs = floor(ncorrect2/2)+1;
remain = rem(ncorrect2, 2);
lateral_s = abs(espectro_compl2(1:nfreqs));
if remain ~= 0
    lateral_s(2:end) = 2*lateral_s(2:end);
else
    lateral_s(2:(end-1)) = 2*lateral_s(2:(end-1));
end
espectro2 = lateral_s;


median_f0 = (f01 + f02) / 2;
%----------------------Desfasagem--------------%

espectro_c1 = fft(data(1 ,:), nsamples)/nsamples;
espectro_c2 = fft(data(2 ,:), nsamples)/nsamples;
delta_f = fs/nsamples;
f_1 = round(f01 / delta_f) + 1;
f_2 = round(f02 / delta_f) + 1;
S_1 = espectro_c1(f_1);
S_2 = espectro_c2(f_2);

ang = angle(S_1) - angle(S_2);
ang = (180*ang)/pi;
if  ang < -180
    ang = ang + 360;
end
%-----------------------Valor Eficaz--------------------------%
RMS1 = sqrt(sum(correctdata1.^2)/ncorrect1);
RMS2 = sqrt(sum(correctdata2.^2)/ncorrect2);
%-----------------Cálculo da impedância-----------------------%
mod = (RMS1/RMS2)*abs(R);
fase = ang + angle(R);
%----------------------------------PLOT-----------------------------------%
% 
% t = 0 : 1 : ncorrect1;
%    data_t = t / fs;
%    data_final1 = data(1:length(t));
%    data_final2 = data(1:length(t));
% 
%    plot(data_final1, data_final2);
% 
%    xlabel('Tempo[s]');
%    ylabel('Amplitude[V]');
%    title(sprintf(['Frequência do sinal 1 = %.4f [Hz]; Frequência do sinal 2 = %.4f [Hz]; Desfasagem : %.4f [graus]; ' ...
%        'Valor eficaz do sinal 1 = %.4f [V]; Valor eficaz do sinal 2 = %.4f [V]; Nº de amostras= %.4d; Fs = %.4f [Hz];' ...
%        'Alcance = [%.4f, %.4f]; Módulo da Impedância = %.4f; Ângulo = %.4f [rad];'], f01, f02, ang, RMS1, RMS2, nsamples, fs, -range, ...
%        range, mod, fase));
%    legend('Sinais');
if ncorrect1 >= 5
    tmax = 5* (fs/median_f0);
    id_t = 0 : 1 : tmax;
    dados_t = id_t /fs;
    dados_5per1 = data(1 ,1:length(id_t));
    dados_5per2 = data(2 ,1:length(id_t));
else
    id_t = 0 : 1 : ncorrect;
    dados_t = id_t /fs;
    dados_5per1 = data(1, 1:length(id_t));
    dados_5per2 = data(2, 1:length(id_t));
end
    plot(dados_t ,dados_5per1, dados_t ,dados_5per2);
    xlabel('Tempo[s]');
    ylabel('Amplitude[V]');
    title(sprintf('Frequência 1 = %.4f [Hz]; Frequência 2 = %.4f [Hz]; Desfasagem = %.4d [º];\n Valor eficaz 1 = %.4f [V]; Valor eficaz 2 = %.4f [V]; Nº de amostras = %d ;\n Fs = %.4f [Hz]; Alcance = [%.1f , %.1f ]\n'...
        , f01, f02, ang, RMS1, RMS2, nsamples, fs, -range, range)); 
    legend('Sinais Adquirido');