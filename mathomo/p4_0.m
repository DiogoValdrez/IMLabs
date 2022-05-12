fs = ;
range = ;
nsamples = ;

devlist("ni")
d = daq("ni");
d.Rate = fs;
ch = addinput(d,"dev3..", "ai1..", "Voltage");
ch.Range = [-range range];
%d.Channels
%d.Channels.TerminalConfig = "SingleEnded"; % é preciso isto?
data = read(d, nsamples);%, "OutputFormat..?", "Matrix");
%plot(data)

%-----------------------ESTIMATIVA DA FREQUENCIA--------------------------%
%funciona daqui para baixo

%DFT
espectro = fft(data, nsamples)/nsamples;
%calcula resolução e número de freqs diferentes positivas
delta_f = fs/nsamples;
nfreqs = floor(nsamples/2)+1;

n = 2*pi/nsamples;
%indice do maximo do especrto
[~,index_M] = max(abs(espectro(2:nfreqs)));
%escolhe o vizinho maior
if (abs(espectro(index_M-1)) > abs(espectro(index_M+1)))
L = index_M-1;
else
L = index_M;
end

%IpDFT
U_1 = real(espectro(L));
U_2 = real(espectro(L+1));
V_1 = imag(espectro(L));
V_2 = imag(espectro(L+1));
L=L-1;
K_Opt = ( (V_2-V_1)*sin(n * L) + (U_2-U_1)*cos(n * L) ) / (U_2-U_1);
Z_1 = V_1*(K_Opt - cos(n * L))/sin(n * L) + U_1;
Z_2 = V_2*(K_Opt - cos(n * (L+1)))/sin(n * (L+1)) + U_2;
lamda = (1/n) * acos( (Z_2*cos(n * (L+1)) - Z_1*cos(n * L)) / (Z_2 - Z_1) );
%f0 é a frequência estimada
f0 = lamda * delta_f




