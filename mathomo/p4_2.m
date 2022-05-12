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

%-----------------------Valor Médio--------------------------%
%funciona daqui para baixo

nsamplesperiod = fs/f0;
nperiods = nsamples/nsamplesperiod;
ncompleteperiods = floor(abs(nperiods));
naverage = round(ncompleteperiods * nsamplesperiod);
media = sum(data)/naverage
