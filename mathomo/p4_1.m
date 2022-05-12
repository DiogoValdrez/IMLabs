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

%-----------------------Valor Eficaz--------------------------%
%to get naverage use p4_2

data1 = data(1:naverage);
sef = (sum(data1.^2))/naverage;
eficaz = sqrt(sef)
