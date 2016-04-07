windw4 = scf(4);
windw0 = scf(0);
factor = 1/50
freqs = csvRead("complex.csv")
freqs(:,2) = freqs(:,2) * factor
filter = csvRead("lowPass.csv")
filtered = csvRead("filteredFreqs.csv")
filtered(:,2) = filtered(:,2) * factor
plot(freqs(:,1),freqs(:,2),'r')
plot2d(filter(:,1),filter(:,2))
plot(filtered(:,1),filtered(:,2),'b')
scf(windw4)
in = csvRead("test.csv")
plot(in(:,1),in(:,2),'r')
