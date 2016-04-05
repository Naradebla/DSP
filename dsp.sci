factor = 1/50
freqs = csvRead("DFT.csv")
freqs(:,2) = freqs(:,2) * factor
filter = csvRead("lowPass.csv")
filtered = csvRead("filteredFreqs.csv")
filtered(:,2) = filtered(:,2) * factor
plot(freqs(:,1),freqs(:,2))
plot(filter(:,1),filter(:,2))
plot(filtered(:,1),filtered(:,2))