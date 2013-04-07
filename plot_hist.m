figure
hold on
axis([0 30 0 3500])
for i = 1:78
plot(HIST_RGB(i,:,3),'b')

title('Histogramm Blau alle Füße','FontSize',16)
xlabel('Farbwert','FontSize',12)
ylabel('Häufigkeit','FontSize',12)

end