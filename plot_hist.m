figure
hold on
axis([0 30 0 3500])
for i = 1:78
plot(HIST_RGB(i,:,3),'b')

title('Histogramm Blau alle F��e','FontSize',16)
xlabel('Farbwert','FontSize',12)
ylabel('H�ufigkeit','FontSize',12)

end