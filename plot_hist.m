close all
figure
hold on
axis([1 30 0 3500])
for i = 1:78
plot(HIST_RGB(i,:,3),'b')

end

plot([20 20], [0 3500],'-r')
title('Histogramm Blau alle Füße','FontSize',16)
xlabel('Farbwert','FontSize',12)
ylabel('Häufigkeit','FontSize',12)




