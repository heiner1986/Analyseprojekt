clear all;
close all;
clc;

% Koerpermaﬂe
m_gr = 200;
m_gr_sd  = 20;

m_gew = 95;
m_gew_sd  = 15;

w_gr = 150;
w_gr_sd  = 20;

w_gew = 50;
w_gew_sd  = 15;



%% Testdaten
x(1:20,1) = m_gr + m_gr_sd.*randn(20,1);
x(21:40,1) = w_gr + w_gr_sd.*randn(20,1);

x(1:20,2) = m_gew + m_gew_sd.*randn(20,1);
x(21:40,2) = w_gew + w_gew_sd.*randn(20,1);

%% Training
xtraining(1:50,1) = m_gr + m_gr_sd.*randn(50,1);
xtraining(51:100,1) = w_gr + w_gr_sd.*randn(50,1);

xtraining(1:50,2) = m_gew + m_gew_sd.*randn(50,1);
xtraining(51:100,2) = w_gew + w_gew_sd.*randn(50,1);

group(1:50,1) = 1;
group(51:100,1) = 2;

[class,err,Posterior,logp,coeff] = classify(x,xtraining,group)

orggroup(1:20,1) = 1;
orggroup(21:40,1) = 2;

figure
hold on

nsubjects = length(x);
clr = [0 0 1; 1 0 0];
for n = 1:nsubjects
    plot(x(n,1),x(n,2),'o','Color',clr(class(n),:));
    plot(x(n,1),x(n,2),'x','Color',clr(orggroup(n),:));
end

    