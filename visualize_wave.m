function visualize_wave(wavebasis)
% Superposition des racines avec les courbes de dispersion
figure(Aff);hold on;          % affichage de (W,K)
plot(wn.^2,kpn.^2,'ob','MarkerSize',8);
plot(wn.^2,kmn.^2,'or','MarkerSize',8);
plot(wc.^2,0,'sm','MarkerSize',8,'MarkerFaceColor','m');
axis tight;box on;

figure(Aff+1);                  % affichage de (w,k)
subplot(2,1,1);hold on;     %-- partie reelle
plot(w,real(kp),'b','LineWidth',2);
plot(w,real(km),'r','LineWidth',2);
plot(wn,real(kpn),'ob','MarkerSize',8);
plot(wn,real(kmn),'or','MarkerSize',8);
plot(wc,0,'sm','MarkerSize',8,'MarkerFaceColor','m');
hh=legend('$k^{+}$','$k^{-}$','$k^{+}_{n}$','$k^{-}_{n}$','$\omega_c$');
ylabel('$Re(k)$','FontSize',16,'Interpreter','latex');
axis tight;box on;
set(gca,'FontSize',20)
set(hh,'FontSize',20,'Interpreter','latex','Location','northwest')
subplot(2,1,2);hold on;     %-- partie imaginaire
plot(w,imag(kp),'b','LineWidth',2);
plot(w,imag(km),'r','LineWidth',2);
plot(wn,imag(kpn),'ob','MarkerSize',8);
plot(wn,imag(kmn),'or','MarkerSize',8);
plot(wc,0,'sm','MarkerSize',8,'MarkerFaceColor','m');
hh=legend('$k^{-}$','$k^{-}_{n}$','$\omega_c$');
xlabel('$\omega$','FontSize',16,'Interpreter','latex');
ylabel('$Im(k)$','FontSize',16,'Interpreter','latex');
axis tight;box on;
set(gca,'FontSize',20)
set(hh,'FontSize',20,'Interpreter','latex','Location','northeast')
end

