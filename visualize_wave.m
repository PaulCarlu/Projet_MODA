function visualize_wave(wavebasis)

% Superposition des racines avec les courbes de dispersion
figure(1);hold on;          % affichage de (W,K)
plot(wn.^2,kpn.^2,'ob');
plot(wn.^2,kmn.^2,'or');
plot(wc.^2,0,'sm');

figure(Aff+1);                  % affichage de (w,k)
subplot(2,1,1);hold on;     %-- partie reelle
plot(w,real(kp),'b','LineWidth',2);
plot(w,real(km),'r','LineWidth',2);
plot(wn,real(kpn),'ob','MarkerSize',8);
plot(wn,real(kmn),'or','MarkerSize',8);
plot(wc,0,'sm','MarkerSize',8,'MarkerFaceColor','m');

end

