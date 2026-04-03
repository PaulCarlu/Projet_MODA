function visualize_wave(phineg,phipos,D,uL,ui,uR,model)

% Sous matrices élémentaires de D original
DiL=D(ui,uL); Di = D(ui,ui); DiR=D(ui,uR);

phi_intern = -(Di\DiL)*phineg -(Di\DiR)*phipos;

PHI_uncond = [phineg ; phi_intern ; phipos];
phi_uncond=PHI_uncond(:,10);
visualize_mode(model.mesh,phi_uncond,'animate', true)
end