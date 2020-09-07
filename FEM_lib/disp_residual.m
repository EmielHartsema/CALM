function disp_residual(fighandle,myCFD)

figure(fighandle)
semilogy(1:length(myCFD.Residual.rx), myCFD.Residual.rx)
hold on
semilogy(1:length(myCFD.Residual.ry), myCFD.Residual.ry)
hold off

legend('Residual Ux', 'Residual Uy')
title('Residual plot')
xlabel('Itteration')
ylabel('scaled residual')
end

