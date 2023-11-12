function makemov_im(dt, intt, intv)

v = VideoWriter('movie_test.avi');
open(v);

fre = load('fre.dat');
fim = load('fim.dat');

eff = load('e.dat');
t = eff(:,1);

N = size(fim,2)-1;
f = fre(:,2:end) + 1i*fim(:,2:end);
z = fre(:,1);

mine = min(eff(:,2));
maxe = max(eff(:,2));
 
 
for k=1:intv:N
    
    last_index_in_t = 1+intt*(k-1);
    
    FigHandle = figure();    
    plot(z, abs(f(:,k)), 'LineWidth',2);               
    title(sprintf('\\tau: %.4f' , t(last_index_in_t)));
    set(gca, 'FontSize', 16)
    xlabel('\zeta')
    ylabel('|F|')
    xlim([0 z(end)])
    
    InsetHandle = figure();
    plot(t(1:last_index_in_t), eff(1:last_index_in_t,2), 'r');               
    title('\eta_\perp');
%     xlim([0 t(1+intt*(N-1))])
    ylim([mine 0.8])
    
    [NewAxHandle, NewFigHandle]=inset(FigHandle,InsetHandle);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    saveas(NewFigHandle, sprintf('%08.4f.bmp', (k-1)*dt), 'bmp');
    
    close(FigHandle); close(InsetHandle); close(NewFigHandle);
    
%     cla    
    disp(k)
end

close(v);
close(gcf);
end

