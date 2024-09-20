function harmolen_plot_distributions(vp,vo)
    vm = max(abs([vp;vo]));
    figure(2)
    subplot(3,1,1)
        hist(vp(:,1),[-1:2*(10/size(vo,1))^.5:1]*vm(1)),hold on,hist(vo(:,1),[-1:2*(10/size(vo,1))^.5:1]*vm(1)),hold off
        a=get(gca,'children');
        if numel(a)>1
            set(a(1),'FaceAlpha',0.2)
            set(a(2),'FaceAlpha',0.4)
        end
        xlabel('velocity v_x')
        ylabel('counts | #')
    subplot(3,1,2)
        hist(vp(:,2),[-1:2*(10/size(vo,1))^.5:1]*vm(2)),hold on,hist(vo(:,2),[-1:2*(10/size(vo,1))^.5:1]*vm(2)),hold off
        a=get(gca,'children');
        if numel(a)>1
            set(a(1),'FaceAlpha',0.2)
            set(a(2),'FaceAlpha',0.4)
        end
        xlabel('velocity v_y')
        ylabel('counts | #')
    subplot(3,1,3)
        hist(vp(:,3),[-1:2*(10/size(vo,1))^.5:1]*vm(3)),hold on,hist(vo(:,3),[-1:2*(10/size(vo,1))^.5:1]*vm(3)),hold off
        a=get(gca,'children');
        if numel(a)>1
            set(a(1),'FaceAlpha',0.2)
            set(a(2),'FaceAlpha',0.4)
        end
        xlabel('velocity v_z')
        ylabel('counts | #')
    drawnow
    
end
