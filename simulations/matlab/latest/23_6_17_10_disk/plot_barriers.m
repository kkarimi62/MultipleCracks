function plot_barriers()

global   r_b c_bx c_by

            centers=[c_bx' c_by'];
            viscircles(centers,r_b,'Color',"b");
            pbaspect([1 1 1]);
            saveas(gcf,'png/barriers.png')

            hold on
end