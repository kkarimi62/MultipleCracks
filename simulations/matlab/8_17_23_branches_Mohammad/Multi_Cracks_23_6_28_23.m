close all; clc; clear all;
format short e

global  n_c C S_xx S_yy S_xy Sigma_xx Sigma_yy Sigma_xy n_grid u_x u_y U V ....
         X Y force_ratio eps sig_e  t position stress stress1 stress2 x_cod_up ....
        y_cod_up x_cod_down y_cod_down max_iterations recomb_L iter_max ...
        n_branch K_IC G Eeff Length_unit sig_amp n_steps 

    
control_data()
counter = 0;
tic
%%%%%%%%%%%%%%%%%%%%%  Establish initial distribution of dislocations %%%%
for i_q=1:n_c 
   for k_q=1:n_branch

        switch(C(i_q).B(k_q).mode)
            
            case 'pure'
while C(i_q).B(k_q).K >= K_IC
counter = counter+1;
force_ratio=ones(n_c,n_branch);
if counter>1
update(i_q,k_q);
end
% close
    iter=0;
    max_iterations=iter_max*C(i_q).B(k_q).n_d;
        while force_ratio(i_q,k_q) >eps && iter<max_iterations
        iter=iter+1;


                for j_q=1:2*C(i_q).B(k_q).n_d
                    C(i_q).B(k_q).D(j_q).sigma=zeros(3,3);
                    % C(i_q).B(k_q).D(j_q).sigma_disk=disk_field_pressure(C(i_q).B(k_q).D(j_q).R_p);
                    C(i_q).B(k_q).D(j_q).sigma_uniaxial=sig_e;            %uniaxial
                    C(i_q).B(k_q).D(j_q).disp=zeros(1,3);
                    cod_up=C(i_q).B(k_q).D(j_q).R_p;
                    cod_down=cod_up;

                   for i_p=1:n_c 
                     for k_p=1:n_branch   
                        for j_p=1:2*C(i_p).B(k_p).n_d 
                            switch(C(i_p).B(k_p).mode)
                                case 'pure'
                           R_pq=C(i_q).B(k_q).D(j_q).R_p-C(i_p).B(k_p).D(j_p).R_p;
                           QQ=rotation(i_p,j_p,k_p,0);
                           TT=rotation(i_p,j_p,k_p,1);
                           r_pq=norm(R_pq);
                           annihilation=C(i_q).B(k_q).D(j_q).orientation+C(i_p).B(k_p).D(j_p).orientation;
                               if (r_pq<C(i_q).B(k_q).d/recomb_L || r_pq<C(i_p).B(k_p).d/recomb_L) && annihilation==0
                                  C(i_q).B(k_q).D(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D(j_p).b=[0 0 0]; 
                               elseif r_pq==0 
                                    'self force=0';
                               else
                                [stress,displacement]=dislocation_field(R_pq*QQ);
                                stress=norm(C(i_q).B(k_q).D(j_q).b)*TT*stress*TT';
                                C(i_q).B(k_q).D(j_q).sigma= C(i_q).B(k_q).D(j_q).sigma+stress;
                               end
                           C(i_q).B(k_q).D(j_q).recomb=ceil(norm(C(i_q).B(k_q).D(j_q).b));
                           C(i_p).B(k_p).D(j_p).recomb=ceil(norm(C(i_p).B(k_p).D(j_p).b));
                            
                        case 'mixed'
                           R_pq1=C(i_q).B(k_q).D(j_q).R_p-C(i_p).B(k_p).D1(j_p).R_p;
                           R_pq2=C(i_q).B(k_q).D(j_q).R_p-C(i_p).B(k_p).D2(j_p).R_p;
                           QQ=rotation(i_p,j_p,k_p,0);
                           r_pq1=norm(R_pq1);
                           r_pq2=norm(R_pq2);
                           annihilation1=C(i_q).B(k_q).D(j_q).orientation+C(i_p).B(k_p).D1(j_p).orientation;
                               if (r_pq1<C(i_q).B(k_q).d/recomb_L || r_pq1<C(i_p).B(k_p).d/recomb_L) && annihilation1==0
                                  C(i_q).B(k_q).D(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D1(j_p).b=[0 0 0]; 
                               elseif r_pq1==0 
                                    'self force=0';
                               else
                                [stress1,displacement1]=dislocation_field(R_pq1*QQ);
                                stress1=norm(C(i_q).B(k_q).D(j_q).b)*QQ*stress1*QQ';
                                C(i_q).B(k_q).D(j_q).sigma= C(i_q).B(k_q).D(j_q).sigma+stress1;
                               end
                           C(i_q).B(k_q).D(j_q).recomb1=ceil(norm(C(i_q).B(k_q).D(j_q).b));
                           C(i_p).B(k_p).D1(j_p).recomb=ceil(norm(C(i_p).B(k_p).D1(j_p).b));

                            annihilation2=C(i_q).B(k_q).D(j_q).orientation+C(i_p).B(k_p).D2(j_p).orientation;
                               if (r_pq2<C(i_q).B(k_q).d/recomb_L || r_pq2<C(i_p).B(k_p).d/recomb_L) && annihilation2==0
                                  C(i_q).B(k_q).D(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D2(j_p).b=[0 0 0]; 
                               elseif r_pq2==0 
                                    'self force=0';
                               else
                                [stress2,displacement2]=dislocation_field(R_pq2*QQ);
                                stress2=norm(C(i_q).B(k_q).D(j_q).b)*QQ'*stress2*QQ;
                                C(i_q).B(k_q).D(j_q).sigma= C(i_q).B(k_q).D(j_q).sigma+stress2;
                               end
                           C(i_q).B(k_q).D(j_q).recomb2=ceil(norm(C(i_q).B(k_q).D(j_q).b));
                           C(i_p).B(k_p).D2(j_p).recomb=ceil(norm(C(i_p).B(k_p).D2(j_p).b));
                            end

                        end
                      end  
                    end
        [t(i_q,j_q,k_q),position(i_q,j_q,k_q)]=get_motion(i_q,j_q,k_q,0);
        [x_cod_up(i_q,j_q,k_q), y_cod_up(i_q,j_q,k_q),x_cod_down(i_q,j_q,k_q), ...
            y_cod_down(i_q,j_q,k_q)]=COD(i_q,j_q,k_q,0);
                end
        max_climb=max(abs(t(i_q,2:(2*C(i_q).B(k_q).n_d-1),k_q)));
                 if C(i_q).B(k_q).n_d==1
                 force_ratio(i_q,k_q)=eps;
                 else
        force_ratio(i_q,k_q)=max_climb/max(abs(t(i_q,:,k_q)),[],"all");

                 end
            convergence = force_ratio(i_q,k_q)
            plot_forces(i_q,k_q,0);
            
        end

        C(i_q).B(k_q).K = [sqrt(abs(t(i_q,1,k_q))*Length_unit*G*Eeff),sqrt(abs(t(i_q,2*C(i_q).B(k_q).n_d,k_q))*Length_unit*G*Eeff)];
        KI_th = sig_e(2,2)*G*sqrt(pi*C(i_q).B(k_q).a*Length_unit)

        t_max = max(t(i_q,:,k_q),[],"all");

        if C(i_q).B(k_q).K(1) > K_IC && C(i_q).B(k_q).K(2) > K_IC
            
            delr =  (C(i_q).B(k_q).d*t(i_q,1,k_q))/(n_steps*max(sig_amp,t_max));
            C(i_q).B(k_q).D(1).R_p = C(i_q).B(k_q).D(1).R_p + 20*delr*C(i_q).B(k_q).e_p;
            delr = C(i_q).B(k_q).d*t(i_q,2*C(i_q).B(k_q).n_d,k_q)/(n_steps*max(sig_amp,t_max));
            C(i_q).B(k_q).D(2*C(i_q).B(k_q).n_d).R_p = C(i_q).B(k_q).D(2*C(i_q).B(k_q).n_d).R_p +20*delr*C(i_q).B(k_q).e_p;
        end
        figure(500)
        plot(counter,C(i_q).B(k_q).a,'o','LineWidth',2)
        hold on

        figure(400)
        plot(counter,C(i_q).B(k_q).K(1)*10^-6,'o','LineWidth',2)
        hold on

end
        
                for j=2:2*C(i_q).B(k_q).n_d
                    if norm( C(i_q).B(k_q).D(j).b )==0
                        x_cod_up(i_q,j,k_q)=x_cod_up(i_q,j-1,k_q);
                        y_cod_up(i_q,j,k_q)=y_cod_up(i_q,j-1,k_q);
                        x_cod_down(i_q,j,k_q)=x_cod_down(i_q,j-1,k_q);
                        y_cod_down(i_q,j,k_q)=y_cod_down(i_q,j-1,k_q);
                    end          
                end

            case 'mixed'
                
while C(i_q).B(k_q).K >= K_IC
counter = counter+1;
force_ratio=ones(n_c,n_branch,2);
close
    iter=0;
    max_iterations=iter_max*C(i_q).B(k_q).n_d;
        while force_ratio(i_q,k_q,1) >eps && iter<max_iterations
        iter=iter+1;


                for j_q=1:2*C(i_q).B(k_q).n_d
                    C(i_q).B(k_q).D1(j_q).sigma=zeros(3,3);
                    % C(i_q).B(k_q).D1(j_q).sigma_disk=disk_field_pressure(C(i_q).B(k_q).D1(j_q).R_p);
                    C(i_q).B(k_q).D1(j_q).sigma_uniaxial=sig_e;            %uniaxial
                    C(i_q).B(k_q).D1(j_q).disp=zeros(1,3);

                    cod_up1=C(i_q).B(k_q).D1(j_q).R_p;
                    cod_down1=cod_up1;


                   for i_p=1:n_c 
                     for k_p=1:n_branch   
                        for j_p=1:2*C(i_p).B(k_p).n_d 

                            switch(C(i_p).B(k_p).mode)
                                case 'pure'
                           R_pq1=C(i_q).B(k_q).D1(j_q).R_p-C(i_p).B(k_p).D(j_p).R_p;

                           QQ=rotation(i_p,j_p,k_p,0);
                           r_pq1=norm(R_pq1);

                           annihilation1=C(i_q).B(k_q).D1(j_q).orientation+C(i_p).B(k_p).D(j_p).orientation;
                               if (r_pq1<C(i_q).B(k_q).d/recomb_L || r_pq1<C(i_p).B(k_p).d/recomb_L) && annihilation1==0
                                  C(i_q).B(k_q).D1(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D(j_p).b=[0 0 0]; 
                               elseif r_pq1==0 
                                    'self force=0';
                               else
                                [stress1,displacement1]=dislocation_field(R_pq1*QQ);
                                stress1=norm(C(i_q).B(k_q).D1(j_q).b)*QQ*stress1*QQ';
                                C(i_q).B(k_q).D1(j_q).sigma= C(i_q).B(k_q).D1(j_q).sigma+stress1;
                               end
                           C(i_q).B(k_q).D1(j_q).recomb=ceil(norm(C(i_q).B(k_q).D1(j_q).b));
                           C(i_p).B(k_p).D(j_p).recomb=ceil(norm(C(i_p).B(k_p).D(j_p).b));

                        case 'mixed'
                           R_pq11=C(i_q).B(k_q).D1(j_q).R_p-C(i_p).B(k_p).D1(j_p).R_p;
                           % R_pq12=C(i_q).B(k_q).D1(j_q).R_p-C(i_p).B(k_p).D2(j_p).R_p;
                           % R_pq21=C(i_q).B(k_q).D2(j_q).R_p-C(i_p).B(k_p).D1(j_p).R_p;
                           
                           QQ=rotation(i_p,j_p,k_p,0);
                           TT =rotation(i_p,j_p,k_q,1);
                           r_pq11=norm(R_pq11);
                           % r_pq12=norm(R_pq12);
                           % r_pq21=norm(R_pq21);

                           annihilation11=C(i_q).B(k_q).D1(j_q).orientation+C(i_p).B(k_p).D1(j_p).orientation;
                               if (r_pq11<C(i_q).B(k_q).d/recomb_L || r_pq11<C(i_p).B(k_p).d/recomb_L) && annihilation11==0
                                  C(i_q).B(k_q).D1(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D1(j_p).b=[0 0 0]; 
                               elseif r_pq11==0 
                                    'self force=0';
                               else
                                [stress11,displacement11]=dislocation_field(R_pq11*QQ);
                                stress11=norm(C(i_q).B(k_q).D1(j_q).b)*TT*stress11*TT';
                                C(i_q).B(k_q).D1(j_q).sigma= C(i_q).B(k_q).D1(j_q).sigma+stress11;
                               end
                           C(i_q).B(k_q).D1(j_q).recomb=ceil(norm(C(i_q).B(k_q).D1(j_q).b));
                           C(i_p).B(k_p).D1(j_p).recomb=ceil(norm(C(i_p).B(k_p).D1(j_p).b));

         
                           %  annihilation12=C(i_q).B(k_q).D1(j_q).orientation+C(i_p).B(k_p).D2(j_p).orientation;
                           %     if (r_pq12<C(i_q).B(k_q).d/recomb_L || r_pq12<C(i_p).B(k_p).d/recomb_L) && annihilation12==0
                           %        C(i_q).B(k_q).D1(j_q).b=[0 0 0];
                           %        C(i_p).B(k_p).D2(j_p).b=[0 0 0]; 
                           %     elseif r_pq12==0 
                           %          'self force=0';
                           %     else
                           %      [stress12,displacement12]=dislocation_field(R_pq12*QQ);
                           %      stress12=norm(C(i_q).B(k_q).D1(j_q).b)*QQ*stress12*QQ;
                           %      C(i_q).B(k_q).D1(j_q).sigma= C(i_q).B(k_q).D1(j_q).sigma+stress12;
                           %     end
                           % C(i_q).B(k_q).D1(j_q).recomb=ceil(norm(C(i_q).B(k_q).D1(j_q).b));
                           % C(i_p).B(k_p).D2(j_p).recomb=ceil(norm(C(i_p).B(k_p).D2(j_p).b));
                           % 
                           % 
                           %   annihilation21=C(i_q).B(k_q).D2(j_q).orientation+C(i_p).B(k_p).D1(j_p).orientation;
                           %     if (r_pq21<C(i_q).B(k_q).d/recomb_L || r_pq21<C(i_p).B(k_p).d/recomb_L) && annihilation21==0
                           %        C(i_q).B(k_q).D2(j_q).b=[0 0 0];
                           %        C(i_p).B(k_p).D1(j_p).b=[0 0 0]; 
                           %     elseif r_pq21==0 
                           %          'self force=0';
                           %     else
                           %      [stress21,displacement21]=dislocation_field(R_pq21*QQ);
                           %      stress21=norm(C(i_q).B(k_q).D2(j_q).b)*QQ*stress21*QQ;
                           %      C(i_q).B(k_q).D2(j_q).sigma= C(i_q).B(k_q).D2(j_q).sigma+stress21;
                           %     end
                           % C(i_q).B(k_q).D2(j_q).recomb=ceil(norm(C(i_q).B(k_q).D2(j_q).b));
                           % C(i_p).B(k_p).D1(j_p).recomb=ceil(norm(C(i_p).B(k_p).D1(j_p).b)); 

                end
                     end
                          end
                                 end
        [t(i_q,j_q,k_q,1),position(i_q,j_q,k_q,1)]=get_motion(i_q,j_q,k_q,1);
        [x_cod_up(i_q,j_q,k_q,1), y_cod_up(i_q,j_q,k_q,1),x_cod_down(i_q,j_q,k_q,1), ...
            y_cod_down(i_q,j_q,k_q,1)]=COD(i_q,j_q,k_q,1);

                end
        max_climb1=max(abs(t(i_q,2:(2*C(i_q).B(k_q).n_d-1),k_q,1)));

                 if C(i_q).B(k_q).n_d==1
                 force_ratio(i_q,k_q,:)=eps;
                 else
        force_ratio(i_q,k_q,1)=max_climb1/max(abs(t(i_q,:,k_q,1)),[],"all");

                 end
            convergence1 = force_ratio(i_q,k_q,1)
            plot_forces(i_q,k_q,1);
            
        end
iter = 0;
        while force_ratio(i_q,k_q,2) >eps && iter<max_iterations
iter = iter + 1

for j_q=1:2*C(i_q).B(k_q).n_d

                    C(i_q).B(k_q).D2(j_q).sigma=zeros(3,3);
                    % C(i_q).B(k_q).D2(j_q).sigma_disk=disk_field_pressure(C(i_q).B(k_q).D2(j_q).R_p);
                    C(i_q).B(k_q).D2(j_q).sigma_uniaxial=sig_e;           %uniaxial
                    C(i_q).B(k_q).D2(j_q).disp=zeros(1,3);

                    cod_up2=C(i_q).B(k_q).D2(j_q).R_p;
                    cod_down2=cod_up2;

                   for i_p=1:n_c 
                     for k_p=1:n_branch   
                        for j_p=1:2*C(i_p).B(k_p).n_d 

                            switch(C(i_p).B(k_p).mode)
                                case 'pure'
                           R_pq2=C(i_q).B(k_q).D2(j_q).R_p-C(i_p).B(k_p).D(j_p).R_p;

                           QQ=rotation(i_p,j_p,k_p,0);
                           r_pq2=norm(R_pq2);


                         annihilation2=C(i_q).B(k_q).D2(j_q).orientation+C(i_p).B(k_p).D(j_p).orientation;
                               if (r_pq2<C(i_q).B(k_q).d/recomb_L || r_pq2<C(i_p).B(k_p).d/recomb_L) && annihilation2==0
                                  C(i_q).B(k_q).D2(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D(j_p).b=[0 0 0]; 
                               elseif r_pq2==0 
                                    'self force=0';
                               else
                                [stress2,displacement2]=dislocation_field(R_pq2*QQ);
                                stress2=norm(C(i_q).B(k_q).D2(j_q).b)*QQ*stress2*QQ';
                                C(i_q).B(k_q).D2(j_q).sigma= C(i_q).B(k_q).D2(j_q).sigma+stress2;
                               end
                           C(i_q).B(k_q).D2(j_q).recomb=ceil(norm(C(i_q).B(k_q).D2(j_q).b));
                           C(i_p).B(k_p).D(j_p).recomb=ceil(norm(C(i_p).B(k_p).D(j_p).b));
                            
                        case 'mixed'
                           % R_pq12=C(i_q).B(k_q).D1(j_q).R_p-C(i_p).B(k_p).D2(j_p).R_p;
                           % R_pq21=C(i_q).B(k_q).D2(j_q).R_p-C(i_p).B(k_p).D1(j_p).R_p;
                           R_pq22=C(i_q).B(k_q).D2(j_q).R_p-C(i_p).B(k_p).D2(j_p).R_p;
                           
                           QQ=rotation(i_p,j_p,k_p,0);
                           TT =rotation(i_p,j_p,k_q,1);
                           % r_pq12=norm(R_pq12);
                           % r_pq21=norm(R_pq21);
                           r_pq22=norm(R_pq22);

                         
                            annihilation22=C(i_q).B(k_q).D2(j_q).orientation+C(i_p).B(k_p).D2(j_p).orientation;
                               if (r_pq22<C(i_q).B(k_q).d/recomb_L || r_pq22<C(i_p).B(k_p).d/recomb_L) && annihilation22==0
                                  C(i_q).B(k_q).D2(j_q).b=[0 0 0];
                                  C(i_p).B(k_p).D2(j_p).b=[0 0 0]; 
                               elseif r_pq22==0 
                                    'self force=0';
                               else
                                [stress22,displacement22]=dislocation_field(R_pq22*QQ);
                                stress22=norm(C(i_q).B(k_q).D2(j_q).b)*TT*stress22*TT';
                                C(i_q).B(k_q).D2(j_q).sigma= C(i_q).B(k_q).D2(j_q).sigma+stress22;
                               end
                           C(i_q).B(k_q).D2(j_q).recomb=ceil(norm(C(i_q).B(k_q).D2(j_q).b));
                           C(i_p).B(k_p).D2(j_p).recomb=ceil(norm(C(i_p).B(k_p).D2(j_p).b));

                           %  annihilation12=C(i_q).B(k_q).D1(j_q).orientation+C(i_p).B(k_p).D2(j_p).orientation;
                           %     if (r_pq12<C(i_q).B(k_q).d/recomb_L || r_pq12<C(i_p).B(k_p).d/recomb_L) && annihilation12==0
                           %        C(i_q).B(k_q).D1(j_q).b=[0 0 0];
                           %        C(i_p).B(k_p).D2(j_p).b=[0 0 0]; 
                           %     elseif r_pq12==0 
                           %          'self force=0';
                           %     else
                           %      [stress12,displacement12]=dislocation_field(R_pq12*QQ);
                           %      stress12=norm(C(i_q).B(k_q).D1(j_q).b)*QQ*stress12*QQ;
                           %      C(i_q).B(k_q).D1(j_q).sigma= C(i_q).B(k_q).D1(j_q).sigma+stress12;
                           %     end
                           % C(i_q).B(k_q).D1(j_q).recomb=ceil(norm(C(i_q).B(k_q).D1(j_q).b));
                           % C(i_p).B(k_p).D2(j_p).recomb=ceil(norm(C(i_p).B(k_p).D2(j_p).b));
                           % 
                           % 
                           %   annihilation21=C(i_q).B(k_q).D2(j_q).orientation+C(i_p).B(k_p).D1(j_p).orientation;
                           %     if (r_pq21<C(i_q).B(k_q).d/recomb_L || r_pq21<C(i_p).B(k_p).d/recomb_L) && annihilation21==0
                           %        C(i_q).B(k_q).D2(j_q).b=[0 0 0];
                           %        C(i_p).B(k_p).D1(j_p).b=[0 0 0]; 
                           %     elseif r_pq21==0 
                           %          'self force=0';
                           %     else
                           %      [stress21,displacement21]=dislocation_field(R_pq21*QQ);
                           %      stress21=norm(C(i_q).B(k_q).D2(j_q).b)*QQ*stress21*QQ;
                           %      C(i_q).B(k_q).D2(j_q).sigma= C(i_q).B(k_q).D2(j_q).sigma+stress21;
                           %     end
                           % C(i_q).B(k_q).D2(j_q).recomb=ceil(norm(C(i_q).B(k_q).D2(j_q).b));
                           % C(i_p).B(k_p).D1(j_p).recomb=ceil(norm(C(i_p).B(k_p).D1(j_p).b)); 

                end
                     end
                          end
                                 end
        [t(i_q,j_q,k_q,2),position(i_q,j_q,k_q,2)]=get_motion(i_q,j_q,k_q,2);

        [x_cod_up(i_q,j_q,k_q,2), y_cod_up(i_q,j_q,k_q,2),x_cod_down(i_q,j_q,k_q,2), ...
            y_cod_down(i_q,j_q,k_q,2)]=COD(i_q,j_q,k_q,2);
                end
        max_climb2=max(abs(t(i_q,2:(2*C(i_q).B(k_q).n_d-1),k_q,2)));

                 if C(i_q).B(k_q).n_d==1
                 force_ratio(i_q,k_q,:)=eps;
                 else
        force_ratio(i_q,k_q,2)=max_climb2/max(abs(t(i_q,:,k_q,2)),[],"all");

                 end
            convergence2 = force_ratio(i_q,k_q,2)
            plot_forces(i_q,k_q,2);
            
        end
        C(i_q).B(k_q).KI = [sqrt(abs(t(i_q,1,k_q,1))*Length_unit*G*Eeff),sqrt(abs(t(i_q,2*C(i_q).B(k_q).n_d,k_q,1))*Length_unit*G*Eeff)]
        C(i_q).B(k_q).KII = [sqrt(abs(t(i_q,1,k_q,2))*Length_unit*G*Eeff),sqrt(abs(t(i_q,2*C(i_q).B(k_q).n_d,k_q,2))*Length_unit*G*Eeff)]
        C(i_q).B(k_q).K = [sqrt(C(i_q).B(k_q).KI(1)^2+C(i_q).B(k_q).KII(1)^2),sqrt(C(i_q).B(k_q).KI(2)^2+C(i_q).B(k_q).KII(2)^2)];

        KI_th = C(i_q).B(k_q).sig_e_rot(2,2)*G*sqrt(pi*C(i_q).B(k_q).a*Length_unit)
        KII_th = C(i_q).B(k_q).sig_e_rot(1,2)*G*sqrt(pi*C(i_q).B(k_q).a*Length_unit)

        t_max1 = max(t(i_q,:,k_q,1),[],"all");
        t_max2 = max(t(i_q,:,k_q,2),[],"all");

        if C(i_q).B(k_q).K(1) > K_IC && C(i_q).B(k_q).K(2) > K_IC
            
            delr1 = (C(i_q).B(k_q).d*t(i_q,1,k_q,1))/(n_steps*max(sig_amp,t_max1));
            delr2 = (C(i_q).B(k_q).d*t(i_q,1,k_q,2))/(n_steps*max(sig_amp,t_max2));
            C(i_q).B(k_q).D1(1).R_p = C(i_q).B(k_q).D1(1).R_p + 50*delr1*C(i_q).B(k_q).e_p;
            C(i_q).B(k_q).D2(1).R_p = C(i_q).B(k_q).D2(1).R_p + 50*delr1*C(i_q).B(k_q).e_p;
            delr1 = C(i_q).B(k_q).d*t(i_q,2*C(i_q).B(k_q).n_d,k_q,1)/(n_steps*max(sig_amp,t_max1));
            delr2 = C(i_q).B(k_q).d*t(i_q,2*C(i_q).B(k_q).n_d,k_q,2)/(n_steps*max(sig_amp,t_max2));
            C(i_q).B(k_q).D1(2*C(i_q).B(k_q).n_d).R_p = C(i_q).B(k_q).D1(2*C(i_q).B(k_q).n_d).R_p +50*delr1*C(i_q).B(k_q).e_p;
            C(i_q).B(k_q).D2(2*C(i_q).B(k_q).n_d).R_p = C(i_q).B(k_q).D2(2*C(i_q).B(k_q).n_d).R_p +50*delr2*C(i_q).B(k_q).e_p;
        end
end
        
                for j=2:2*C(i_q).B(k_q).n_d
                    if norm( C(i_q).B(k_q).D1(j).b )==0
                        x_cod_up(i_q,j,k_q,1)=x_cod_up(i_q,j-1,k_q,1);
                        y_cod_up(i_q,j,k_q,1)=y_cod_up(i_q,j-1,k_q,1);
                        x_cod_down(i_q,j,k_q,1)=x_cod_down(i_q,j-1,k_q,1);
                        y_cod_down(i_q,j,k_q,1)=y_cod_down(i_q,j-1,k_q,1);
                    end
                    if norm( C(i_q).B(k_q).D2(j).b )==0
                        x_cod_up(i_q,j,k_q,2)=x_cod_up(i_q,j-1,k_q,2);
                        y_cod_up(i_q,j,k_q,2)=y_cod_up(i_q,j-1,k_q,2);
                        x_cod_down(i_q,j,k_q,2)=x_cod_down(i_q,j-1,k_q,2);
                        y_cod_down(i_q,j,k_q,2)=y_cod_down(i_q,j-1,k_q,2);
                    end          
                end

   end
                
        end
end


    toc;
            calculate_COD();
            plot_COD();
            % plot_barriers();
%%%%%%%%%%%%%%%%%%%%%%%%%%Stress Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=generate_grid();
for i_p=1:n_c
 for k_p=1:n_branch   
    for j_p=1:2*C(i_p).B(k_p).n_d

        switch(C(i_p).B(k_p).mode)

            case 'pure'
        C(i_p).B(k_p).D(j_p).sigma=zeros(3,3);
      for II=1:n_grid
            for JJ=1:n_grid 
                R_pq(1)=C(i_p).B(k_p).D(j_p).R_p(1)-X(II,JJ);
                R_pq(2)=C(i_p).B(k_p).D(j_p).R_p(2)-Y(II,JJ);
                R_pq(3)=0;
                r_pq=sqrt(R_pq(1)^2+R_pq(2)^2);
                QQ=rotation(i_p,j_p,k_p,0);
                TT=rotation(i_p,j_p,k_p,1);
                 if r_pq==0
                            'self force=0';
                 else
                     r_vec=[X(II,JJ) Y(II,JJ) 0];
                     sigma_disk=disk_field_pressure(r_vec);
                    [stress,displacement]=dislocation_field(R_pq*QQ);
                    stress=norm(C(i_p).B(k_p).D(j_p).b)*TT*stress*TT';
                    displacement=norm(C(i_p).B(k_p).D(j_p).b)*displacement*QQ;
                   % S= C(i_p).B(k_p).D(j_p).sigma+stress+sigma_disk;
                   S= C(i_p).B(k_p).D(j_p).sigma+stress;
                   % for uniaxial
                   S_xx(II,JJ)=S(1,1)+sig_e(1,1);
                   S_yy(II,JJ)=S(2,2)+sig_e(2,2);
                   S_xy(II,JJ)=S(1,2)+sig_e(1,2);
                   % for disk
                   % S_xx(II,JJ)=S(1,1);
                   % S_yy(II,JJ)=S(2,2);
                   % S_xy(II,JJ)=S(1,2);
                   u_x=displacement(1);
                   u_y=displacement(2);
                 end
            end
      end
              Sigma_xx=Sigma_xx+S_xx;
              Sigma_yy=Sigma_yy+S_yy;
              Sigma_xy=Sigma_xy+S_xy;
              U=U+u_x;
              V=V+u_y;

            case 'mixed'


  C(i_p).B(k_p).D1(j_p).sigma=zeros(3,3);
  C(i_p).B(k_p).D2(j_p).sigma=zeros(3,3);

     for II=1:n_grid
            for JJ=1:n_grid 
                R_pq1(1)=C(i_p).B(k_p).D1(j_p).R_p(1)-X(II,JJ);
                R_pq1(2)=C(i_p).B(k_p).D1(j_p).R_p(2)-Y(II,JJ);
                R_pq1(3)=0;
                r_pq1=sqrt(R_pq1(1)^2+R_pq1(2)^2);

                R_pq2(1)=C(i_p).B(k_p).D2(j_p).R_p(1)-X(II,JJ);
                R_pq2(2)=C(i_p).B(k_p).D2(j_p).R_p(2)-Y(II,JJ);
                R_pq2(3)=0;
                r_pq2=sqrt(R_pq2(1)^2+R_pq2(2)^2);
                QQ=rotation(i_p,j_p,k_p,0);
                TT=rotation(i_p,j_p,k_p,1);

                 if r_pq1==0 && r_pq2==0
                            'self force=0';
                 else
                     r_vec=[X(II,JJ) Y(II,JJ) 0];
                     sigma_disk=disk_field_pressure(r_vec);
                    [stress1,displacement1]=dislocation_field(R_pq1*QQ);
                    stress1=norm(C(i_p).B(k_p).D1(j_p).b)*TT*stress1*TT';
                    displacement1=norm(C(i_p).B(k_p).D1(j_p).b)*displacement1*QQ;

                    [stress2,displacement2]=dislocation_field(R_pq2*QQ);
                    stress2=norm(C(i_p).B(k_p).D2(j_p).b)*TT*stress2*TT';
                    displacement2=norm(C(i_p).B(k_p).D2(j_p).b)*displacement2*QQ;
                   % S= C(i_p).B(k_p).D1(j_p).sigma+C(i_p).B(k_p).D2(j_p).sigma+stress1+stress2+sigma_disk;
                   S= C(i_p).B(k_p).D1(j_p).sigma+C(i_p).B(k_p).D2(j_p).sigma+stress1+stress2;
                   % for uniaxial
                   S_xx(II,JJ)=S(1,1)+sig_e(1,1);
                   S_yy(II,JJ)=S(2,2)+sig_e(2,2);
                   S_xy(II,JJ)=S(1,2)+sig_e(1,2);
                   % for disk
                     % S_xx(II,JJ)=S(1,1);
                     % S_yy(II,JJ)=S(2,2);
                     % S_xy(II,JJ)=S(1,2);
                     u_x=displacement1(1)+displacement2(1);
                     u_y=displacement1(2)+displacement2(2);
                 end
            end
      end
 
              Sigma_xx=Sigma_xx+S_xx;
              Sigma_yy=Sigma_yy+S_yy;
              Sigma_xy=Sigma_xy+S_xy;
              U=U+u_x;
              V=V+u_y;
        end
              
    end
  end   
end
 
plot_field(Sigma_xx,Sigma_yy,Sigma_xy,U,V)
