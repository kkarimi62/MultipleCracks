function update(i,k)

global C nu sig_e n_d_min



C(i).B(k).a = pdist([C(i).B(k).D(2*C(i).B(k).n_d).R_p;C(i).B(k).D(1).R_p],'euclidean');
if C(i).B(k).phi == 0
        C(i).B(k).b_mag=2*C(i).B(k).a*(1-nu)*sig_e(2,2)/n_d_min;
        else
            C(i).B(k).b_mag=2*C(i).B(k).a*(1-nu)*sig_e(1,1)/n_d_min;    
end

 for j=1:2*C(i).B(k).n_d
        C(i).B(k).D(j).b=-C(i).B(k).b_mag*C(i).B(k).b_direction;
 end

end