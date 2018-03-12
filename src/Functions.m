classdef Functions 
    properties
    end

    methods 
        function omega_fdot = update_omega( obj, data, omega, delta, h)
            for i = 1:data.nodes(3)
                if ismember(i,data.gen)
                    react = 0;
                    for j = 1:data.nodes(3)
                        if ~ismember(j,data.gen)
                            delta_ij = delta(i) - delta(j);
                            react = react + data.B_ij(i,j)*sin(delta_ij); 
                        end
                    end
                        omega_dot = -(data.D_i(i)/data.M_i(i))*omega(i)-(1/data.M_i(i))*(data.real_power(i) + react);
                        omega_fdot(i) = omega_dot;%*h + omega(i);
                else
                        omega_fdot(i) = 0;
                end
            end
        end
        function delta_fdot = update_delta(obj, data, omega, delta, neu, h)
            react  = 0;
            react1 = 0;
            for i = 2:data.nodes(3)
                if ismember(i,data.gen)
                    delta_dot = omega(i) - omega(1);
                    delta_fdot(i) = delta_dot;%*h + delta(i);
                else 
                    for j = 1:data.nodes(3)
                        if ismember(j,data.gen)
                            delta_ij = delta(i) - delta(j);
                            react = react + data.B_ij(i,j)*sin(delta_ij);
                        else
                            delta_ij = delta(i) - delta(j);
                            react1 = react1 + data.B_ij(i,j)*sin(delta_ij)*neu(i,j);
                        end
                    end
                    delta_dot = -(data.real_power(i) + react + react1) - omega(1);
                    delta_fdot(i) = delta_dot;%*h + delta(i);
                end         
            end
        end
        function neu_fdot = update_neu(obj, data, idx, a, b, delta, neu, h)        
            for i = 1:data.nodes(2)
                mi = data.network_data.branch(i,idx.FROM_BUS);
                mj = data.network_data.branch(i,idx.TO_BUS);
                delta_ij = delta(mi) - delta(mj);
                lamda = data.B_ij(mi,mj)*(1-cos(delta_ij))/data.W_ij(mi,mj);
                t = power_flow_f(a,b,neu(mi,mj));
                neu_dot = (t - lamda);
                neu_fdot(mi,mj) = neu_dot;%*h + neu(mi,mj);
                if isnan(neu(mi,mj)) || neu(mi,mj) == inf || neu(mi,mj) < 0
                    neu(mi,mj) = 0;
                end 
            end

        end
        
        
   end
end
    

