% edges: 
% sel,sta,tip,ip,iq,np,name,r,x,b,g,ktr,kti,n_anc,bd,pl_ip,ql_ip,i_max
% nodes: 
% sel,sta,tip,ny,name,uhom,nsx,pn,qn,pg,qg,vzd,qmin,qmax,bshr,gshr,vras,delta,qnr,pnr
% polin:
% nsx,p0,p1,p2,frec,q0,q1,q2,frecq,p3,p4,q3,q4,umin

nodes = readtable('networks/ues_vostok_nodes.csv');
edges = readtable('networks/ues_vostok_edges.csv');
polin = readtable('networks/ues_vostok_polin.csv');

network = preprocess_network(nodes, edges);
[network, edges] = preprocess_data(network, edges, polin);

for i = 1:1%size(network, 2)
    disp(['Island ', num2str(i), ' (', ...
        num2str(height(network{1, i}.all)), ' nodes)']);
        
    [y, r] = admittance_matrix(network{1, i}.all, edges);
    network{1, i} = newton(0.1, 20, network{1, i}, y, r);
    
    disp('-----------------------------');
end
%------------------------------------------

function network = newton(precision, max_iters, network, y, r)
    norm_S = 0;

    sy = sparse(y);
    sr = sparse(r);

    for iter = 1:max_iters
        [P, Q] = balance(network, "all", sy, sr);

        new_norm_S = sqrt(sum(P)^2 + sum(Q)^2);
        print_iteration_info(iter, P, Q, network, y, norm_S / new_norm_S);
        norm_S = new_norm_S;
        
        S = [P; Q];
        
        if abs(S(:)) < precision
            disp("PF calculated successfully!");
            
%             [~, gQ] = balance(network, ...
%                 @(node)(0), ...
%                 @(node)(network.generators.contains(node)), ...
%                 y, r);
%             
%            network.all.qg(nodes.tip == "Ген") = -gQ;
            
            [bP, bQ] = balance(network, "base", y, r);
            
%             network.all.pg(network.base_node) = -bP;
%             network.all.qg(network.base_node) = -bQ;
            
            disp(['Base node: ', num2str(-bP), ...
                ' ', num2str(-bQ)]);
            
            %nodes = commit_slc_result(nodes, slc);
            
            return;
        end
        
        sJ = sparse_jacobi(network, sy, sr);
        %J = jacobi(network, y, r); 
        dX = -sJ\S;
        
        network = increment_guess(network, dX);       
    end
    
    disp("Max iterations count reached!");
end
%------------------------------------------

function network = increment_guess(network, dX)
    v_count = height(network.all) - sum(network.node_type(:) == 1) - 1;
    dV = dX(1:v_count);
    dDelta = dX(v_count+1:end);
    
    vi = 1;
    di = 1;
    for i = 1:height(network.all)
        if network.node_type(i) == 2
            network.v(i) = network.v(i) + dV(vi);
            network.delta(i) = network.delta(i) + dDelta(di);
            vi = vi + 1;
            di = di + 1;
        end
        
        if network.node_type(i) == 1
            network.delta(i) = network.delta(i) + dDelta(di);
            di = di + 1;
        end
    end
end
%------------------------------------------

function print_iteration_info(iter, P, Q, network, y, Rk)
    if iter == 1
        fprintf('%5s %20s %8s %4s %4s %4s %4s %5s %8s %4s\n', ...
            'Iter', 'Max mismatch (P, Q)', 'Nodes', ...
            '>V', 'Node', '<V', 'Node', ...
            'Angle', 'Edge', 'Rk');
    end
    
    [max_pmismatch, max_pmismatch_id] = max(abs(P));
    [max_qmismatch, max_qmismatch_id] = max(abs(Q));
    
    dV = network.v./network.rat_v;
    [max_dv, max_dv_id] = max(abs(dV));
    [min_dv, min_dv_id] = min(abs(dV));
    
    max_angle = 0;
    max_angle_i = 0;
    max_angle_j = 0;
    
%     for i = 1:size(y, 1)
%         for j = 1:size(y, 2)
%             if y(i, j) == 0
%                 continue;
%             end
%             
%             angle = abs(network.delta(i) - network.delta(j));
%             if angle >= max_angle
%                 max_angle = angle;
%                 max_angle_i = network.all.ny(i);
%                 max_angle_j = network.all.ny(j);
%             end
%         end 
%     end
    
    fprintf(...
    '%5d %9.2f %9.2f %4d %4d %4.2f %4d %4.2f %4d %5.2f %4d %4d %4.2f\n', ...
        iter, max_pmismatch, max_qmismatch, ...
        network.all.ny(max_pmismatch_id), network.all.ny(max_qmismatch_id), ...
        max_dv,  network.all.ny(max_dv_id), min_dv, network.all.ny(min_dv_id), ...
        max_angle, max_angle_i, max_angle_j, Rk);
end
%------------------------------------------

function nodes = commit_slc_result(nodes, slc)
    for i = 1:height(nodes)
        [pn, qn] = slc_load(nodes, i, slc, 1);
        nodes.pnr(i) = pn;
        nodes.qnr(i) = qn;
    end
end
%------------------------------------------

function [pn, qn] = slc_load(network, node, d)
    c = [1 1 0 0 0 0 1 0 0 0 0 1 1 0];
    if network.slc(node, 1) == 1
        c = network.slc(node, 2:end); 
    end
    
    if d == 1
        v = network.v(node)/network.rat_v(node);
        pn = real(network.s_load(node)) ...
            * (c(1) + c(2)*v + c(3)*(v^2));   
        qn = imag(network.s_load(node)) ...
            * (c(6) + c(7)*v + c(8)*(v^2)); 
    elseif d == 2
        v_nom = network.rat_v(node);
        pn = real(network.s_load(node)) ...
            * (c(2)/v_nom + (2*c(3)*network.v(node)/(v_nom)^2));   
        qn = imag(network.s_load(node)) ...
            * (c(7)/v_nom + (2*c(8)*network.v(node)/(v_nom)^2));
    end
end
%------------------------------------------

function [P, Q] = balance(network, nodes_filter, sy, sr)
    P = [];
    Q = [];

    [nzi, nzj, nzy] = find(sy);

    prev_i = 0;
    for k = 1:size(nzy)
        j = nzi(k);
        i = nzj(k);

        if nodes_filter == "all"
            calc_p = network.node_type(i) ~= 0;
            calc_q = network.node_type(i) == 2;
        elseif nodes_filter == "gen"
            calc_p = network.node_type(i) == 2;
            calc_q = network.node_type(i) == 2;
        elseif nodes_filter == "base"
            calc_p = network.node_type(i) == 0;
            calc_q = network.node_type(i) == 0;
        end

        if i ~= prev_i
            if calc_p
                P(end+1) = 0;
            end
            
            if calc_q
                Q(end+1) = 0;
            end
            prev_i = i;
        end

        vi = network.v(i);

        if i == j
            if network.slc(i, 1) == 1
                [pn, qn] = slc_load(network, i, 1);
            else
                pn = real(network.s_load(i));
                qn = imag(network.s_load(i));
            end

            if calc_p
                P(end) = P(end) - real(nzy(k)) * vi^2;    %#ok<*AGROW>
                P(end) = P(end) + real(network.s_gen(i));
                P(end) = P(end) - pn;
            end

            if calc_q
                Q(end) = Q(end) + imag(nzy(k)) * vi^2;     %#ok<*NASGU>
                Q(end) = Q(end) + imag(network.s_gen(i));
                Q(end) = Q(end) - qn;
            end
        else
            vj = network.v(j);
            di = network.delta(i);
            dj = network.delta(j);
            
            if calc_p
                P(end) = P(end) + ...
                    vi * vj ...
                    * abs(nzy(k)) * sin(di - dj - pi/2 + sr(i, j));
            end

            if calc_q
                Q(end) = Q(end) + ...
                    vi * vj ...
                    * abs(nzy(k)) * cos(di - dj + pi/2 + sr(i, j));
            end
        end
    end

    P = P.';
    Q = Q.';
end
%------------------------------------------

function J = sparse_jacobi(network, sy, sr)
    nodes_count = height(network.all);
    pq_count = nodes_count - sum(network.node_type(:) == 1) - 1;
    dPV = zeros(nodes_count - 1, pq_count);
    dPd = zeros(nodes_count - 1);
    dQV = zeros(pq_count);
    dQd = zeros(pq_count, nodes_count - 1);

    [nzi, nzj, nzy] = find(sy);

    prev_i = 0;
    shift_i = 0;
    for k = 1:size(nzy)
        j = nzi(k);
        i = nzj(k); 

        if i ~= prev_i
            if network.node_type(i) == 1
                shift_i = shift_i + 1;
            end

            prev_i = i;
        end

        if network.node_type(i) == 0
            continue;
        end

        igens_proceed = 0;
        jgens_proceed = 0;
        iflag = 0;
        jflag = 0;
        for l = 1:size(network.node_type)
            if (l >= i)
                iflag = 1;
            end
            if (l >= j)
                jflag = 1;
            end

            if iflag && jflag
                break;
            end

            if (network.node_type(l)) == 1
                if ~iflag
                    igens_proceed = igens_proceed + 1;
                end
                if ~jflag
                    jgens_proceed = jgens_proceed + 1;
                end
            end
        end

        if i == j
            [pn, qn] = slc_load(network, i, 2);
            
            if network.node_type(i) == 2
                dPV(i, i - igens_proceed) = dPV(i, i - igens_proceed) ...
                    - 2 * real(nzy(k)) * network.v(i) - pn;

                dQV(i - shift_i, i - igens_proceed) = ...
                    dQV(i - shift_i, i - igens_proceed) ...
                    + 2 * imag(nzy(k)) * network.v(i) - qn;
            end
        else
            di = network.delta(i);
            dj = network.delta(j);
            p_angle = di - dj - pi/2 + sr(i, j);
            q_angle = di - dj + pi/2 + sr(i, j);

            dPd(i, i) = dPd(i, i) + network.v(i) * network.v(j) ...
                * abs(nzy(k)) * cos(p_angle);

            if network.node_type(j) ~= 0
                dPd(i, j) = -network.v(i) * network.v(j) ...
                    * abs(nzy(k)) * cos(p_angle);
            end
            
            if network.node_type(i) == 2
                dQd(i - shift_i, i) = dQd(i - shift_i, i) ...
                    - network.v(i) * network.v(j) ...
                    * abs(nzy(k)) * sin(q_angle);

                if network.node_type(j) ~= 0
                    dQd(i - shift_i, j) = network.v(i) * network.v(j) ...
                        * abs(nzy(k)) * sin(q_angle);
                end
            end

            if network.node_type(i) == 2
                dPV(i, i - igens_proceed) = dPV(i, i - igens_proceed) ...
                    + network.v(j) * abs(nzy(k)) * sin(p_angle);
            end

            if network.node_type(j) == 2
                dPV(i, j - jgens_proceed) = network.v(i) ...
                    * abs(nzy(k)) * sin(p_angle);
            end

            if network.node_type(i) == 2
                
                dQV(i - shift_i, i - igens_proceed) = ...
                    dQV(i - shift_i, i - igens_proceed) + network.v(j) ...
                    * abs(nzy(k)) * cos(q_angle);

                if network.node_type(j) == 2
                    dQV(i - shift_i, j - jgens_proceed) = network.v(i) ...
                        * abs(nzy(k)) * cos(q_angle);
                end
            end

        end
    end

    J = sparse([dPV dPd; dQV dQd]);
end

% ip = vetv.ip;
% iq = vetv.iq;
% G = graph(ip, iq);
% plot(G,'Layout','force');
% 
% adjG = full(adjacency(G));
% incG = full(incidence(G));
% incG(balance_node, :) = []; 
%------------------------------------------

% density balance
%     for i = 1:height(network.all)
% 
%         [pn, qn] = slc_load(network, i, 1);
%         vi = network.v(i);
%         
%             if nodes_filter == "all"
%                 calc_p = network.node_type(i) ~= 0;
%                 calc_q = network.node_type(i) == 2;
%             elseif nodes_filter == "gen"
%                 calc_p = network.node_type(i) == 2;
%                 calc_q = network.node_type(i) == 2;
%             elseif nodes_filter == "base"
%                 calc_p = network.node_type(i) == 0;
%                 calc_q = network.node_type(i) == 0;
%             end
%         
%         if calc_p
%             P(end+1) = -real(y(i, i)) * vi^2;    %#ok<*AGROW>
%             P(end) = P(end) + real(network.s_gen(i)); 
%             P(end) = P(end) - pn;
%         end
% 
%         if calc_q
%             Q(end+1) = imag(y(i, i)) * vi^2;     %#ok<*NASGU>
%             Q(end) = Q(end) + imag(network.s_gen(i));
%             Q(end) = Q(end) - qn;
%         end
% 
%         for j = 1:height(network.all)
%             if j == i || y(i, j) == 0
%                 continue;
%             end
%             
%             vj = network.v(j);
%             di = network.delta(i);
%             dj = network.delta(j);
%             
%             if calc_p
%                 P(end) = P(end) + ...
%                     vi * vj ...
%                     * abs(y(i, j)) * sin(di - dj - pi/2 + r(i, j));
%             end
% 
%             if calc_q
%                 Q(end) = Q(end) + ...
%                     vi * vj ...
%                     * abs(y(i, j)) * cos(di - dj + pi/2 + r(i, j));
%             end
%         end
%     end


% function J = jacobi(network, y, r)
%     n = height(network.all);
%     dPV = zeros(n);
%     dPd = zeros(n);
%     dQV = zeros(n);
%     dQd = zeros(n);
%     
%     for i = 1:n      
%         [pn, qn] = slc_load(network, i, 2);
%         
%         dPV(i, i) = -2 * real(y(i, i)) * network.v(i) - pn;
%         dQV(i, i) = 2 * imag(y(i, i)) * network.v(i) - qn;
%         
%         for j = 1:n
%             if j == i || y(i, j) == 0
%                 continue;
%             end
%             
%             di = network.delta(i);
%             dj = network.delta(j);
%             p_angle = di - dj - pi/2 + r(i, j);
%             q_angle = di - dj + pi/2 + r(i, j);
% 
%             dPV(i, i) = dPV(i, i) + network.v(j) ...
%                 * abs(y(i, j)) * sin(p_angle);
%     
%             dPd(i, i) = dPd(i, i) + network.v(i) * network.v(j) ...
%                 * abs(y(i, j)) * cos(p_angle);
%                                
%             dQV(i, i) = dQV(i, i) + network.v(j) ...
%                 * abs(y(i, j)) * cos(q_angle);
%                 
%             dQd(i, i) = dQd(i, i) - network.v(i) * network.v(j) ...
%                 * abs(y(i, j)) * sin(q_angle);
%                 
%             dPV(i, j) = network.v(i) ...
%                 * abs(y(i, j)) * sin(p_angle);
% 
%             dPd(i, j) = -network.v(i) * network.v(j) ...
%                 * abs(y(i, j)) * cos(p_angle);
%                                       
%             dQV(i, j) = network.v(i) ...
%                 * abs(y(i, j)) * cos(q_angle);
%                         
%             dQd(i, j) = network.v(i) * network.v(j) ...
%                 * abs(y(i, j)) * sin(q_angle);
%                    
%         end
%     end
%     
%     J = escape_jacobi_params(network, dPV, dPd, dQV, dQd);
% end
% %------------------------------------------
% 
% function J = escape_jacobi_params(network, dPV, dPd, dQV, dQd)
%    
%     rows_erased = containers.Map(["dPV", "dPd", "dQV", "dQd"], [0, 0, 0, 0]);
%     cols_erased = containers.Map(["dPV", "dPd", "dQV", "dQd"], [0, 0, 0, 0]);
%     
%     for i = 1:height(network.all)
%         if network.node_type(i) == 0
%             dPV(i - rows_erased("dPV"), :) = [];
%             rows_erased("dPV") = rows_erased("dPV") + 1;
%             dPd(i - rows_erased("dPd"), :) = [];
%             rows_erased("dPd") = rows_erased("dPd")+ 1;
%             
%             dPV(:, i - cols_erased("dPV")) = [];
%             cols_erased("dPV") = cols_erased("dPV") + 1;
%             dPd(:, i - cols_erased("dPd")) = [];
%             cols_erased("dPd") = cols_erased("dPd") + 1;
%             
%             dQV(i - rows_erased("dQV"), :) = [];
%             rows_erased("dQV") = rows_erased("dQV") + 1;
%             dQd(i - rows_erased("dQd"), :) = [];
%             rows_erased("dQd") = rows_erased("dQd") + 1;
%             
%             dQV(:, i - cols_erased("dQV")) = [];
%             cols_erased("dQV") = cols_erased("dQV") + 1;
%             dQd(:, i - cols_erased("dQd")) = [];
%             cols_erased("dQd") = cols_erased("dQd") + 1;
%         end
%         
%         if network.node_type(i) == 1
%             dPV(:, i - cols_erased("dPV")) = [];
%             cols_erased("dPV") = cols_erased("dPV") + 1;
%             dQV(:, i - cols_erased("dQV")) = [];
%             cols_erased("dQV") = cols_erased("dQV") + 1;
%             
%             dQV(i - rows_erased("dQV"), :) = [];
%             rows_erased("dQV") = rows_erased("dQV") + 1;
%             dQd(i - rows_erased("dQd"), :) = [];
%             rows_erased("dQd") = rows_erased("dQd") + 1;
%         end        
%     end
% 
%     J = [dPV dPd; dQV dQd];
% end

