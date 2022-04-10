function [y, r] = admittance_matrix(nodes, edges)
    nodes_count = height(nodes);

    y = zeros(nodes_count, nodes_count);
    r = zeros(nodes_count, nodes_count);

    for k = 1:nodes_count - 1       
        y(k, k) = (nodes.gshr(k) - 1i*nodes.bshr(k)) * 10^(-6);
    end
        
    for k = 1:height(edges)
        if edges.sta(k) ~= 0
            continue;
        end
        
        i = nodes(nodes.ny == edges.ip(k), :).Id;
        j = nodes(nodes.ny == edges.iq(k), :).Id;
        
        % edge is not in island
        if isempty(i) || isempty(j)
            continue;
        end     
        
        z = edges.r(k) + 1i*edges.x(k);
               
        kt = edges.ktr(k) + 1i*edges.kti(k);
        kt_amp = abs(kt);
               
        y_edge = (edges.g(k) - 1i*edges.b(k)) * 10^(-6);
        
        if edges.tip(k) == "ЛЭП" || edges.tip(k) == "Выкл"
            y(i, j) = y(i, j) + (-1 / z);
            y(j, i) = y(i, j);

            y(i, i) = y(i, i) + (1 / z) + y_edge / 2;
            y(j, j) = y(j, j) + (1 / z) + y_edge / 2;
        elseif edges.tip(k) == "Тр-р"
            y(i, j) = y(i, j) + (-1 / (z * kt_amp));
            y(j, i) = y(i, j);

            y(i, i) = y(i, i) + ((1 / z)*(1 - 1/kt_amp) + y_edge);
            y(j, j) = y(j, j) + (1 / (z * kt_amp))*(1/kt_amp - 1);
                
            y(i, i) = y(i, i) + (1 / (z * kt_amp));
            y(j, j) = y(j, j) + (1 / (z * kt_amp));
        end
        
        r(i, j) = -angle(y(i, j)) + angle(1/kt);
        r(j, i) = r(i, j);
    end

%     y = sparse(y);
%     r = sparse(r);
end