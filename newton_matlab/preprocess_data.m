function [network, edges] = preprocess_data(network, edges, slc)
    for isl = 1:size(network, 2)
        island.all = sortrows(network{isl}, 'Id');
        
        nodes_replace_nan_cols = ["gshr", "bshr", "pg", "pn", "qg", "qn", ...
            "qnr", "pnr"];
        edges_replace_nan_cols = ["r", "x", "g", "b", "ktr", "kti"];
        
        nodes_count = height(island.all);
        island.generators = java.util.HashSet;
        island.s_gen = zeros(nodes_count, 1);
        island.s_load = zeros(nodes_count, 1);
        island.v = zeros(nodes_count, 1);
        island.rat_v = zeros(nodes_count, 1);
        island.delta = zeros(nodes_count, 1);
        island.node_slc = java.util.HashMap;
        island.slc = java.util.HashMap;
        
        for node = 1:nodes_count
            %
            for col = nodes_replace_nan_cols
                if isnan(island.all.(col)(node))
                    island.all.(col)(node) = 0;
                end  
            end
            
            if isnan(island.all.vzd(node))
                island.all.vras(node) = island.all.uhom(node);
            else
                island.all.vras(node) = island.all.vzd(node);
            end
            
            island.all.delta(node) = 0;
            
            if island.all.tip(node) == "Ген"
                island.generators.add(node);

                island.all.qg(node) = 0;
            elseif island.all.tip(node) == "База"
                island.base_node = node;

                island.all.pg(node) = 0;
                island.all.qg(node) = 0;
            else
                island.all.vras(node) = island.all.uhom(node);
            end
            %
            
            island.v(node) = island.all.vras(node);
            island.rat_v(node) = island.all.uhom(node);
            island.delta(node) = island.all.delta(node);
            
            island.s_gen(node) = ...
                island.all.pg(node) + 1i*island.all.qg(node);
            island.s_load(node) = ...
                island.all.pn(node) + 1i*island.all.qn(node); 
            
            %
            node_slc = slc(slc.nsx == island.all.nsx(node), :);
            if island.all.nsx(node) == 1
                c = [0.83 -0.3 0.47; ...
                    3.7 -7 4.3];
                
                island.slc.put(island.all.nsx(node), c);
                island.node_slc.put(node, c);  
            elseif island.all.nsx(node) == 2
                c = [0.83 -0.3 0.47; ...
                    4.9 -10.1 6.2];
                
                island.slc.put(island.all.nsx(node), c);
                island.node_slc.put(node, c);                 
            end
            
            if height(node_slc) ~= 0
                c = [node_slc.p0 node_slc.p1 node_slc.p2; ...
                    node_slc.q0 node_slc.q1 node_slc.q2];
                
                island.slc.put(island.all.nsx(node), c);
                island.node_slc.put(node, c);
            end 
            
            %
        end
        
        network{isl} = island;
    end   
    
    
    %
    for edge = 1:height(edges)
        for col = edges_replace_nan_cols
            if isnan(edges.(col)(edge))
                edges.(col)(edge) = 0;
            end  
        end
        
        if edges.r(edge) == 0 && edges.x(edge) == 0
            edges.x(edge) = 0.2;
        end
    end
end