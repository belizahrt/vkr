function network = preprocess_network(nodes, edges)
    % g = make_graph_rastr(nodes, edges);
    
    network = split_islands(nodes, edges);
    
    i = 1;
    while i <= size(network, 2)
        island = network{i};
        
        correct = check_balance_node(island);
        
        if ~correct
            network{i} = {};
            disp('Balance node is not assigned! Island has been removed');
            continue;
        end
        
        % todo: implement zero impedance supernodes
        
        network{i} = assign_nodes_id(island);
        
        i = i + 1;
    end    
end

function g = make_graph_rastr(nodes, edges)
    g = cell(2, height(nodes));
    
    for i = 1:height(nodes)     
        g(1, i) = num2cell(nodes.ny(i));
        
        node_edges = edges(edges.ip == nodes.ny(i), :);
        for j = 1:height(node_edges)
            g{2, i}(j) = num2cell(node_edges.iq(j));
        end
    end
end

function g = bfs_rastr(edges, source, predicator)
    g = [];
    
    visited = containers.Map('KeyType','int32','ValueType','int32');
    queue = [];
    
    queue(end+1) = source;
    visited(source) = 1;
    
    while size(queue) ~= 0
        node = queue(end);
        queue(end) = [];
        
        g(end+1) = node;
        
        iq = edges(edges.ip == node, :).iq;
        ip = edges(edges.iq == node, :).ip;
        childs = cat(1, iq, ip);
        for i = 1:size(childs)
            child = childs(i);
              
            skip = predicator(node, child);
            if isempty(skip) 
                skip = predicator(child, node);
            end
                
            if ~skip
                %visited(child) = 1;
                continue;
            end
            
            if visited.isKey(double(child)) ~= 1
                queue(end+1) = child;
                visited(child) = 1;
            end
        end
    end
end

function islands = split_islands(nodes, edges)
    islands = {};

    while height(nodes) ~= 0
        
        g = bfs_rastr(edges, nodes.ny(1), ...
            @(ip, iq) ...
            (edges(edges.ip == ip & edges.iq == iq, :).sta == 0));
        
        island = table();
        for i = 1:size(g, 2)
            filter = nodes.ny == g(i);
            node = nodes(filter, :);
            
            island = [island; node];
            
            nodes(filter, :) = [];
        end
        
        if height(island) > 1
            islands{end+1} = island;
        end
    end
end

function correct = check_balance_node(island)
    correct = true;
    
    balance_node_count = 0;
    for i = 1:height(island)
        if island.tip(i) == "База"
            balance_node_count = balance_node_count + 1;
        end
    end
    
    if balance_node_count ~= 1
        correct = false;
    end
end

function [nodes, edges] = equate_zero_impedance_nodes(nodes, edges)

end

function nodes = get_zero_impedance_nodes(edges, source)
    
end

function island = assign_nodes_id(island)
    nodes_count = height(island);
    id = 1;
    for i = 1:nodes_count
        if island.tip(i) == "База"
            island.Id(i) = nodes_count;
        else
            island.Id(i) = id;
            id = id + 1;
        end
    end       
end