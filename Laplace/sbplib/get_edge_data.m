function edge_data = get_edge_data(edges_data,n1,n2)
edge_data_idx = find(all(edges_data(:,1:2) == [n1,n2],2));
if isempty(edge_data_idx)
    edge_data_idx = find(all(edges_data(:,1:2) == [n2,n1],2));
end
if isempty(edge_data_idx)
    error("can not find edge in edge data")
end
edge_data = edges_data(edge_data_idx,:);