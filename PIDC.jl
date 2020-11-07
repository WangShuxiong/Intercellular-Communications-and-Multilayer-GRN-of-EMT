using NetworkInference


println("Getting nodes...")
# input path of the data (gene*cell matrix)
data_path = joinpath(dirname(@__FILE__), "folder_name")
data_file_path = joinpath(data_path, "file_name")
nodes = get_nodes(data_file_path)

println("Inferring networks...")
#get the correlations of each genes
pidc_network = InferredNetwork(PIDCNetworkInference(), nodes)
#save the correlations, need to convert to gene-gene correlaion matrix later
write_network_file(joinpath(data_path, "correlations.txt"), pidc_network)
