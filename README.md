# CnCLULESH
# Author: Chenyang Liu

## Tiled cnc version of lulesh


usage : ./lulesh -n \<mesh_size> -b \<block_dim> [-p numprocs][-i max_iterations][-debug]

mesh_size is the number of elements along one dimension of the 3D (cube) mesh.

block_dim must divide evenly into mesh_size, is the size of each tiled block.



