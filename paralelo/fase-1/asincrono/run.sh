mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n $1 $2 $3
