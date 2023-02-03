# parallel_programming_MPI

**Mandelbrot program**
There is a demo Mandlebrot that uses a cluster to compute an image of a Mandlebrot fractal. Compile and run the program with the following command:
mpiCC -o mandelbrot_parallel mandelbrot_parallel.cpp -lmpi
mpirun -np (num processors) ./mandelbrot_parallel out_file.ppm
You need to download the file "out_file.ppm" in your own local machine to see the output. Run the program with the numbers of processors as 2, 4, 8, 12, 16, and 24

**Parallel Sum**
parallel program to compute the sum of 1000 numbers read in from the file rand_data.txt by distributing the computation in blocks of 100 numbers over the processors.
compile and run the code with
mpicc -o parsum parsum.c -lmpi
mpirun -np 10 ./parsum
(You may want to capture the output of the run with mpirun -p all -np 10 . /parsum | tee log.txt)

**PI Calculation**
a simple program to determine the value of pi. The program uses the method to evaluate the integral of 4/(1+x*x) between 0 and 1.
you run the code with 5, 10, 20, 40, 60 or more processors, respectively, to see the speedups
