#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h> 
#include <vector>
#include <time.h>
#include <mpi.h>

using namespace std;
const int imgX=3000; //  horizontal image resolution 
const int imgY=3000; //  vertical image resolution
const int iter_n=3000; // max. iteration number
const double yMin= -0.135;   // Mandelbrot scene y - range 
const double yMax= -0.115;
const double xMin= -0.79;  // Mandelbrot scene x -range
const double xMax= -0.77;
int img_array[imgX][imgY] = {0};  // our MAndelbrot set values array
int img_line[imgY] = {0};

// convergation function - base of the Mandelbrot set value generation
// it will get two parameters (x and y coordinates) and will give an iteration count in return
int converges(double cx,double cy)
{
    int n=0;
    double zx=0;
    double new_zx=0;
    double zy=0;
    // we iterate until max. iteration count iter_n, or until z^2 (complex!) runs over 4 - this means, our series will run to infinity, so it's not part of the set
    while( (n<iter_n) && (zx*zx + zy*zy)<4 )
    {
        // z * z => new_zx = (zx*zx - zy*zy)  new_zy = (zx*zy + zx*zy)   // we work with complex numbers
        // z*z + c = zx^2 - zy^2 +cx   +  i(zx*zy*2 + cy) 
        new_zx = zx*zx - zy*zy + cx;
        zy = 2*zx*zy + cy;
        zx = new_zx;
        n++;
    }
    return n;
}


int main(int argc, char **argv)
{
    //variables for MPI communication:
    int id, nproc;
    MPI_Status status;
    int answer[2];
    double question;

    double resX=0;  // Resolution of our iteration steps
    double resY=0; // this will be calculated by (yMax-yMin) / imgY later..
    double cx=xMin; // calculation will start at this point and we will change this dynamically
    double cy=yMin;
    double s_time, e_time; // we will show some timing data  

    // Initialize MPI:
    MPI_Init(&argc, &argv);
    // Get my rank:
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Get the total number of processors:
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    MPI_Barrier(MPI_COMM_WORLD); //for precize timing


    if(id == 0)
    { //Master
        if(argc<2)
        {
            cout << "Usage:" <<endl;
            cout << argv[0] <<" out_file.ppm" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1);
        }

        ofstream myfile;  // we will write to this file	
        string filename1(argv[1]);
        // filename1 = "mandel.ppm";
        char *fileName1 = (char*)filename1.c_str();

        //prepare the step resolution
        resX = (xMax-xMin) / imgX;
        resY = (yMax-yMin) / imgY;

        s_time=MPI_Wtime(); // we get a time value at this point

        // we do the calculation for every point of our complex plane, thus on our 2D image with appropriate steps
        for(int i=0;i<imgX;i++)
        {
            MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
            &status);
            // answer[0] -- from whom
            // answer[1] -- '-1' or the X coordinate

            if(answer[1]>=0){ //not the first answer
            MPI_Recv(&img_array[answer[1]][0], imgY, MPI_INT, answer[0],
            2, MPI_COMM_WORLD, &status);
            }
            MPI_Send(&i, 1, MPI_INT, answer[0], 3, MPI_COMM_WORLD);
            MPI_Send(&cx, 1, MPI_DOUBLE, answer[0], 4, MPI_COMM_WORLD);

            cx=cx+resX;
        }

        //the remaining answers:
        int term=-1;
        for(int i=1;i<nproc;++i)
        {
            MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&img_array[answer[1]][0], imgY, MPI_INT, answer[0], 2, MPI_COMM_WORLD, &status);

            //sending the termination signal
            MPI_Send(&term, 1, MPI_INT, answer[0], 3, MPI_COMM_WORLD);
        }

        // we get another time at this point, so we can calculate the elapsed time for the calculation
        e_time=MPI_Wtime();

        cout << "Time elapsed during calculation: " << e_time-s_time << " secs." << endl;;

        // file IO
        myfile.open(fileName1);
        myfile << "P3\r\n";
        myfile << imgX;
        myfile << " ";
        myfile << imgY;
        myfile << "\r\n";
        myfile << "255\r\n";

        // we have to colour our dataset. Actually, the members of the Mandelbrot set are used to be the same colour (black?) and have from point of visualisations view no interest.
        // the outer points are represented by assigning colours to their iteration steps and this generates vivid forms and colors
        for(int i=0;i<imgX;i++)
        {
            for(int j=0;j<imgY;j++)
            {

                if( (img_array[i][j] < 256) )  // we go from black to red in this range
                {
                    myfile << img_array[i][j] << " 0 0"; // (int)(84*pow(img_array[i][j],0.2)) << " 0 0"; //myfile << img_array[i][j] << " 0 0";
                }
                else if( img_array[i][j] < 512)  // we go from red to yellow in this range
                {
                    myfile << "255 " << img_array[i][j]-256 << " 0";
                }
                else if( img_array[i][j] < 768)  // we go from yellow to white in this range
                {
                    myfile << "255 255 " << img_array[i][j]-512;
                }
/*              
                // we could refine our palette for more resolution, more iteration-step images
                else if( img_array[i][j] < 1024)
                {
                myfile << 1024-img_array[i][j] << " 255 255"; 
                }
                else if( img_array[i][j] < 1280)
                {
                myfile << "0 "  << 1280-img_array[i][j] << " 255"; 
                }
                else if( img_array[i][j] < 1536)
                {
                myfile << "0 0 "  << 1536-img_array[i][j]; 
                }
*/
                else  // everything else is black
                {
                    myfile << "0 0 0 ";
                } 
                myfile << " ";
            }
            myfile << "\r\n";
        }

        myfile.close();  // we close our file

        e_time=MPI_Wtime(); // and give another elapsed time info (IO included)

        cout << "Time elapsed total: " << e_time-s_time << " secs \r\n";
    }
    
    else //Slave
    { 
        //prepare the step resolution
        resX = (xMax-xMin) / imgX;
        resY = (yMax-yMin) / imgY;

        int i;
        answer[0]=id;
        answer[1]=-1;
        MPI_Send(answer, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);

        while(1)
        {
            MPI_Recv(&i, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
            if(i<0) break; //got termination command!
            answer[1]=i;

            MPI_Recv(&question, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);

            cy=yMin;  // at every new step in X direction, we start at the first Y value 

            for(int j=0;j<imgY;j++)
            {
                img_line[j]=converges(question,cy);
                cy=cy+resY;
            }
            MPI_Send(answer, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(img_line, imgY, MPI_INT, 0, 2, MPI_COMM_WORLD);
        }
    }

    // Terminate MPI:
    MPI_Finalize();

    return 0;
}