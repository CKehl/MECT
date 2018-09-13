The code contained in this package is described in 

Mean shift based clustering in high dimensions: 
A texture classification example. 
B. Georgescu, I. Shimshoni, P. Meer
in the proceedings of ICCV 2003.
Please refer to the paper for details

The package includes the following files:
fams.cpp
fams.h

and project files for MS Visual C++.


Program parameters:
K L k data_file_name input_dir [ [-j jump] | -p percent] | -h width [-f epsilon Kmin Kjump]

K,L are the LSH parameters. When they are both zero the LSH data
    structure is not built and the linear algorithm is run.

k   is the number of neighbors used in the construction of the pilot density.

data_file_name and input_dir describe the location of the input.
    Input files end with extension txt. So if for example the file is
    C:\temp\d45.txt . The data_file_name=d45 and the input_dir is C:\temp\

jump and percent are two ways to choose from which points to start the
    mean shift procedure. 
    The default is to perform the procedure on every point.
    jump means once every jump a point is selected. if jump is 5 then
    points 1,6,11,16... are chosen. 
    percent means that percent/100*n (n the number of points) are
    chosen at random with replacement. 
    Obviously, each run of the program will yield a different set.
    Only one of these options can be chosen.

width enables the user to run the fixed bandwidth mean shift
    procedure. The width*d (d the dimension of the data) 
    is the distance under the L1 norm used as the fixed bandwidth.

When -f option is used the program computes automatically the optimal K and L.
    the optimal K is searched between Kmin and K (the first parameter) with step Kjump
    the optimal L is searched between 1 and L (second parameter)
    epsilon represents the allowed error (in experiments epsilon = 0.05)
    

File formats:
Files are all ASCII files.
The input file has a one line header with two
number_of_points dimension
Then the input points are given as lines including d numbers.

The pilot process produces the neighborhood size for each point
If the input file is for example d040_1S.txt and the number of neighbors is k
the pilot file will be pilot_d040_1S_k.txt
If the file exists the program will assume that it was created by a
previous run and use the information 
in it and will not rerun the pilot creation procedure.
This enables the users to provide their own neighborhood sizes.

The output files are for example:
The result of MS on the selected points is in out_data_file_name.txt
The result of joined modes is in modes_data_file_name.txt

Example:
We are given a file d040_1S.txt
fams 30 46 200 d040_1S ./ -f 0.05 10 2
Finds K,L and optionally runs mean shift on all points.
To run mean shift on same data with K=24 and L=35
fams 24 35 200 d040_1S ./

NOTES:

- Internally the program scales the data between 0 and 2^16-1, therefore
  the precision is limited by range_of_data/2^16 

- Running mean shift on all the points is not necessary when only the modes
  are needed

- In the file fams.h there are several constants that influence the speed and
  accuracy of the program. For example you can set the bandwithd on which two 
  modes are joined or the number of trials on which the test is run when 
  finding K and L.