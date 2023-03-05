## Here we shall add code notes.

### General Notes

    N denotes the particle number (effective).

    J denotes the cell number.

    L denotes the lenght of the integration region.

### 1D notes

The main variable is a vector its first 2N components are the N positions and the N momenta (or velocities in the non-relativistic version). After that the vector contains the J values of the electric field.


u = [x,v,E]=[$x_1,x_2$, ..., $x_N$, $v_1$, $v_2$ ..., $v_N$, $E_1$, ..., $E_J$]

length(u) = 2N+J

There are functions to change to pairs, that is:

u_r =[[$x_1,v_1$], [$x_2,v_2$], ..., [$x_N,v_N], E_1, ..., E_J]


### 2D or more notes 

Here the vector u is ordered as phase-space points and then the fields.
In 2D:

u_r =[[$x_1,y_1,px_1,px_1$], [$x_2,y_2,px_2,px_2$], ..., [$x_N,y_N,px_N,py_N$], [$Ex_{11},Ey_{11}, Ex_{12}, Ey_{12}..., Ex_{J,J},Ey_{JJ}$], [$B_{11}, B_{12}, ..., B_{JJ}$]

Recall that here $B = B_z$. 

Here length(u) = 4N + 3J[1]*J[2]

In 3D:

u_r =[[$x_1,y_1,z_1,px_1,px_1pz_1$], [$x_2,y_2,z_2,px_2,px_2,pz_2$], ..., [$x_N,y_N,z_N,px_N,py_N,pz_N$], [$Ex_{111},Ey_{111},Ez_{111} Ex_{112}, Ey_{112}, Ez_{112}..., Ex_{J[1],J[2],J[3]},Ey_{J[1],J[2],J[3]}$], $Bx_{111}, By_{111}, Bz_{111}, ..., Bz_{J[1],J[2],J[3]}$]

length(u) = 6N + 6J[1]*J[2]*J[3]

There are functions to retrieve each particle and each field. We then can view the fields as matrices of vectors. We try to do this with views.

For instance the chunck of data for the E vector in D=2 when J=(20,30) can be view as an array of sizes (2,20,30) or as a matrix of 2-vector of size (20,30)

E_array = reshape(u[range_E],(D,J...))


E_v = nestedview(E_array,1) is an array of vectors.








