//We might want to seperate this into 3 files (rocketPARAMS, simulationsPARAMS, etc)

const double PI = 3.14159;
const double E = 2.71828;
const int DIM = 3; // simulator dimension
const int X = 0; //enumerate accessing values from pos/vel/acc arrays 
const int Y = 1;
const int Z = 2;

//ROCKET PARAMETERS //want to migrate parameters from main.cpp to consts hpp so its easier to change

//In theory we could change this so that we can precisely control the thrust vector at each specified time instant


const int HOVER = 0; //maintain a accerlation of zero
const float MAX_THRUST = 50000; //maintain a constant maximium thrust in positive Z direction


//SIMULATION PARAMETERS

double dt = .00001; //.00001 seconds per physics tick, may be changed to a variable time step if needed
const int loopPerSample = 10000; //sample rate for simulation 
const double SECS_PER_ITR = .001; //1 ms for each iteration in flight profile
const double TIME_FINAL = 60*100; //(s) 30 minute



/*
 * x = yaw (positive x is out of the page)
 * y = pitch (positive y is to the right of the page)
 * z = roll (positive z is up to the page)
 *
 * https://en.wikipedia.org/wiki/Rotation_matrix Use general form to convert from global axis to rocket axis 
 * https://math.stackexchange.com/questions/2895880/inversion-of-rotation-matrix Use for inverse rotation matrix
 * 
 * TODO (1) Transform Fg to rocket frame to solve for the moment created by Fg, use global Fg to solve for global acc
 * // did it need to verify it works
 * (2) Transform Fdrag to rocket axis: do this by transforming global velocity of air and rocket to rocket frame (need to align with surface area to solve)
 * next transform Fdrag back to global axis by using the transpose of the transformation matrix to solve for global acc.
 * Use rocket frame Fdrag to solve for moment created by Fdrag
 * 
 * (3) Use transposed transformation matrix on thrust to go from rocket frame to global frame, use this to calculate global acc. Use rocket frame thrust to calculate
 * moment created by thrust.
 */