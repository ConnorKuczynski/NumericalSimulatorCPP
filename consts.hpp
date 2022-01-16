const double PI = 3.14159;
const double E = 2.71828;
const int DIM = 3; // simulator dimension
const int X = 0; //enumerate accessing values from pos/vel/acc arrays 
const int Y = 1;
const int Z = 2;

//THRUST MODES 

//add more as needed

const int HOVER = 0; //maintain a accerlation of zero
const int MAX_THRUST = 1; //maintain a constant maximium thrust in positive Z direction


//SIMULATION PARAMETERS

const double SECS_PER_ITR = .001; //1 ms for each iteration in flight profile
const double TIME_FINAL = 60*30; //(s) 1 minute