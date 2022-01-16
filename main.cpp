
#include <iostream>
#include <cmath>
#include "main.hpp"
using namespace std;




/*
 * stores all the physical parameters of the rocket
 * list of param:
 * dry mass, fuel mass, oxidizer mass, length, width, height, (idea is to model each rocket part as a simple shape where we know the moment of inertia
 *       For example, we seperate the fuel tanks from the dry body, and model them as a sphere (or the best approximate shape), while body is a rectangle)
 * 
 * Summing the center of mass of each individual shape divided by total mass is equal to the center of mass of entire rocket 
 * 
 * we will model the initial rocket as a rectangular prism with a moment of inertia equal to 
 * 
 * Momentum of inertia will be stored relative to the P of the entire rocket
 */
class Environment {
public:
  const double GM = 3.986004418*pow(10,14); 
  const double R = 6371927; //radius of Earth in San Diego (meters)
  double* air;
  Environment(double* a) {
    this->air = a;
  }
  double* getAir()
  {
    return air;
  }
  double getG(double h)
  {
    return GM/pow((R + h),2);
  }
  /*
   * density (kg/m) of atmosphere as a function of height (meters)
   * https://www.grc.nasa.gov/www/k-12/rocket/atmosmet.html
   */
  double getP(int h)
  {
    if (h < 11000) {
      return (101.29*pow(((15.04-.00649*h + 273.1)/288.08),5.256))/(.2869*((15.04-.00649*h + 273.1)));
    }
    else if (h > 25000) {
      return 2.488*pow(((-131.21+.00299*h + 273.1)/216.6),-11.388)/(.2869*((-131.21+.00299*h + 273.1)));   
    }
    else {
      return 22.65*pow(E,(1.73-.000157*h))/(.2869*((-56.46 + 273.1)));
    }
  }
};
//all coordinate by default are measured relative to P of rocket
//to specify local coordinate system use suffix Local

//Shape 

float* Shape::getSA() {
  return sa; 
}
double* Shape::getI() {
  return I;
}
float* Shape::getOrigin() {
  return P;
}
double* Shape::updateCOVLocal() {
  for (int i = 0; i < DIM; i++) {
    covLocal[i] = cov[i] - P[i];
  } 
  return covLocal;
}
double* Shape::getCOVLocal() {
  return covLocal;
}
double* Shape::getCOV() {
  return cov;
}
double* Shape::updateI(float) { throw std::invalid_argument("Shape updateI() was run instead of a subclass shape"); }
void Shape::updateCOM(float*, float, float, float) { throw std::invalid_argument("Shape updateCOM() was run instead of a subclass shape"); }
void Shape::updateFilledVolume(Tank* t, double dltaM) { throw std::invalid_argument("Shape updateFilledVolume() was run instead of a subclass shape"); };
double* Shape::updateIFuel(Tank* t) {  throw std::invalid_argument("Shape updateIFuel() was run instead of a subclass shape"); }
  //  virtual double* Shape::updateI(float) { throw std::invalid_argument("Shape updateI() was run instead of a subclass shape");};
   // virtual void Shape::updateCOM(float*, float, float, float) { throw std::invalid_argument("Shape updateCOM() was run instead of a subclass shape"); };
   // virtual void Shape::updateFilledVolume(Tank* t, double dltaM);
   // virtual double* Shape::updateIFuel(Tank* t);
    

/*
 * length = X, width = Y, height = Z
 */

//END Shape

//RectangularPrism

int RectagularPrism::getAreaLW()
{
  return length*width;
}
int RectagularPrism::getLength()
{
  return length;
}
int RectagularPrism::getWidth()
{
  return width;
}
int RectagularPrism::getHeight()
{
  return height;
}
void RectagularPrism::updateFilledVolume(Tank* t, double dltaM) {
  t->updateFilledVolumeHelper(this, dltaM);
}
double* RectagularPrism::updateIFuel(Tank* t) {
  return t->updateIFuelHelper(this);
}
double* RectagularPrism::updateI(float mass)
{
  double* r = cov;

  I[X] = mass*((1/12.)*(pow(width,2) + pow(height,2)) + pow(r[Y],2) + pow(r[Z],2));
  I[Y] = mass*((1/12.)*(pow(length,2) + pow(height,2)) + pow(r[X],2) + pow(r[Z],2));
  I[Z] = mass*((1/12.)*(pow(length,2) + pow(width,2)) + pow(r[X],2) + pow(r[Y],2));

  return I;
}
void RectagularPrism::updateCOM(float* com, float fuelH, float fuelM, float dry_mass) 
{
  for (int i = 0; i < DIM; i++) {
    com[i] = 0;
    if (i == Z) { 
      com[Z] += ((fuelH/2 + P[Z])*fuelM)/(fuelM + dry_mass); 
      com[Z] += ((height/2. + P[Z])*dry_mass)/(fuelM + dry_mass);
    }
    else {
      com[i] += ((P[i])*fuelM)/(fuelM + dry_mass);
      com[i] += ((P[i])*dry_mass)/(fuelM + dry_mass);
    }
    
  }
  
}

//End RectangularPrism

//Sphere
float Sphere::getRadius() {
  return radius;
}
void Sphere::updateFilledVolume(Tank* t, double dltaM) {
  t->updateFilledVolumeHelper(this, dltaM);
}
float Sphere::getVolume() {
  return (4/3.)*PI*pow(radius,3);
}
double* Sphere::updateIFuel(Tank* t) {
  return t->updateIFuelHelper(this);
}
//hollow sphere
double* Sphere::updateI(float mass) {
  double* r = cov;
  I[X] = mass*((2/3.)*pow(radius,2) + pow(r[Y],2) + pow(r[Z],2));
  I[Y] = mass*((2/3.)*pow(radius,2) + pow(r[X],2) + pow(r[Z],2));
  I[Z] = mass*((2/3.)*pow(radius,2) + pow(r[Y],2) + pow(r[X],2));

  return I;
}
void Sphere::updateCOM(float* com, float fuelH, float fuelM, float dry_mass) 
{
  for (int i = 0; i < DIM; i++) {
    com[i] = 0;
    if (i == Z) { 
      float ang = 2*acos(1-fuelH/radius);
      com[Z] += ((radius - 4*radius*pow(sin(ang/2),3)/(3*(ang-sin(ang))) + P[Z])*fuelM)/(fuelM+dry_mass);
      com[Z] += ((radius/2. + P[Z])*dry_mass)/(fuelM + dry_mass);
    }
    else {
      com[i] += ((P[i] + radius)*dry_mass)/(fuelM + dry_mass);
      com[i] += ((P[i])*fuelH)/(fuelM + dry_mass);
    }
    
  }
  
}

//End Sphere

//Tank
/*
* currently we are assuming that the fuel stays perpendicular to the z axis of the rocket which obviously will need to changed in the future.
*/
//calculate momentum of inertia caused by fuel, calculate base on container of fuel
double* Tank::updateIFuel(Shape* s) {
  return (s->updateIFuel(this));
}
double* Tank::updateIFuelHelper(RectagularPrism* s)
{
  float Z_coor = com[Z];
  double* r = s->getCOV();
  int w = s->getWidth();
  int l = s->getLength();
  fuelI[X] = fuelM*((1/12.)*(pow(w,2) + pow(fuelH,2)) + pow(r[Y],2) + pow(Z_coor,2));
  fuelI[Y] = fuelM*((1/12.)*(pow(l,2) + pow(fuelH,2)) + pow(r[X],2) + pow(Z_coor,2));
  fuelI[Z] = fuelM*((1/12.)*(pow(w,2) + pow(l,2)) + pow(r[X],2) + pow(r[Y],2));

  return fuelI;
}
//http://blitiri.blogspot.com/2014/05/mass-moment-of-inertia-of-hemisphere.html
double* Tank::updateIFuelHelper(Sphere* s)
{
  int radius = s->getRadius();
  float* origin = s->getOrigin();
  float* r = com;
  float vol = s->getVolume();
  //calculate X and Y is same as finding filled sphere but find ratio of filled to full and then use parallel axis to account for center of mass not being center of volume
  //next shift to pivot point of the rocket. For Z direction solve integral for finding sphere but instead of -R to R calculate 0 to fuelH.
  //Account for difference of center of volume of sphere and center of mass and shift to axis of rocket pivot
  fuelI[X] = fuelM*((volumeFilled/vol)*(2/5.)*pow(radius,2) + pow(r[Y]-origin[Y],2) + pow(r[Z]-origin[Z],2) + (pow(r[Y],2) + pow(r[Z],2)));
  fuelI[Y] = fuelM*((volumeFilled/vol)*(2/5.)*pow(radius,2) + pow(r[X]-origin[X],2) + pow(r[Z]-origin[Z],2) + pow(r[X],2) + pow(r[Z],2));
  fuelI[Z] = fuelM*(PI/volumeFilled*(pow(fuelH,5)/5 - (2/3.)*pow(fuelH,3)*pow(radius,2) + fuelH*pow(radius,4)) + pow(r[Y]-origin[Y],2) + pow(r[X]-origin[X],2) + pow(r[Y],2) + pow(r[X],2));

  return fuelI;
}
//change in mass is in kg 
int Tank::updateMass(double dt)
{
  if (fuelM < 0) { return dry_mass; }
  double dltaM = this->mass_vel_out*dt;
  updateFilledVolume(s, dltaM);
  s->updateCOM(com, fuelH, fuelM, dry_mass);
  updateIFuel(s);
  fuelM -= dltaM;

  return (dry_mass+fuelM);
}
float Tank::getDryMass()
{
  return dry_mass;
}
float Tank::getFuelMass()
{
  return fuelM;
}
float Tank::getMass()
{
  return dry_mass + fuelM;
}
float* Tank::getCOM()
{
  s->updateCOM(com, fuelH, fuelM, dry_mass);
  return com;
}
Shape* Tank::getShape() {
  return s;
} 
void Tank::updateFilledVolume(Shape* s, double dltaM) {
  s->updateFilledVolume(this, dltaM);
}
void Tank::updateFilledVolumeHelper(RectagularPrism* s, double dltaM) {
  float dltaV = dltaM*densityFuel;
  float dltaH = dltaV/s->getAreaLW();
  volumeFilled -= dltaV;
  fuelH -= dltaH;
}
double* Tank::getIFuel() {
  return fuelI;
}
float Tank::getFuelHeight(){
  return fuelH;
}
void Tank::updateFilledVolumeHelper(Sphere* s, double dltaM) {
  float dltaV = dltaM*densityFuel;
  float dltaH = dltaV/(PI*fuelH*(2*s->getRadius()-fuelH));
  volumeFilled -= dltaV;
  fuelH -= dltaH;
}

    
    

    
    
//End Tank



double* Control::updateThrust(Dynamics::State* s, Rocket* r, float fuelMass, double dt) {
  //cout << "TIME ELAPSED: " << s->timeElapsed << " " << time_itr << endl;
  
  updateProfile(dt); 
  
  if (fuelMass < 0) {
    for (int i = 0; i < DIM; i++) {
      thrust[i] = 0;
    }
    return thrust;
  }
  //cout << time_itr << "/" << (int)(TIME_FINAL/SECS_PER_ITR) << endl;
  if (flight_profile[time_itr] == HOVER)
  {
    thrust[Z] += -s->acc[Z]*s->mass;
  }
  else if (flight_profile[time_itr] == MAX_THRUST) {
    thrust[Z] = r->getMaxThrust();
  }
  else {
    thrust[Z] = 0;
    cout << "BAD NEWS";
  }
  return thrust;
}
int* Control::initProfile() {

  flight_profile = new int[(int)(TIME_FINAL/SECS_PER_ITR)];
  for (int i = 0; i < TIME_FINAL/SECS_PER_ITR; i++) {
    flight_profile[i] = MAX_THRUST;
    /*
    if (i < (TIME_FINAL/SECS_PER_ITR) / 2) {
      flight_profile[i] = MAX_THRUST;
    }
    else {
      flight_profile[i] = HOVER;
    }
    */
  }
  return flight_profile;
}
void Control::updateProfile(double dt) {
  if (secs >= SECS_PER_ITR) 
  { 
    time_itr++;
    secs = dt; 
  }
  else { secs += dt; }
}
double* Control::getThrust() {
  return thrust;
}


//dry mass does not include dry mass of tanks
    /*
    * updatePARAM vs getPARAM (sometimes you just want to get the PARAM vs getting and updating it) this is done to avoid recalculating values that we already know
    *
    */
void Rocket::updateComDry()
{
  for (int i = 0; i < DIM; i++) {

    comDry[i] = cov[i]*dry_mass;
  }
}
Tank** Rocket::getTanks()
{
  return tanks;
}
double* Rocket::updateIDry()
{
  for (int j = 0; j < numShapes; j++)
  {
    //cout << "Shapes (X,Y,Z)" << endl;
    double* IShapes = s[j]->updateI(massParts[j]); 
    //cout << IShapes[X] << " " << IShapes[Y] << " " << IShapes[Z] << endl;
    //cout << "Dry mass tanks" << endl;
    double* ITanksDry = tanks[j]->getShape()->updateI(tanks[j]->getDryMass());
    //cout << ITanksDry[X] << " " << ITanksDry[Y] << " " << ITanksDry[Z] << endl;
    //cout << "Fuel height" << tanks[j]->getFuelHeight() << endl; 
    for (int i = 0; i < DIM; i++)
    {
      IDry[i] = IShapes[i];
      IDry[i] += ITanksDry[i];
    }
  }
  return IDry;
}
double* Rocket::updateI() {     
  for (int j = 0; j < numTanks; j++)
  {
    
    double* ITanks = tanks[j]->getIFuel();
    for (int i = 0; i < DIM; i++)
    {
      I[i] = IDry[i];
      I[i] += ITanks[i];  
    }
    //cout << I[X] << " " << I[Y] <<  " " << I[Z] << endl;
  }
  
  
  return I;
}
double* Rocket::updateCOM()
{
  
  for (int j = 0; j < DIM; j++)
  {
    com[j] = comDry[j];
  }
  
  for (int i = 0; i < numTanks; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      com[j] += (tanks[i]->getMass()*(tanks[i]->getCOM()[j]));
      com[j] /= tot_mass; 
    }    
  }
  return com;
}
double* Rocket::updateCOP()
{
  return cop;
}
double* Rocket::getI()
{
  return I;
}
int Rocket::updateMass(double dt)
{
  tot_mass = dry_mass;
  for(int i = 0; i < numTanks; i++)
  {
    Tank* tank = tanks[i];
    tot_mass += tank->updateMass(dt);
  }
  updateI();
  return tot_mass;
}
/*
* return Coeffient of drag based on the geometry of the surface of each direction of the rocket
*/
float* Rocket::getCD() {
  return cd;
}
double* Rocket::getCOM() {
  return com;
}
double* Rocket::getCOP() {
  return cop;
}
float* Rocket::getSA() {
  return sa;
}
int Rocket::getMass(){
  return tot_mass;
}
const int* Rocket::getP() {
  return P;
}
const int Rocket::getNumParts() {
  return numShapes;
}
const int Rocket::getNumTanks() {
  return numTanks;
}
const double Rocket::getMaxThrust() {
  return max_thrust;
}


 /* 
void RK4(double[] initialVec, double[] vec, double (*diffEq)(float, double, double), float dt) {
    int[] n = int[DIM];
 
    float[] k1, k2, k3, k4, k5;
 
    for (int i = 0; i < DIM; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1[i] = //dt*diffEq(0, initialVec);
        k2[i] = dt*diffEq(.5*dt, initialVec + 0.5*k1[i], );
        k3[i] = dt*diffEq(.5*dt, initialVec + 0.5*k2[i]);
        k4[i] = dt*diffEq(dt, initialVec + k3[i]);
 
        vec[i] += (1.0/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
} */


/*
 * ORDER of calculating values
 * ALWAYS update F before M 
 * 
 *  updateFdrag
 *  updateFg
 *  updateFThrust
 * 
 *  updateNetF
 * 
 *  updateMdrag
 *  updateMg
 *  updateMThrust
 * 
 *  updateNetM
 *  updateAcc
 * 
 *  updateNetM
 *  updateAngAcc
 *  
 * 
 *  clearNet
 */



void Dynamics::State::updatePrev(State* curr)
{
  ang = curr->ang; //global frame variables 
  pos = curr->pos;
  vel = curr->vel;
  timeElapsed = curr->timeElapsed;
  ang_vel = curr->ang_vel;
  ang_acc = curr->ang_acc;
  acc = curr->acc;
  mass = curr->mass;
}
string Dynamics::State::to_string() 
{
  string output = "(X,Y,Z) \n";
  string posStr = "POSITION: (";
  string velStr = "VELOCITY: (";
  string accStr = "ACCERLATION: (";
  string angStr = "ANGLE: (";
  string ang_velStr = "ANGULAR VELOCITY: (";
  string ang_accStr = "ANGULAR ACCERLATION: (";
  for (int i = 0; i < DIM; i++) {
    if (i < DIM - 1) { 
      posStr += (std::to_string(pos[i]) + ","); 
      velStr += (std::to_string(vel[i]) + ","); 
      accStr += (std::to_string(acc[i]) + ","); 
      angStr += (std::to_string(ang[i]) + ","); 
      ang_velStr += (std::to_string(ang_vel[i]) + ","); 
      ang_accStr += (std::to_string(ang_acc[i]) + ","); 
    }
    else { 
      posStr += (std::to_string(pos[i]) + ")\n"); 
      velStr += (std::to_string(vel[i]) + ")\n"); 
      accStr += (std::to_string(acc[i]) + ")\n"); 
      angStr += (std::to_string(ang[i]) + ")\n"); 
      ang_velStr += (std::to_string(ang_vel[i]) + ")\n"); 
      ang_accStr += (std::to_string(ang_acc[i]) + ")\n"); 
    }
  } 
  output = posStr + velStr + accStr + angStr + ang_velStr + ang_accStr;
  return output;
}
  string Dynamics::to_string()
  {
    string output = s->to_string();
    string Istr = "MOMENT OF INERTIA: (";
    string COMstr = "CENTER OF MASS: (";
    string massStr = "Mass: (" + std::to_string(r->getMass()) + ")\n";
    for (int i = 0; i < DIM; i++) {
      if (i < DIM - 1) { 
        Istr += (std::to_string(r->getI()[i]) + ","); 
        COMstr += (std::to_string(r->getCOM()[i]) + ","); 
      }
      else { 
        Istr += (std::to_string(r->getI()[i]) + ")\n"); 
        COMstr += (std::to_string(r->getCOM()[i]) + ")\n"); 
      }
    } 
    output += (Istr + COMstr + massStr);
    return output;
  }
  void Dynamics::RK4(double dt)
  {
    int n[DIM];
  
    double k1Pos[DIM], k2Pos[DIM], k3Pos[DIM], k4Pos[DIM];
    double k1Ang[DIM], k2Ang[DIM], k3Ang[DIM], k4Ang[DIM];  

//initial force calculations
    Dynamics* state0 = new Dynamics(*this); //deep clone current state
    double* forces[3] = {mDrag,mThrust,mG};

    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
   // state0->updateFThrust();
    
    state0->updateNetF(forces);
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();

    for (int i = 0; i < DIM; i++) {
      k1Pos[i] = dt*dPosdt(0*dt, state0->getVel()[i], state0->getAcc()[i]); //dltaPos
      k1Ang[i] = dt*dPosdt(0*dt, state0->getAngVel()[i], state0->getAngAcc()[i]); //dltaAng
    }
    

    //this order is important we update position based off .5*dltaPos
    //then we update vel because its already applied from dPos/dt
    state0->updatePos(.5, k1Pos);
    state0->updateVel(.5*dt);
    state0->updateAng(.5, k1Ang);
    state0->updateAngVel(.5*dt);

    state0->r->updateMass(.5*dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF(forces);
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();
    state0->clearNet();

    for (int i = 0; i < DIM; i++) {
      k2Pos[i] = dt*dPosdt(.5*dt, getVel()[i], getAcc()[i]);
      k2Ang[i] = dt*dPosdt(.5*dt, getAngVel()[i], getAngAcc()[i]);
    }

    state0->updatePos(.5, k2Pos);
    state0->updateVel(.5*dt);
    state0->updateAng(.5, k2Ang);
    state0->updateAngVel(.5*dt);

    state0->r->updateMass(.5*dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF(forces);
    state0->updateAcc();
    
    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();
    state0->clearNet();

    for (int i = 0; i < DIM; i++) {
      k3Pos[i] = dt*dPosdt(.5*dt, getVel()[i], getAcc()[i]);
      k3Ang[i] = dt*dPosdt(.5*dt, getAngVel()[i], getAngAcc()[i]);
    }

    state0->updatePos(dt, k3Pos);
    state0->updateVel(dt);
    state0->updateAng(dt, k3Ang);
    state0->updateAngVel(dt);

    state0->r->updateMass(dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF(forces);
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();
    state0->clearNet();

    for (int i = 0; i < DIM; i++) {
      k4Pos[i] = dt*dPosdt(dt, getVel()[i], getAcc()[i]);
      k4Ang[i] = dt*dPosdt(dt, getAngVel()[i], getAngAcc()[i]);

      this->s->pos[i] += (1.0/6.0)*(k1Pos[i] + 2*k2Pos[i] + 2*k3Pos[i] + k4Pos[i]);
      this->s->ang[i] += (1.0/6.0)*(k1Ang[i] + 2*k2Ang[i] + 2*k3Ang[i] + k4Ang[i]);
    }

    
  }
  void Dynamics::updateState(double dt) {
    double* forces[3] = {fDrag,fThrust,fG};
    double* moments[3] = {mDrag,mThrust,mG};
    

    this->updateFdrag();
    this->updateFg();
    this->updateFThrust(s->prev, dt);
    this->updateNetF(forces);
    this->updateAcc();

    //this->updateMdrag();
    //this->updateMg();
    //this->updateMThrust();
    //this->updateInerMom();
    //this->updateNetM(moments);
    //this->updateAngAcc();
    s->timeElapsed += dt;
  }
  void Dynamics::euler(double dt) {
    this->updateState(dt);

    for (int i = 0; i < DIM; i++)
    {
      s->pos[i] += dPosdt(dt, s->vel[i], s->acc[i]);
      s->vel[i] += s->acc[i]*dt;
    }
    for (int i = 0; i < DIM; i++)
    {
      s->ang[i] += dPosdt(dt, s->ang_vel[i], s->ang_acc[i]);
      s->ang_vel[i] += s->ang_acc[i]*dt;
    }
    s->mass = this->r->updateMass(dt);
    this->r->updateCOM();

    s->prev->updatePrev(s);
    clearNet();
  }

  double Dynamics::dPosdt(double dt, double vel, double acc)
  {
    return .5*acc*pow(dt,2) + vel*dt;
  }
  /*
   * this may get replaced/moved in the future but this succesfully calculates the cross_p, assumes DIM = 3
   */ 
  void Dynamics::crossProduct(double dist[], double force[], double moment[])
  {
    moment[0] = dist[1] * force[2] - dist[2] * force[1];
    moment[1] = dist[2] * force[0] - dist[0] * force[2];
    moment[2] = dist[0] * force[1] - dist[1] * force[0];
  }

  /*
   * @param double[] air, model air resistance as a vector 
   * vel is the velocity of the rocket, need to find relative wind vector
   * cd is coffiecient of drag of the rocket in each vector direction / it is proportional to surface area of that part of the rocket
   * TODO account for air resistance that flow around the rocket
   */
  void Dynamics::updateFdrag()
  {
    double P = e->getP(s->pos[Z]); 
    for (int i = 0; i < DIM; i++)
    {
      if (-s->vel[i] + e->air[i] > 0) {
        fDrag[i] = .5*pow((-s->vel[i] + e->air[i]),2)*P*r->getCD()[i]*r->getSA()[i];
      }
      else {
        fDrag[i] = -.5*pow((-s->vel[i] + e->air[i]),2)*P*r->getCD()[i]*r->getSA()[i];
      }
      //cout << "Environment air DIM: " << i << " " << fDrag[i] << endl;
    } 
  }
  /*
   * the moment created from the force of drag
   */
  void Dynamics::updateMdrag()
  {
    double dist[DIM];
    for (int i = 0; i < DIM; i++) {
      dist[i] = r->getCOP()[i];
    }
    crossProduct(dist, fDrag, mDrag);
  }



  void Dynamics::updateFg()
  {
    fG[Z] = -e->getG(s->pos[Z])*r->getMass(); 
  }
  void Dynamics::updateMg()
  {
    crossProduct(r->getCOM(), fG, mG);
  }

   /*
   * given by the control block
   */
  void Dynamics::updateFThrust(Dynamics::State* s, double dt)
  {
    float fuelMass = 0;
    for (int i = 0; i < r->getNumTanks(); i++)
    {
      fuelMass += r->getTanks()[i]->getFuelMass();
    }
    double* newThrust = 0;

    newThrust = c->updateThrust(s, r, fuelMass, dt);

    
    for (int i = 0; i < DIM; i++)
    {
        fThrust[i] = newThrust[i];
    }
  }
  void Dynamics::updateMThrust()
  {
    crossProduct(r->getCOM(), fThrust, mThrust);
  }

  //moment caused by movement
  void Dynamics::updateInerMom()
  {
    crossProduct(r->getCOM(), netF, inerM);
  }

  void Dynamics::updateNetF(double** forces)
  {
    for (int i = 0; i < NUM_FORCES; i++) {
      for (int j = 0; j < DIM; j++) {
        netF[j] += forces[i][j]; 
      }
    }

  }
  void Dynamics::updateNetM(double** moments)
  {
    for (int i = 0; i < NUM_FORCES; i++) {
      for (int j = 0; j < DIM; j++) {
        netM[j] += moments[i][j]; 
      }
    }
  } 

  void Dynamics::clearNet()
  {
    for (int i = 0; i < DIM; i++){
      netF[i] = 0;
      netM[i] = 0;
    }
  }
  void Dynamics::updateAcc()
  {
    for (int i = 0; i < DIM; i++){
      s->acc[i] = netF[i]/r->getMass();
    }
  }
  void Dynamics::updateVel(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->vel[i] += dt*s->acc[i];
    }
  }
  double* Dynamics::getVel() {
    return s->vel;
  }
  double* Dynamics::getAcc() {
    return s->acc;
  }
  double* Dynamics::getAngVel() {
    return s->ang_vel;
  }
  double* Dynamics::getAngAcc() {
    return s->ang_acc;
  }
  void Dynamics::updatePos(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->pos[i] += dt*s->vel[i];
    }
  }

  void Dynamics::updateVel(double dt, double acc[]) {
    for (int i = 0; i < DIM; i++)
    {
      s->vel[i] += dt*acc[i];
    }
  }
  void Dynamics::updatePos(double dt, double vel[]) {
    for (int i = 0; i < DIM; i++)
    {
      s->pos[i] += dt*vel[i];
    }
  }

  void Dynamics::updateAngAcc()
  {
    for (int i = 0; i < DIM; i++){
      s->ang_acc[i] = (netM[i]- inerM[i])/(r->getI()[i]);
    }
  }
  void Dynamics::updateAngVel(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->ang_vel[i] += dt*s->ang_acc[i];
    }
  }
  void Dynamics::updateAng(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->ang[i] += dt*s->ang_vel[i];
    }
  }
  void Dynamics::updateAng(double dt, double ang_vel[]) {
    for (int i = 0; i < DIM; i++)
    {
      s->ang[i] += dt*ang_vel[i];
    }
  }

Dynamics::State* initState(Rocket* r)
{
  
  double* ang = new double[DIM] {0,0,0}; //global frame
  double* pos = new double[DIM] {0,0,0}; 
  double* vel = new double[DIM] {0,0,0};
  double* ang_vel = new double[DIM] {0,0,0};
  double* acc = new double[DIM] {0,0,0};
  double* ang_acc = new double[DIM] {0,0,0};

  Dynamics::State* s = new Dynamics::State(ang,pos,vel,ang_vel,acc,ang_acc,r->getMass());
  return s;
}
Shape** initStructure()
{
  int length = 10; //meters X
  int height = 10; //meters Y
  int width = 10; //meters Z
  float* P = new float[DIM] {0, 0, 0}; //origin of shape is same as origin of rocket
  double* COV = new double[DIM] {0., 0, height/2.};
  float* surfaceArea = new float[DIM] {(float)(height*width), (float)(length*width), (float)(length*height)};

  Shape** shapes = new Shape*[1] {new RectagularPrism(surfaceArea,COV,P,length*width*height,length,width,height)};
  return shapes;
}
Tank** initTanks() 
{
  int length = 10; //meters X
  int height = 10; //meters Y
  int width = 10; //meters Z
  float* P = new float[DIM] {0, 0, 0}; //origin of tank is same as origin of rocket
  double* COV = new double[DIM] {0, 0, height/2.};
  float* surfaceArea = new float[DIM] {(float)(height*width), (float)(length*width), (float)(length*height)};

  Shape* tankShape = new RectagularPrism(surfaceArea,COV,P,length*width*height,length,width,height);
  float dry_mass = 1000; //kg
  float massF = 10; //kg
  double mass_vel_out = .1; //kg/s
  double fuelD = 0.657; //kg/m^3
  float volFill = 1; //m^3
  float fuelHeight = .5; //meters
  
  Tank** tanks = new Tank*[1] {new Tank(tankShape, dry_mass, massF, mass_vel_out, fuelD, volFill, fuelHeight)}; 
  return tanks;
}
Rocket* initRocket()
{
  Shape** shapes = initStructure();
  int* massParts = new int[1] {1000};
  double* cov = (shapes[0])->getCOV();
  float* coeffDrag = new float[DIM] {1.28, 1.28, 1.28}; //https://www.grc.nasa.gov/www/k-12/airplane/shaped.html
  int d_m = 1000; //kg
  float* surfaceArea = shapes[0]->getSA();
  float m_vel_out = .1; //m/s set as a parameters before, not actually used for anything currently, mass changes are handled in the tanks
  int numTanks = 1;
  int numShape = 1;
  float max_t = 50000; //N
  Tank** t = initTanks();

  Rocket* r = new Rocket(shapes, cov, coeffDrag, d_m, surfaceArea, m_vel_out, t, numTanks, numShape, massParts, max_t); 
  return r;
}


Environment* initEnvironment()
{
  return (new Environment(new double[DIM] {0,0,0}));
}
Control* initControl()
{
  return (new Control(new double[DIM] {0,0,0}));
}
void sampleData(Dynamics* sys, double timeElapsed)
{
  cout << sys->to_string() << endl;
}
int main() {
  double dt = .0001; //.0001 seconds per physics tick
  double timeElapsed = 0;
  
  double secsPerSample = dt*loopPerSample;
  double secs = secsPerSample;
  Rocket* r = initRocket();
  Dynamics::State* s = initState(r);
  Dynamics* sys = new Dynamics(r, initEnvironment(), initControl(), s);

  while (timeElapsed < TIME_FINAL) {
    sys->euler(dt);

    timeElapsed += dt;
    if (secs >= secsPerSample) 
    {
      sampleData(sys, timeElapsed); 
      secs = dt; 
    }
    else { secs += dt; }
  }
  cout << "DONE";
}
